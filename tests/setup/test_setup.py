"""setup.shのlauncher管理と移行警告を確認する / Verify launcher ownership and migration warnings."""

import os
import stat
import subprocess
from pathlib import Path


REPOSITORY_ROOT = Path(__file__).parents[2]
SETUP_SCRIPT = REPOSITORY_ROOT / "setup.sh"


def run_setup(tmp_path, *arguments, legacy=False):
    # Docker buildをstub化してsetupの副作用だけを確認する / Stub Docker to verify setup side effects.
    home = tmp_path / ("legacy-home" if legacy else "new-home")
    home.mkdir(exist_ok=True)
    if legacy:
        legacy_bin = home / "bin"
        legacy_bin.mkdir()
        (legacy_bin / "phylo-mip").write_text("legacy phylo-mip\n", encoding="utf-8")
        (legacy_bin / "merge_data").write_text("legacy merge_data\n", encoding="utf-8")

    fake_bin = tmp_path / "fake-bin"
    fake_bin.mkdir(exist_ok=True)
    docker_log = tmp_path / "docker.log"
    docker = fake_bin / "docker"
    docker.write_text(
        "#!/usr/bin/env bash\n"
        f"printf '%s\\n' \"$*\" >> \"{docker_log}\"\n",
        encoding="utf-8",
    )
    docker.chmod(docker.stat().st_mode | stat.S_IXUSR)

    environment = os.environ.copy()
    environment.update({"HOME": str(home), "PATH": f"{fake_bin}:/usr/bin:/bin"})
    return subprocess.run(
        ["/bin/bash", str(SETUP_SCRIPT), *arguments],
        cwd=REPOSITORY_ROOT,
        env=environment,
        capture_output=True,
        text=True,
        check=False,
    ), home, docker_log


def test_setup_does_not_create_home_launcher_or_edit_shell_rc(tmp_path):
    # 新規環境でlauncherとshell設定を作成しないことを確認する / Verify no launcher or shell rc is created in a new environment.
    result, home, docker_log = run_setup(tmp_path)
    assert result.returncode == 0, result.stderr
    assert not (home / "bin").exists()
    assert not (home / ".bashrc").exists()
    assert not (home / ".bash_profile").exists()
    assert not (home / ".zshrc").exists()
    assert docker_log.read_text(encoding="utf-8").count("build") == 1

    repeated, repeated_home, repeated_log = run_setup(tmp_path)
    assert repeated.returncode == 0, repeated.stderr
    assert not (repeated_home / "bin").exists()
    assert repeated_log.read_text(encoding="utf-8").count("build") == 2


def test_setup_warns_and_preserves_legacy_launchers(tmp_path):
    # legacy launcherを警告のみで保持することを確認する / Verify legacy launchers are warned about and preserved.
    result, home, _ = run_setup(tmp_path, "--update", legacy=True)
    assert result.returncode == 0, result.stderr
    assert "legacy launcher detected" in result.stdout
    assert str(home / "bin" / "phylo-mip") in result.stdout
    assert str(home / "bin" / "merge_data") in result.stdout
    assert (home / "bin" / "phylo-mip").read_text(encoding="utf-8") == "legacy phylo-mip\n"
    assert (home / "bin" / "merge_data").read_text(encoding="utf-8") == "legacy merge_data\n"
    assert not list((home / "bin").glob("*.backup.*"))
    assert not (home / ".bashrc").exists()


def test_repository_launchers_are_canonical_entrypoints(tmp_path):
    # repository直下launcherの存在・実行・usage表示を確認する / Verify canonical launcher presence, execution, and usage.
    environment = os.environ.copy()
    environment["HOME"] = str(tmp_path / "home")
    (tmp_path / "home").mkdir()
    for launcher, expected_text in (("phylo-mip", "Usage:"), ("merge_data", "Insufficient arguments")):
        path = REPOSITORY_ROOT / launcher
        assert path.is_file()
        assert os.access(path, os.X_OK)
        result = subprocess.run(
            [str(path)], cwd=REPOSITORY_ROOT, env=environment,
            capture_output=True, text=True, check=False,
        )
        assert result.returncode == 1
        assert expected_text in result.stdout


def test_repository_launchers_forward_valid_commands_to_docker(tmp_path):
    # 有効な引数をDockerへ渡せることを確認する / Verify valid arguments are forwarded to Docker.
    home = tmp_path / "home"
    fake_bin = tmp_path / "bin"
    home.mkdir()
    fake_bin.mkdir()
    docker_log = tmp_path / "docker.log"
    docker = fake_bin / "docker"
    docker.write_text(
        "#!/usr/bin/env bash\n"
        f"printf '%s\\n' \"$*\" >> \"{docker_log}\"\n",
        encoding="utf-8",
    )
    docker.chmod(docker.stat().st_mode | stat.S_IXUSR)

    input_csv = tmp_path / "input.csv"
    input_csv.write_text("qseqid,sallacc,pident,qseq\nq1,ACC1,99.0,ATGC\n", encoding="utf-8")
    qiime = tmp_path / "qiime.tsv"
    qiime.write_text("Sample\t#OTU ID\tAbundance\nS1\tq1\t1\n", encoding="utf-8")
    phylo = tmp_path / "phylo.csv"
    phylo.write_text("qseqid,accessionID\nq1,ACC1\n", encoding="utf-8")

    environment = os.environ.copy()
    environment.update({"HOME": str(home), "PATH": f"{fake_bin}:/usr/bin:/bin"})
    phylo_result = subprocess.run(
        [str(REPOSITORY_ROOT / "phylo-mip"), str(input_csv), "--onlyp"],
        cwd=REPOSITORY_ROOT, env=environment, capture_output=True, text=True, check=False,
    )
    merge_result = subprocess.run(
        [str(REPOSITORY_ROOT / "merge_data"), "-q", str(qiime), "-p", str(phylo), "-f", "csv", "-o", "merged.csv"],
        cwd=REPOSITORY_ROOT, env=environment, capture_output=True, text=True, check=False,
    )

    assert phylo_result.returncode == 0, phylo_result.stderr
    assert merge_result.returncode == 0, merge_result.stderr
    docker_calls = docker_log.read_text(encoding="utf-8").splitlines()
    assert len(docker_calls) == 2
    assert all("phylo-mip" in call for call in docker_calls)
