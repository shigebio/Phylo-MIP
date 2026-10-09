"""setup.shのlauncher管理と移行警告を確認する / Verify launcher ownership and migration warnings."""

import os
import stat
import subprocess
from pathlib import Path


REPOSITORY_ROOT = Path(__file__).parents[2]
SETUP_SCRIPT = REPOSITORY_ROOT / "setup.sh"


def run_setup(tmp_path, *arguments, legacy=False, docker_available=True, docker_info_status=0, docker_build_status=0):
    # Dockerのinfo/buildをstub化してsetupの判定を確認する / Stub Docker info/build to verify setup decisions.
    home = tmp_path / ("legacy-home" if legacy else "new-home")
    home.mkdir(exist_ok=True)
    if legacy:
        legacy_bin = home / "bin"
        legacy_bin.mkdir()
        (legacy_bin / "phylo-mip").write_text("legacy phylo-mip\n", encoding="utf-8")
        (legacy_bin / "merge_data").write_text("legacy merge_data\n", encoding="utf-8")
        (home / ".bashrc").write_text("export PATH=\"$HOME/bin:$PATH\"\n# user setting\n", encoding="utf-8")

    fake_bin = tmp_path / "fake-bin"
    fake_bin.mkdir(exist_ok=True)
    docker_log = tmp_path / "docker.log"
    if docker_available:
        docker = fake_bin / "docker"
        docker.write_text(
            "#!/usr/bin/env bash\n"
            "case \"$1\" in\n"
            f"  info) printf '%s\\n' info >> \"{docker_log}\"; exit {docker_info_status} ;;\n"
            f"  build) printf '%s\\n' build >> \"{docker_log}\"; exit {docker_build_status} ;;\n"
            "  *) exit 0 ;;\n"
            "esac\n",
            encoding="utf-8",
        )
        docker.chmod(docker.stat().st_mode | stat.S_IXUSR)
    else:
        # DockerをPATHから除外しつつdirnameだけ提供する / Hide Docker while keeping dirname available.
        (fake_bin / "dirname").symlink_to("/usr/bin/dirname")

    environment = os.environ.copy()
    path = f"{fake_bin}:/usr/bin:/bin" if docker_available else str(fake_bin)
    environment.update({"HOME": str(home), "PATH": path})
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
    assert docker_log.read_text(encoding="utf-8").splitlines() == ["info", "build"]

    repeated, repeated_home, repeated_log = run_setup(tmp_path)
    assert repeated.returncode == 0, repeated.stderr
    assert not (repeated_home / "bin").exists()
    assert repeated_log.read_text(encoding="utf-8").splitlines() == ["info", "build", "info", "build"]


def test_setup_warns_and_preserves_legacy_launchers(tmp_path):
    # legacy launcherを警告のみで保持することを確認する / Verify legacy launchers are warned about and preserved.
    result, home, _ = run_setup(tmp_path, legacy=True)
    assert result.returncode == 0, result.stderr
    assert "legacy launcher detected" in result.stdout
    assert str(home / "bin" / "phylo-mip") in result.stdout
    assert str(home / "bin" / "merge_data") in result.stdout
    assert (home / "bin" / "phylo-mip").read_text(encoding="utf-8") == "legacy phylo-mip\n"
    assert (home / "bin" / "merge_data").read_text(encoding="utf-8") == "legacy merge_data\n"
    assert not list((home / "bin").glob("*.backup.*"))
    assert (home / ".bashrc").read_text(encoding="utf-8") == "export PATH=\"$HOME/bin:$PATH\"\n# user setting\n"


def test_setup_fails_when_docker_cli_is_missing(tmp_path):
    # Docker CLI未導入時にbuildせず失敗することを確認する / Verify setup fails before build when Docker CLI is missing.
    result, _, docker_log = run_setup(tmp_path, docker_available=False)
    assert result.returncode != 0
    assert "Docker CLI is not available" in result.stderr
    assert "Setup complete." not in result.stdout
    assert not docker_log.exists()


def test_setup_fails_when_docker_engine_is_unavailable(tmp_path):
    # Docker engine接続失敗時にbuildせず失敗することを確認する / Verify setup fails before build when the engine is unavailable.
    result, _, docker_log = run_setup(tmp_path, docker_info_status=1)
    assert result.returncode != 0
    assert "Docker engine is not accessible" in result.stderr
    assert "docker info" in result.stderr
    assert docker_log.read_text(encoding="utf-8").splitlines() == ["info"]
    assert "Setup complete." not in result.stdout


def test_setup_fails_when_docker_build_fails(tmp_path):
    # Docker build失敗時に非0終了・成功表示なしとなることを確認する / Verify build failure returns non-zero without success output.
    result, _, docker_log = run_setup(tmp_path, docker_build_status=7)
    assert result.returncode == 7
    assert "Docker image build failed" in result.stderr
    assert docker_log.read_text(encoding="utf-8").splitlines() == ["info", "build"]
    assert "Setup complete." not in result.stdout


def test_setup_succeeds_only_after_docker_info_and_build(tmp_path):
    # Docker info/build成功後だけ成功表示することを確認する / Verify success is reported only after info and build succeed.
    result, _, docker_log = run_setup(tmp_path)
    assert result.returncode == 0
    assert docker_log.read_text(encoding="utf-8").splitlines() == ["info", "build"]
    assert "Setup complete." in result.stdout


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
