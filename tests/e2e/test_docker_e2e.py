"""Docker wrapper 実行後に、ホスト側へ出力ディレクトリが戻ることを確認する。"""

import os
import pty
import shutil
from pathlib import Path

import pytest


@pytest.mark.skipif(
    os.environ.get("RUN_DOCKER_E2E") != "1" or shutil.which("docker") is None,
    reason="Docker E2E requires RUN_DOCKER_E2E=1 and a Docker executable",
)
@pytest.mark.e2e
def test_docker_e2e_writes_output_on_host(tmp_path):
    # Docker外側の成果物を確認する / Verify that Docker writes artifacts on the host.
    project_root = Path(__file__).parents[2]
    input_path = tmp_path / "docker_input.csv"
    input_path.write_text(
        "qseqid,sallacc,pident,qseq\nq1,ACC001,99.5,ATGC\n",
        encoding="utf-8",
    )

    # The project wrapper intentionally uses ``docker run -it``.  Give it a
    # pseudo-terminal so the E2E test exercises the wrapper rather than
    # failing before Docker starts with "the input device is not a TTY".
    if os.name == "nt":
        pytest.skip("Docker wrapper E2E requires a POSIX pseudo-terminal")
    status = pty.spawn([str(project_root / "phylo-mip"), str(input_path), "--onlyp"])
    assert os.waitstatus_to_exitcode(status) == 0
    output_dirs = list(tmp_path.glob("phylomip_output_*"))
    assert output_dirs
    assert (output_dirs[0] / "taxonomy").exists()
    assert any((output_dirs[0] / "taxonomy").iterdir())
