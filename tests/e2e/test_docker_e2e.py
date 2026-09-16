import os
import shutil
import subprocess
from pathlib import Path

import pytest


@pytest.mark.skipif(
    os.environ.get("RUN_DOCKER_E2E") != "1" or shutil.which("docker") is None,
    reason="Docker E2E requires RUN_DOCKER_E2E=1 and a Docker executable",
)
def test_docker_e2e_writes_output_on_host(tmp_path):
    # Docker外側の成果物を確認する / Verify that Docker writes artifacts on the host.
    project_root = Path(__file__).parents[2]
    input_path = tmp_path / "docker_input.csv"
    input_path.write_text(
        "qseqid,sallacc,pident,qseq\nq1,ACC001,99.5,ATGC\n",
        encoding="utf-8",
    )

    result = subprocess.run(
        [str(project_root / "phylo-mip"), str(input_path), "--onlyp"],
        cwd=project_root,
        capture_output=True,
        text=True,
        check=False,
        timeout=300,
    )

    assert result.returncode == 0, result.stderr
    output_dirs = list(tmp_path.glob("phylomip_output_*"))
    assert output_dirs
    assert (output_dirs[0] / "taxonomy").exists()
    assert any((output_dirs[0] / "taxonomy").iterdir())