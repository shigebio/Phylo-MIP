import csv
import re
import subprocess
import sys
from pathlib import Path


SCRIPT_PATH = Path(__file__).parents[2] / "app" / "Phylo-MIP.py"


def run_onlyp(input_path, *arguments):
    bootstrap = (
        "import runpy, sys, types; "
        "ete3 = types.ModuleType('ete3'); ete3.NodeStyle = object; "
        "sys.modules['ete3'] = ete3; "
        "sys.argv = sys.argv[1:]; "
        "runpy.run_path(sys.argv[0], run_name='__main__')"
    )
    return subprocess.run(
        [sys.executable, "-c", bootstrap, str(SCRIPT_PATH), str(input_path), "--onlyp", *arguments],
        capture_output=True,
        text=True,
        check=False,
    )


def write_input(path):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["qseqid", "accessionID", "class", "pident", "qseq"])
        writer.writerow(["q1", "ACC001", "Insecta", "99.5", "ATGC"])


def test_output_directory_structure_and_explicit_basename(tmp_path):
    # timestamp付きoutput directoryと--o basenameを確認する / Verify timestamped output and explicit basename.
    input_path = tmp_path / "input.csv"
    write_input(input_path)

    result = run_onlyp(input_path, "--o", "named")

    assert result.returncode == 0, result.stderr
    output_dirs = list(tmp_path.glob("phylomip_output_*/"))
    assert len(output_dirs) == 1
    assert re.fullmatch(r"phylomip_output_\d{8}_\d{6}", output_dirs[0].name)
    assert {path.name for path in output_dirs[0].iterdir()} == {
        "alignment", "bptp", "mptp", "phylogeny", "taxonomy"
    }
    assert (output_dirs[0] / "taxonomy" / "named.fasta").exists()


def test_output_directory_default_basename_uses_current_input_name(tmp_path):
    # --o未指定時の現行basenameを固定する / Pin the current basename behavior without --o.
    input_path = tmp_path / "sample.csv"
    write_input(input_path)

    result = run_onlyp(input_path)

    assert result.returncode == 0, result.stderr
    taxonomy_dirs = list(tmp_path.glob("phylomip_output_*/taxonomy"))
    assert len(taxonomy_dirs) == 1
    assert (taxonomy_dirs[0] / "sample.csv_output.fasta").exists()