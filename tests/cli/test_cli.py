"""CLI の必須引数・ヘルプ表示と、主要オプションの出力を確認する。"""

import csv
import subprocess
import sys
from pathlib import Path


SCRIPT_PATH = Path(__file__).parents[2] / "app" / "Phylo-MIP.py"


def run_cli(*arguments):
    bootstrap = (
        "import runpy, sys, types; "
        "ete3 = types.ModuleType('ete3'); ete3.NodeStyle = object; "
        "sys.modules['ete3'] = ete3; "
        "sys.argv = sys.argv[1:]; "
        "runpy.run_path(sys.argv[0], run_name='__main__')"
    )
    return subprocess.run(
        [sys.executable, "-c", bootstrap, str(SCRIPT_PATH), *arguments],
        capture_output=True,
        text=True,
        check=False,
    )


def test_cli_requires_input_csv():
    # 入力CSVを必須引数として検証する / Verify that the input CSV is required.
    result = run_cli()

    assert result.returncode == 2
    assert "the following arguments are required: input_csv" in result.stderr


def test_cli_help_accepts_primary_options():
    # 主要CLIオプションがヘルプに表示されることを確認する / Verify that primary CLI options are advertised in help.
    result = run_cli("--help")

    assert result.returncode == 0
    assert "--top" in result.stdout
    assert "--onlyp" in result.stdout
    assert "--tree" in result.stdout
    assert "--method {NJ,ML}" in result.stdout


def test_onlyp_writes_expected_fasta_and_csv(tmp_path):
    # --onlyp のFASTA/CSV出力とclassフィルタを固定する / Pin the --onlyp FASTA/CSV output and class filtering behavior.
    input_path = tmp_path / "input.csv"
    with input_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["qseqid", "accessionID", "class", "pident", "qseq"])
        writer.writerow(["q1", "ACC001", "Insecta", "99.5", "ATGC"])
        writer.writerow(["q2", "ACC002", "Mammalia", "98.0", "GGTA"])

    result = run_cli(str(input_path), "--onlyp", "--class", "Insecta", "--o", "regression")

    assert result.returncode == 0, result.stderr
    output_dirs = list(tmp_path.glob("phylomip_output_*/taxonomy"))
    assert len(output_dirs) == 1
    taxonomy_dir = output_dirs[0]
    assert (taxonomy_dir / "regression_filtered.csv").exists()
    fasta = (taxonomy_dir / "regression.fasta").read_text(encoding="utf-8")
    assert ">q1_ACC001_Insecta" in fasta
    assert "ATGC" in fasta
    assert ">q2_ACC002_Mammalia" in fasta
    assert "GGTA" in fasta
