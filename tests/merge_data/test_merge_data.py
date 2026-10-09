"""TSV/CSV の区切り判定と、OTU ID を基準にした merge 結果を確認する。"""

import importlib.util
from pathlib import Path


MODULE_PATH = Path(__file__).parents[2] / "app" / "merge_data.py"


def load_merge_data_module():
    spec = importlib.util.spec_from_file_location("merge_data", MODULE_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_detect_delimiter_and_debug_read(tmp_path, capsys):
    # TSV区切り判定とdebug読込を確認する / Verify TSV delimiter detection and debug reading.
    module = load_merge_data_module()
    tsv_path = tmp_path / "data.tsv"
    tsv_path.write_text("qseqid\tvalue\nq1\tA\n", encoding="utf-8")

    assert module.detect_delimiter(tsv_path) == "\t"
    headers, rows = module.read_file_into_dict(tsv_path, "qseqid", debug=True)

    assert headers == ["qseqid", "value"]
    assert rows == {"q1": ["q1", "A"]}
    assert "Read 1 rows" in capsys.readouterr().out


def test_merge_files_matches_otu_ids_and_preserves_unmatched_rows(tmp_path):
    # qseqid照合とunmatched行の保持を確認する / Verify qseqid matching and preservation of unmatched rows.
    module = load_merge_data_module()
    fixture_dir = Path(__file__).parents[1] / "fixtures"
    qiime_path = tmp_path / "qiime.tsv"
    phylo_path = tmp_path / "phylo.csv"
    qiime_path.write_bytes((fixture_dir / "merge_data" / "qiime.tsv").read_bytes())
    phylo_path.write_bytes((fixture_dir / "merge_data" / "phylo.csv").read_bytes())

    module.merge_files(qiime_path, phylo_path, "merged.csv", "csv")

    output_path = tmp_path / "merged.csv"
    rows = output_path.read_text(encoding="utf-8").splitlines()
    assert rows[0] == "Sample,#OTU ID,qseqid,accessionID,class,order,family,taxonomic_name,pident,qseq,source,Abundance"
    assert rows[1] == "sample-a,q1,q1,ACC001,Insecta,Lepidoptera,FamilyA,Species one,99.5,ATGC,NCBI,10"
    assert rows[2] == "sample-b,missing,,,,,,,,,,4"


def test_merge_files_debug_reports_matching_ids_and_writes_output(tmp_path, capsys):
    # debug情報のID集計と、debug modeでのmerge出力生成を確認する。
    module = load_merge_data_module()
    fixture_dir = Path(__file__).parents[1] / "fixtures"
    qiime_path = tmp_path / "qiime.tsv"
    phylo_path = tmp_path / "phylo.csv"
    qiime_path.write_bytes((fixture_dir / "merge_data" / "qiime.tsv").read_bytes())
    phylo_path.write_bytes((fixture_dir / "merge_data" / "phylo.csv").read_bytes())

    module.merge_files(qiime_path, phylo_path, "debug_merged.csv", "csv", debug=True)

    output = capsys.readouterr().out
    assert "QIIME file has 2 OTU IDs" in output
    assert "Phylo-MIP file has 2 qseqids" in output
    assert "Number of common IDs: 1" in output
    assert "Percentage of QIIME IDs matched: 50.00%" in output
    assert (tmp_path / "debug_merged.csv").read_text(encoding="utf-8").splitlines()[1].startswith(
        "sample-a,q1,q1,ACC001"
    )
