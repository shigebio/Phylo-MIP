import importlib.util
from pathlib import Path


MODULE_PATH = Path(__file__).parents[1] / "app" / "merge_data.py"


def load_merge_data_module():
    spec = importlib.util.spec_from_file_location("merge_data", MODULE_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_detect_delimiter_and_debug_read(tmp_path, capsys):
    module = load_merge_data_module()
    tsv_path = tmp_path / "data.tsv"
    tsv_path.write_text("qseqid\tvalue\nq1\tA\n", encoding="utf-8")

    assert module.detect_delimiter(tsv_path) == "\t"
    headers, rows = module.read_file_into_dict(tsv_path, "qseqid", debug=True)

    assert headers == ["qseqid", "value"]
    assert rows == {"q1": ["q1", "A"]}
    assert "Read 1 rows" in capsys.readouterr().out


def test_merge_files_matches_otu_ids_and_preserves_unmatched_rows(tmp_path):
    module = load_merge_data_module()
    fixture_dir = Path(__file__).parent / "fixtures"
    qiime_path = tmp_path / "qiime.tsv"
    phylo_path = tmp_path / "phylo.csv"
    qiime_path.write_bytes((fixture_dir / "qiime.tsv").read_bytes())
    phylo_path.write_bytes((fixture_dir / "phylo.csv").read_bytes())

    module.merge_files(qiime_path, phylo_path, "merged.csv", "csv")

    output_path = tmp_path / "merged.csv"
    rows = output_path.read_text(encoding="utf-8").splitlines()
    assert rows[0] == "Sample,#OTU ID,qseqid,accessionID,class,order,family,taxonomic_name,pident,qseq,source,Abundance"
    assert rows[1] == "sample-a,q1,q1,ACC001,Insecta,Lepidoptera,FamilyA,Species one,99.5,ATGC,NCBI,10"
    assert rows[2] == "sample-b,missing,,,,,,,,,,4"