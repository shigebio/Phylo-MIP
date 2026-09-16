"""taxonomy API の fallback/retry、pident 境界、CSV/FASTA metadata 整合を確認する。"""

import csv
import io
import json

import pandas as pd
import pytest
from Bio import SeqIO

from conftest import FIXTURES, load_phylo_module, write_blast_csv
from unittest.mock import Mock


def ncbi_record(organism="Mock species", taxonomy="cellular organisms; Eukaryota; Animalia; Arthropoda; Insecta; Ordera; Familya; Genusa"):
    return [{"GBSeq_organism": organism, "GBSeq_taxonomy": taxonomy}]


@pytest.mark.parametrize(
    ("pident", "expected"),
    [
        (98.00, "Mock species"),
        (97.99, "Genus GBIF"),
        (95.00, "Genus GBIF"),
        (94.99, "Family GBIF"),
        (92.00, "Family GBIF"),
        (91.99, "Order GBIF"),
        (85.00, "Order GBIF"),
        (84.99, "Low_Identity_Match"),
    ],
)
def test_pident_threshold_boundaries_are_stable(tmp_path, monkeypatch, pident, expected):
    # pident境界値と現行rank選択を確認する / Verify pident boundaries and current rank selection.
    input_path = tmp_path / "input.csv"
    write_blast_csv(input_path, [{"qseqid": "q1", "sallacc": "ACC1", "pident": pident, "qseq": "ATGC"}])
    module = load_phylo_module(input_path, monkeypatch, "--onlyp")
    module["fetch_ncbi_data"] = lambda accession: ncbi_record()
    module["get_gbif_taxonomic_info"] = lambda species: {
        "source": "GBIF", "species": "Mock species", "genus": "Genus GBIF",
        "family": "Family GBIF", "order": "Order GBIF", "class": "Insecta",
    }

    _, csv_entry = module["process_row"](0, module["df"].iloc[0])

    assert csv_entry[5] == expected


def test_gbif_success_and_fasta_csv_metadata_alignment(tmp_path, monkeypatch):
    # GBIF成功時のtaxonomyとCSV/FASTA metadata対応を確認する / Verify GBIF taxonomy and CSV/FASTA metadata alignment.
    input_path = tmp_path / "input.csv"
    write_blast_csv(input_path, [{"qseqid": "q1", "sallacc": "ACC1", "pident": 99.5, "qseq": "AT-GC"}])
    module = load_phylo_module(input_path, monkeypatch, "--onlyp")
    module["fetch_ncbi_data"] = lambda accession: ncbi_record()
    module["get_gbif_taxonomic_info"] = lambda species: {
        "source": "GBIF", "species": "Species GBIF", "genus": "Genus GBIF",
        "family": "Family GBIF", "order": "Order GBIF", "class": "Insecta",
    }

    fasta_entry, csv_entry = module["process_row"](0, module["df"].iloc[0])
    csv_path = tmp_path / "taxonomy.csv"
    fasta_path = tmp_path / "taxonomy.fasta"
    module["save_csv"](csv_path, [csv_entry])
    module["save_fasta"](fasta_path, [fasta_entry])

    parsed_csv = pd.read_csv(csv_path)
    fasta_record = next(SeqIO.parse(fasta_path, "fasta"))
    assert {"qseqid", "accessionID", "taxonomic_name", "qseq", "source"} <= set(parsed_csv.columns)
    assert len(parsed_csv) == 1
    assert parsed_csv.loc[0, "qseqid"] == "q1"
    assert parsed_csv.loc[0, "accessionID"] == "ACC1"
    assert parsed_csv.loc[0, "taxonomic_name"] == "Species GBIF"
    assert str(fasta_record.seq) == "ATNGC"
    assert fasta_record.id.startswith("q1_ACC1_Species_GBIF_99_50")


def test_gbif_failure_falls_back_to_ncbi(tmp_path, monkeypatch):
    # GBIF失敗時のNCBI fallbackを確認する / Verify NCBI fallback when GBIF fails.
    input_path = tmp_path / "input.csv"
    write_blast_csv(input_path, [{"qseqid": "q1", "sallacc": "ACC1", "pident": 99.0, "qseq": "ATGC"}])
    module = load_phylo_module(input_path, monkeypatch, "--onlyp")
    module["fetch_ncbi_data"] = lambda accession: ncbi_record()
    module["get_gbif_taxonomic_info"] = lambda species: None

    _, csv_entry = module["process_row"](0, module["df"].iloc[0])

    assert csv_entry[5] == "Mock species"
    assert csv_entry[8] == "NCBI"


def test_gbif_missing_required_taxonomy_uses_fallback(tmp_path, monkeypatch):
    # GBIF必要rank欠損時の現行fallbackを確認する / Verify fallback for incomplete GBIF taxonomy.
    input_path = tmp_path / "input.csv"
    write_blast_csv(input_path, [{"qseqid": "q1", "sallacc": "ACC1", "pident": 99.0, "qseq": "ATGC"}])
    module = load_phylo_module(input_path, monkeypatch, "--onlyp")
    module["fetch_ncbi_data"] = lambda accession: ncbi_record()

    class Response:
        status_code = 200

        @staticmethod
        def json():
            return {"species": "Incomplete", "genus": "Genus", "family": "Family"}

    monkeypatch.setattr(module["requests"], "get", lambda *args, **kwargs: Response())
    monkeypatch.setattr(module["time"], "sleep", lambda seconds: None)

    _, csv_entry = module["process_row"](0, module["df"].iloc[0])

    assert csv_entry[5] == "Mock species"
    assert csv_entry[8] == "NCBI"


def test_ncbi_failure_returns_current_uncertain_taxonomy_behavior(tmp_path, monkeypatch):
    # NCBI失敗時の現行出力を確認する / Verify current output when NCBI retrieval fails.
    input_path = tmp_path / "input.csv"
    write_blast_csv(input_path, [{"qseqid": "q1", "sallacc": "ACC1", "pident": 99.0, "qseq": "ATGC"}])
    module = load_phylo_module(input_path, monkeypatch, "--onlyp")
    module["fetch_ncbi_data"] = lambda accession: None

    _, csv_entry = module["process_row"](0, module["df"].iloc[0])

    assert csv_entry[5] == "Uncertain_taxonomy"
    assert csv_entry[8] == "NCBI Failed"


def test_fetch_ncbi_retries_then_returns_records(tmp_path, monkeypatch):
    # NCBI retry回数と成功結果を確認する / Verify NCBI retries and eventual success.
    input_path = tmp_path / "input.csv"
    write_blast_csv(input_path, [{"qseqid": "q1", "sallacc": "ACC1", "pident": 99.0, "qseq": "ATGC"}])
    module = load_phylo_module(input_path, monkeypatch, "--onlyp")
    attempts = []

    def efetch(**kwargs):
        attempts.append(kwargs["id"])
        if len(attempts) < 3:
            raise RuntimeError("temporary failure")
        class Handle:
            def close(self):
                pass

        return Handle()

    monkeypatch.setattr(module["Entrez"], "efetch", efetch)
    monkeypatch.setattr(module["Entrez"], "read", lambda handle: ncbi_record())
    monkeypatch.setattr(module["time"], "sleep", lambda seconds: None)

    assert module["fetch_ncbi_data"]("ACC1") == ncbi_record()
    assert attempts == ["ACC1", "ACC1", "ACC1"]


def test_gbif_request_failure_returns_none_after_retries(tmp_path, monkeypatch):
    # GBIF通信失敗時のretryとNone返却を確認する / Verify GBIF retries and None after request failures.
    input_path = tmp_path / "input.csv"
    write_blast_csv(input_path, [{"qseqid": "q1", "sallacc": "ACC1", "pident": 99.0, "qseq": "ATGC"}])
    module = load_phylo_module(input_path, monkeypatch, "--onlyp")
    attempts = []

    def get(*args, **kwargs):
        attempts.append(args[0])
        raise module["requests"].exceptions.RequestException("offline")

    monkeypatch.setattr(module["requests"], "get", get)
    monkeypatch.setattr(module["time"], "sleep", lambda seconds: None)

    assert module["get_gbif_taxonomic_info"]("Mock species") is None
    assert len(attempts) == 3


def test_gbif_response_is_parsed_and_retries_until_success(phylo_module, monkeypatch):
    response_data = json.loads((FIXTURES / "taxonomy/mock_gbif.json").read_text())
    get = Mock(side_effect=[
        Mock(status_code=503), Mock(status_code=200, json=lambda: response_data)
    ])
    sleep = Mock()
    monkeypatch.setattr(phylo_module["requests"], "get", get)
    monkeypatch.setattr(phylo_module["time"], "sleep", sleep)

    assert phylo_module["get_gbif_taxonomic_info"]("Mock species") == response_data
    assert get.call_count == 2
    get.assert_called_with("https://api.gbif.org/v1/species/match?name=Mock species", timeout=30)
    sleep.assert_called_once_with(1)


@pytest.mark.parametrize("response", [{}, {"order": None}])
def test_gbif_missing_order_retries_three_times(phylo_module, monkeypatch, response):
    get = Mock(return_value=Mock(status_code=200, json=lambda: response))
    sleep = Mock()
    monkeypatch.setattr(phylo_module["requests"], "get", get)
    monkeypatch.setattr(phylo_module["time"], "sleep", sleep)
    assert phylo_module["get_gbif_taxonomic_info"]("Mock species") is None
    assert get.call_count == 3
    assert [call.args for call in sleep.call_args_list] == [(1,), (1,), (1,)]


@pytest.mark.parametrize("pident,expected", [
    (99, "Mock species"), (96, "Genusa"), (93, "Familya"), (86, "Ordera")
])
def test_ncbi_fallback_uses_positions_in_mock_lineage(phylo_module, monkeypatch, pident, expected):
    # This records positional extraction, not biological rank correctness.
    phylo_module["fetch_ncbi_data"] = lambda accession: ncbi_record()
    get = Mock(side_effect=phylo_module["requests"].exceptions.RequestException("offline"))
    monkeypatch.setattr(phylo_module["requests"], "get", get)
    monkeypatch.setattr(phylo_module["time"], "sleep", lambda seconds: None)
    row = pd.Series({"qseqid": "q1", "sallacc": "ACC1", "pident": pident, "qseq": "ATGC"})
    _, entry = phylo_module["process_row"](0, row)
    assert entry[2:6] == ["Insecta", "Ordera", "Familya", expected]
    assert entry[8] == "NCBI"
    assert get.call_count == 3


def test_both_services_unavailable_short_circuits_after_ncbi_retries(phylo_module, monkeypatch):
    fetch = Mock(side_effect=RuntimeError("NCBI offline"))
    get = Mock(side_effect=phylo_module["requests"].exceptions.RequestException("GBIF offline"))
    sleep = Mock()
    monkeypatch.setattr(phylo_module["Entrez"], "efetch", fetch)
    monkeypatch.setattr(phylo_module["requests"], "get", get)
    monkeypatch.setattr(phylo_module["time"], "sleep", sleep)

    fasta, entry = phylo_module["process_row"](0, phylo_module["df"].iloc[0])
    assert entry == ["q1", "ACC1", "Unknown", "Unknown", "Unknown", "Uncertain_taxonomy", 99.0, "ATGC", "NCBI Failed"]
    assert fasta == ">q1_ACC1_Uncertain_taxonomy_99_00\nATGC\n"
    assert fetch.call_count == 3
    assert [call.args for call in sleep.call_args_list] == [(1,), (2,), (4,)]
    get.assert_not_called()


def test_taxonomy_pipeline_writes_matching_csv_and_fasta(phylo_module):
    phylo_module["df"] = pd.DataFrame([
        {"qseqid": "q1", "sallacc": "ACC1", "pident": 99.5, "qseq": "AT-GC"},
        {"qseqid": "q2", "sallacc": "ACC2", "pident": 96.0, "qseq": "GGTA"},
    ])
    phylo_module["fetch_ncbi_data"] = lambda accession: ncbi_record()
    data = json.loads((FIXTURES / "taxonomy/mock_gbif.json").read_text())
    phylo_module["get_gbif_taxonomic_info"] = lambda species: data
    fasta_path = phylo_module["process_with_progress"]()
    from pathlib import Path
    csv_path = Path(phylo_module["taxonomy_dir"]) / "taxonomic_data.csv"
    with csv_path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        assert reader.fieldnames == ["qseqid", "accessionID", "class", "order", "family", "taxonomic_name", "pident", "qseq", "source"]
        rows = {row["qseqid"]: row for row in reader}
    assert set(rows) == {"q1", "q2"}
    records = list(SeqIO.parse(fasta_path, "fasta"))
    assert len(records) == 2
    assert {record.id: str(record.seq) for record in records} == {
        "q1_ACC1_Species_GBIF_99_50": "ATNGC", "q2_ACC2_Genus_GBIF_96_00": "GGTA"
    }
    for qseqid, accession, name, sequence in [
        ("q1", "ACC1", "Species GBIF", "ATNGC"), ("q2", "ACC2", "Genus GBIF", "GGTA")
    ]:
        assert rows[qseqid]["accessionID"] == accession
        assert rows[qseqid]["taxonomic_name"] == name
        assert rows[qseqid]["qseq"] == sequence
        assert rows[qseqid]["source"] == "GBIF"
