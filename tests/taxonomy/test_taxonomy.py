import csv
import io
import json

import pandas as pd
import pytest
from Bio import SeqIO

from conftest import load_phylo_module, write_blast_csv


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