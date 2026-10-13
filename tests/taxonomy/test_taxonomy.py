"""taxonomy API の fallback/retry、pident 境界、CSV/FASTA metadata 整合を確認する。"""

import csv
import io
import json
from concurrent.futures import ThreadPoolExecutor
from urllib.error import HTTPError

import pandas as pd
import pytest
from Bio import SeqIO

from conftest import FIXTURES, load_phylo_module, write_blast_csv
from unittest.mock import Mock
from phylomip.taxonomy import CSV_COLUMNS, NCBIRequestLimiter, TaxonomyProcessor
import phylomip.taxonomy as taxonomy_module


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
    assert list(parsed_csv.columns[:8]) == [
        "qseqid", "accessionID", "class", "order", "family", "taxonomic_name", "pident", "qseq"
    ]
    assert set(CSV_COLUMNS) == set(parsed_csv.columns)
    assert len(parsed_csv) == 1
    assert parsed_csv.loc[0, "qseqid"] == "q1"
    assert parsed_csv.loc[0, "accessionID"] == "ACC1"
    assert parsed_csv.loc[0, "taxonomic_name"] == "Species GBIF"
    assert parsed_csv.loc[0, "resolved_source"] == "GBIF"
    assert parsed_csv.loc[0, "gbif_species"] == "Species GBIF"
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
    assert dict(zip(CSV_COLUMNS, csv_entry))["resolved_source"] == "NCBI_lineage_fallback"


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
    assert dict(zip(CSV_COLUMNS, csv_entry))["resolved_source"] == "NCBI_lineage_fallback"


def test_ncbi_failure_returns_current_uncertain_taxonomy_behavior(tmp_path, monkeypatch):
    # NCBI失敗時の現行出力を確認する / Verify current output when NCBI retrieval fails.
    input_path = tmp_path / "input.csv"
    write_blast_csv(input_path, [{"qseqid": "q1", "sallacc": "ACC1", "pident": 99.0, "qseq": "ATGC"}])
    module = load_phylo_module(input_path, monkeypatch, "--onlyp")
    module["fetch_ncbi_data"] = lambda accession: None

    _, csv_entry = module["process_row"](0, module["df"].iloc[0])

    assert csv_entry[5] == "Uncertain_taxonomy"
    assert dict(zip(CSV_COLUMNS, csv_entry))["resolved_source"] == "unresolved"


def test_ncbi_taxid_is_resolved_by_rank_name_and_saved_separately(phylo_module):
    records = [{
        "GBSeq_organism": "NCBI species",
        "GBSeq_taxonomy": "legacy; lineage; values",
        "GBSeq_feature-table": [{
            "GBFeature_quals": [{"GBQualifier_name": "db_xref", "GBQualifier_value": "taxon:123"}]
        }],
    }]
    ncbi_taxonomy = {
        "kingdom": "Animalia", "phylum": "Chordata", "class": "Amphibia",
        "order": "Anura", "family": "Ranidae", "genus": "Rana", "species": "Rana temporaria",
    }
    phylo_module["fetch_ncbi_data"] = lambda accession: records
    phylo_module["fetch_ncbi_taxonomy"] = lambda taxid: ncbi_taxonomy
    phylo_module["get_gbif_taxonomic_info"] = lambda species: None

    _, entry = phylo_module["process_row"](0, phylo_module["df"].iloc[0])

    values = dict(zip(CSV_COLUMNS, entry))
    assert TaxonomyProcessor.extract_ncbi_taxid(records) == "123"
    assert values["ncbi_taxid"] == "123"
    assert values["ncbi_class"] == "Amphibia"
    assert values["resolved_species"] == "Rana temporaria"
    assert values["resolved_source"] == "NCBI"
    assert values["taxonomic_name"] == "Rana temporaria"


def test_ncbi_first_resolution_records_rank_conflicts(phylo_module):
    records = [{
        "GBSeq_organism": "NCBI species",
        "GBSeq_taxonomy": "legacy; lineage; values",
        "GBSeq_feature-table": [{
            "GBFeature_quals": [{"GBQualifier_name": "db_xref", "GBQualifier_value": "taxon:456"}]
        }],
    }]
    phylo_module["fetch_ncbi_data"] = lambda accession: records
    phylo_module["fetch_ncbi_taxonomy"] = lambda taxid: {
        "class": "Amphibia", "order": "Anura", "family": "Ranidae",
        "genus": "Rana", "species": "Rana temporaria",
    }
    phylo_module["get_gbif_taxonomic_info"] = lambda species: {
        "source": "GBIF", "class": "Amphibia", "order": "Anura", "family": "Ranidae",
        "genus": "Pelophylax", "species": "Pelophylax nigromaculatus",
    }

    _, entry = phylo_module["process_row"](0, phylo_module["df"].iloc[0])

    values = dict(zip(CSV_COLUMNS, entry))
    assert values["resolved_source"] == "NCBI"
    assert values["resolved_genus"] == "Rana"
    assert values["resolved_species"] == "Rana temporaria"
    assert values["gbif_genus"] == "Pelophylax"
    assert values["taxonomy_conflict"] == "true"
    assert values["taxonomy_conflict_ranks"] == "genus;species"
    assert "source" not in values


def test_ncbi_taxonomy_failure_uses_gbif_as_whole_record_fallback(phylo_module):
    records = [{
        "GBSeq_organism": "NCBI species",
        "GBSeq_taxonomy": "legacy; lineage; values",
        "GBSeq_feature-table": [{
            "GBFeature_quals": [{"GBQualifier_name": "db_xref", "GBQualifier_value": "taxon:999"}]
        }],
    }]
    phylo_module["fetch_ncbi_data"] = lambda accession: records
    phylo_module["fetch_ncbi_taxonomy"] = lambda taxid: None
    phylo_module["get_gbif_taxonomic_info"] = lambda species: {
        "source": "GBIF", "class": "Insecta", "order": "Coleoptera", "family": "Carabidae",
        "genus": "Carabus", "species": "Carabus auratus",
    }

    _, entry = phylo_module["process_row"](0, phylo_module["df"].iloc[0])

    values = dict(zip(CSV_COLUMNS, entry))
    assert values["ncbi_taxid"] == "999"
    assert values["resolved_source"] == "GBIF"
    assert values["resolved_genus"] == "Carabus"
    assert values["taxonomic_name"] == "Carabus auratus"


def test_ncbi_taxonomy_and_gbif_failure_marks_lineage_fallback(phylo_module):
    records = [{
        "GBSeq_organism": "NCBI species",
        "GBSeq_taxonomy": "cellular organisms; Eukaryota; Animalia; Arthropoda; Insecta; Ordera; Familya; Genusa",
        "GBSeq_feature-table": [{
            "GBFeature_quals": [{"GBQualifier_name": "db_xref", "GBQualifier_value": "taxon:1000"}]
        }],
    }]
    phylo_module["fetch_ncbi_data"] = lambda accession: records
    phylo_module["fetch_ncbi_taxonomy"] = lambda taxid: None
    phylo_module["get_gbif_taxonomic_info"] = lambda species: None

    _, entry = phylo_module["process_row"](0, phylo_module["df"].iloc[0])

    values = dict(zip(CSV_COLUMNS, entry))
    assert values["resolved_source"] == "NCBI_lineage_fallback"
    assert values["resolved_class"] == "Insecta"
    assert values["resolved_genus"] == "Genusa"


def test_missing_rank_is_written_as_blank_without_shifting_later_ranks(phylo_module):
    phylo_module["fetch_ncbi_data"] = lambda accession: [{
        "GBSeq_organism": "NCBI species",
        "GBSeq_taxonomy": "legacy; lineage; values",
        "GBSeq_feature-table": [{
            "GBFeature_quals": [{"GBQualifier_name": "db_xref", "GBQualifier_value": "taxon:789"}]
        }],
    }]
    phylo_module["fetch_ncbi_taxonomy"] = lambda taxid: {
        "class": "Amphibia", "order": None, "family": "Ranidae",
        "genus": "Rana", "species": "Rana temporaria",
    }
    phylo_module["get_gbif_taxonomic_info"] = lambda species: None

    _, entry = phylo_module["process_row"](0, phylo_module["df"].iloc[0])

    values = dict(zip(CSV_COLUMNS, entry))
    assert values["resolved_class"] == "Amphibia"
    assert values["resolved_order"] is None
    assert values["resolved_family"] == "Ranidae"
    assert values["taxonomic_name"] == "Rana temporaria"


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


def test_fetch_ncbi_taxonomy_parses_rank_names_and_uses_cache(tmp_path, monkeypatch):
    processor = TaxonomyProcessor(tmp_path)
    calls = []

    class Handle:
        def close(self):
            pass

    def efetch(**kwargs):
        calls.append(kwargs)
        return Handle()

    taxonomy_response = [{
        "LineageEx": [
            {"Rank": "kingdom", "ScientificName": "Animalia"},
            {"Rank": "class", "ScientificName": "Amphibia"},
            {"Rank": "genus", "ScientificName": "Rana"},
        ],
        "Rank": "species",
        "ScientificName": "Rana temporaria",
    }]
    monkeypatch.setattr(taxonomy_module.Entrez, "efetch", efetch)
    monkeypatch.setattr(taxonomy_module.Entrez, "read", lambda handle: taxonomy_response)
    monkeypatch.setattr(taxonomy_module.time, "sleep", lambda seconds: None)

    first = processor.fetch_ncbi_taxonomy("123")
    second = processor.fetch_ncbi_taxonomy("123")

    assert first["kingdom"] == "Animalia"
    assert first["class"] == "Amphibia"
    assert first["genus"] == "Rana"
    assert first["species"] == "Rana temporaria"
    assert second == first
    assert len(calls) == 1
    assert calls[0] == {"db": "taxonomy", "id": "123", "retmode": "xml"}


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
    assert dict(zip(CSV_COLUMNS, entry))["resolved_source"] == "NCBI_lineage_fallback"
    assert get.call_count == 3


def test_both_services_unavailable_short_circuits_after_ncbi_retries(phylo_module, monkeypatch):
    fetch = Mock(side_effect=RuntimeError("NCBI offline"))
    get = Mock(side_effect=phylo_module["requests"].exceptions.RequestException("GBIF offline"))
    sleep = Mock()
    monkeypatch.setattr(phylo_module["Entrez"], "efetch", fetch)
    monkeypatch.setattr(phylo_module["requests"], "get", get)
    monkeypatch.setattr(phylo_module["time"], "sleep", sleep)

    fasta, entry = phylo_module["process_row"](0, phylo_module["df"].iloc[0])
    assert entry[:8] == ["q1", "ACC1", None, None, None, "Uncertain_taxonomy", 99.0, "ATGC"]
    values = dict(zip(CSV_COLUMNS, entry))
    assert values["resolved_source"] == "unresolved"
    assert values["taxonomy_conflict"] == "false"
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
        assert list(reader.fieldnames[:8]) == [
            "qseqid", "accessionID", "class", "order", "family", "taxonomic_name", "pident", "qseq"
        ]
        assert set(reader.fieldnames) == set(CSV_COLUMNS)
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
        assert rows[qseqid]["resolved_source"] == "GBIF"


def test_ncbi_limiter_is_shared_across_threads(monkeypatch):
    # 全スレッドで同じ待機間隔を共有することを確認する / Verify one limiter is shared by all threads.
    limiter = NCBIRequestLimiter(interval=0.4)
    sleeps = Mock()
    monkeypatch.setattr(taxonomy_module.time, "monotonic", lambda: 0.0)
    monkeypatch.setattr(taxonomy_module.time, "sleep", sleeps)
    with ThreadPoolExecutor(max_workers=2) as executor:
        list(executor.map(lambda _: limiter.wait(), range(2)))
    assert sleeps.call_count == 1
    assert sleeps.call_args.args == (0.4,)


def test_nucleotide_and_taxonomy_requests_use_combined_limiter(tmp_path, monkeypatch):
    # nucleotide/taxonomyを別々に数えず共通制限することを確認する / Verify combined throttling.
    processor = TaxonomyProcessor(tmp_path, ncbi_limiter=NCBIRequestLimiter(interval=0.4))
    sleeps = Mock()
    monkeypatch.setattr(taxonomy_module.time, "monotonic", lambda: 0.0)
    monkeypatch.setattr(taxonomy_module.time, "sleep", sleeps)

    class Handle:
        def close(self):
            pass

    monkeypatch.setattr(taxonomy_module.Entrez, "efetch", Mock(return_value=Handle()))
    monkeypatch.setattr(taxonomy_module.Entrez, "read", Mock(side_effect=[
        ncbi_record(), [{"Rank": "species", "ScientificName": "Mock species"}]
    ]))
    processor.fetch_ncbi_data("ACC1")
    processor.fetch_ncbi_taxonomy("123")
    assert sleeps.call_count == 1
    assert sleeps.call_args.args == (0.4,)


def test_prepare_ncbi_data_deduplicates_accessions_and_taxids(tmp_path):
    # 重複accession/TaxIDの取得を一度に抑えることを確認する / Verify request deduplication.
    processor = TaxonomyProcessor(tmp_path, ncbi_limiter=NCBIRequestLimiter(interval=0))
    rows = pd.DataFrame([
        {"sallacc": "ACC1", "qseqid": "q1"},
        {"sallacc": "ACC1", "qseqid": "q2"},
        {"sallacc": "ACC2", "qseqid": "q3"},
    ])
    nucleotide_calls, taxonomy_calls = [], []

    def fetch_data(accession):
        nucleotide_calls.append(accession)
        records = [{"GBSeq_organism": "Mock", "GBSeq_feature-table": [
            {"GBFeature_quals": [{"GBQualifier_name": "db_xref", "GBQualifier_value": "taxon:1"}]}
        ]}]
        processor.cache[accession] = records
        return records

    def fetch_taxonomy(taxid):
        taxonomy_calls.append(taxid)
        processor.taxonomy_cache[str(taxid)] = {"class": "Insecta"}
        return processor.taxonomy_cache[str(taxid)]

    processor.fetch_ncbi_data = fetch_data
    processor.fetch_ncbi_taxonomy = fetch_taxonomy
    processor.prepare_ncbi_data(rows)
    assert nucleotide_calls == ["ACC1", "ACC2"]
    assert taxonomy_calls == ["1"]


def test_accession_versions_are_deduplicated_and_cached_separately(tmp_path):
    # accession.versionを保持し、異なるversionを別recordとして扱うことを確認する。
    processor = TaxonomyProcessor(tmp_path, ncbi_limiter=NCBIRequestLimiter(interval=0))
    rows = pd.DataFrame([
        {"sallacc": "ACC1.1"}, {"sallacc": "ACC1.2"}, {"sallacc": "ACC1.1"},
    ])
    calls = []
    records = {
        "ACC1.1": [{"GBSeq_accession-version": "ACC1.1", "GBSeq_organism": "Version one"}],
        "ACC1.2": [{"GBSeq_accession-version": "ACC1.2", "GBSeq_organism": "Version two"}],
    }

    def fetch_data(accession):
        calls.append(accession)
        processor.cache[accession] = records[accession]
        return records[accession]

    processor.fetch_ncbi_data = fetch_data
    processor.fetch_ncbi_taxonomy = lambda taxid: None
    processor.prepare_ncbi_data(rows)
    assert calls == ["ACC1.1", "ACC1.2"]
    assert processor.cache["ACC1.1"][0]["GBSeq_organism"] == "Version one"
    assert processor.cache["ACC1.2"][0]["GBSeq_organism"] == "Version two"


def test_prepare_ncbi_data_uses_batch_and_falls_back_for_partial_results(tmp_path):
    # 大規模入力のHistory batchと未取得分の個別fallbackを確認する / Verify batch mode and partial fallback.
    processor = TaxonomyProcessor(tmp_path, ncbi_limiter=NCBIRequestLimiter(interval=0))
    rows = pd.DataFrame([{"sallacc": "ACC{:03d}".format(i)} for i in range(51)])
    batch_accessions, direct_accessions, batch_taxids, direct_taxids = [], [], [], []

    def batch_nucleotide(accessions):
        batch_accessions.append(list(accessions))
        for accession in accessions[:-1]:
            processor.cache[accession] = [{
                "GBSeq_accession-version": accession,
                "GBSeq_feature-table": [{"GBFeature_quals": [
                    {"GBQualifier_name": "db_xref", "GBQualifier_value": "taxon:1"}
                ]}],
            }]
        return {accessions[-1]}

    def direct_nucleotide(accession):
        direct_accessions.append(accession)
        processor.cache[accession] = [{
            "GBSeq_accession-version": accession,
            "GBSeq_feature-table": [{"GBFeature_quals": [
                {"GBQualifier_name": "db_xref", "GBQualifier_value": "taxon:2"}
            ]}],
        }]
        return processor.cache[accession]

    def batch_taxonomy(taxids):
        batch_taxids.append(list(taxids))
        processor.taxonomy_cache["1"] = {"class": "Insecta"}
        return {"2"}

    def direct_taxonomy(taxid):
        direct_taxids.append(taxid)
        processor.taxonomy_cache[str(taxid)] = {"class": "Insecta"}
        return processor.taxonomy_cache[str(taxid)]

    processor._prefetch_nucleotide_batch = batch_nucleotide
    processor._prefetch_taxonomy_batch = batch_taxonomy
    processor.fetch_ncbi_data = direct_nucleotide
    processor.fetch_ncbi_taxonomy = direct_taxonomy
    processor.prepare_ncbi_data(rows)
    assert len(batch_accessions) == 1
    assert direct_accessions == ["ACC050"]
    assert batch_taxids == [["1", "2"]]
    assert direct_taxids == ["2"]


def test_batch_records_are_mapped_to_accession_and_taxid_caches(tmp_path):
    # History batchの応答を各行のcacheへ対応付けることを確認する / Verify batch mapping into caches.
    processor = TaxonomyProcessor(tmp_path, ncbi_limiter=NCBIRequestLimiter(interval=0))
    nucleotide = [{
        "GBSeq_accession-version": "ACC1.1",
        "GBSeq_feature-table": [{"GBFeature_quals": [
            {"GBQualifier_name": "db_xref", "GBQualifier_value": "taxon:123"}
        ]}],
    }]
    taxonomy = [{
        "TaxId": "123",
        "LineageEx": [{"Rank": "class", "ScientificName": "Insecta"}],
        "Rank": "species", "ScientificName": "Mock species",
    }]
    second_nucleotide = {
        "GBSeq_accession-version": "ACC1.2",
        "GBSeq_feature-table": [{"GBFeature_quals": [
            {"GBQualifier_name": "db_xref", "GBQualifier_value": "taxon:456"}
        ]}],
    }
    processor._history_fetch = Mock(side_effect=[[nucleotide[0], second_nucleotide], taxonomy])
    assert processor._prefetch_nucleotide_batch(["ACC1.1", "ACC1.2"]) == set()
    assert processor._prefetch_taxonomy_batch(["123", "456"]) == {"456"}
    assert set(processor.cache) == {"ACC1.1", "ACC1.2"}
    assert processor._record_taxid(processor.cache["ACC1.1"][0]) == "123"
    assert processor._record_taxid(processor.cache["ACC1.2"][0]) == "456"
    assert processor.taxonomy_cache["123"]["class"] == "Insecta"
    assert processor.taxonomy_cache["123"]["species"] == "Mock species"


def test_ncbi_429_honors_retry_after(tmp_path, monkeypatch):
    # 429時にRetry-Afterを優先することを確認する / Verify Retry-After-aware backoff.
    processor = TaxonomyProcessor(tmp_path, ncbi_limiter=NCBIRequestLimiter(interval=0))
    error = HTTPError("https://eutils.ncbi.nlm.nih.gov", 429, "Too Many Requests",
                      {"Retry-After": "7"}, None)

    class Handle:
        def close(self):
            pass

    monkeypatch.setattr(taxonomy_module.Entrez, "efetch", Mock(side_effect=[error, Handle()]))
    monkeypatch.setattr(taxonomy_module.Entrez, "read", Mock(return_value=ncbi_record()))
    sleep = Mock()
    monkeypatch.setattr(taxonomy_module.time, "sleep", sleep)
    assert processor.fetch_ncbi_data("ACC1") == ncbi_record()
    assert sleep.call_args_list == [((7.0,), {})]


def test_ncbi_request_interval_defaults_and_conservative_override(monkeypatch):
    # API key有無の安全側既定値と、長くするoverrideを確認する。
    monkeypatch.delenv("NCBI_API_KEY", raising=False)
    monkeypatch.delenv("NCBI_REQUEST_INTERVAL", raising=False)
    assert NCBIRequestLimiter().interval == 0.4

    monkeypatch.setenv("NCBI_API_KEY", "test-key")
    assert NCBIRequestLimiter().interval == 0.12
    monkeypatch.setenv("NCBI_REQUEST_INTERVAL", "0.8")
    assert NCBIRequestLimiter().interval == 0.8
    monkeypatch.setenv("NCBI_REQUEST_INTERVAL", "0.01")
    assert NCBIRequestLimiter().interval == 0.12


def test_readmes_separate_blast_and_entrez_rate_guidance():
    # BLASTとtaxonomy用Entrezのrate limit説明を混同しないことを確認する。
    from pathlib import Path
    project_root = Path(__file__).parents[2]
    for readme_name in ("README.md", "README-jp.md"):
        text = (project_root / readme_name).read_text()
        blast = text.index("BLAST")
        entrez = text.index("Entrez")
        assert blast < entrez
        assert "0.4" in text
        assert "0.12" in text
        assert "NCBI_BATCH_THRESHOLD=100" in text
