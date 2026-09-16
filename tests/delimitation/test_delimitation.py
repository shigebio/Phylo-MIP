import csv
from pathlib import Path

from conftest import load_phylo_module, write_blast_csv


def load_module(tmp_path, monkeypatch):
    input_path = tmp_path / "input.csv"
    write_blast_csv(input_path, [{"qseqid": "q1", "sallacc": "ACC1", "pident": 99.0, "qseq": "ATGC"}])
    return load_phylo_module(input_path, monkeypatch, "--onlyp")


def test_bptp_parser_maps_species_and_support_values(tmp_path, monkeypatch):
    # bPTP partitionのspecies、support、OTU mappingを確認する / Verify bPTP species, support, and OTU mapping.
    module = load_module(tmp_path, monkeypatch)
    partition = tmp_path / "result.PTPhSupportPartition.txt"
    partition.write_text(
        "Species 1 (support = 0.95)\nq1_ACC1_Species_one_99_00, q2_ACC2_Species_two_98_00\n"
        "Species 2 (support = 0.80)\nq3_ACC3_Species_three_97_00\n",
        encoding="utf-8",
    )

    species, support = module["extract_species_data"](partition)

    assert species == {"q1": "1", "q2": "1", "q3": "2"}
    assert support == {"q1": "0.95", "q2": "0.95", "q3": "0.80"}


def test_delimitation_outputs_are_nonempty_and_update_taxonomy_csv(tmp_path, monkeypatch):
    # bPTP/mPTP成果物とtaxonomy CSVへのpartition反映を確認する / Verify delimitation artifacts and taxonomy CSV updates.
    module = load_module(tmp_path, monkeypatch)
    bptp_dir = tmp_path / "20260916_bPTP_analysis"
    mptp_dir = tmp_path / "20260916_mPTP_analysis"
    bptp_dir.mkdir()
    mptp_dir.mkdir()
    for filename in (
        "20260916_bPTP_species_delimitation.PTPMLPartition.txt",
        "20260916_bPTP_species_delimitation.PTPhSupportPartition.txt",
        "20260916_bPTP_species_delimitation.PTPPartitions.txt",
        "20260916_bPTP_species_delimitation.PTPPartitonSummary.txt",
        "20260916_bPTP_species_delimitation.PTPllh.txt",
    ):
        (bptp_dir / filename).write_text("Species 1\nq1_ACC1_Mock_99_00\n", encoding="utf-8")
    (mptp_dir / "20260916_mPTP_species_delimitation.txt").write_text(
        "Species 3\nq1_ACC1_Mock_99_00\n", encoding="utf-8"
    )
    (mptp_dir / "20260916_mPTP_species_delimitation.svg").write_text("<svg>mock</svg>", encoding="utf-8")
    taxonomy_csv = tmp_path / "taxonomy.csv"
    taxonomy_csv.write_text(
        "qseqid,accessionID,taxonomic_name\nq1,ACC1,Mock\n", encoding="utf-8"
    )

    assert all(path.stat().st_size > 0 for path in bptp_dir.iterdir())
    assert all(path.stat().st_size > 0 for path in mptp_dir.iterdir())
    module["update_csv_with_species_data"](
        taxonomy_csv,
        {"q1": "1"}, {"q1": "0.95"},
        {"q1": "2"}, {"q1": "0.90"},
        {"q1": "3"},
    )

    with taxonomy_csv.open(encoding="utf-8", newline="") as handle:
        row = next(csv.DictReader(handle))
    assert row["bPTP_Bayesian_Partitioned_Species"] == "1"
    assert row["PTPhSupport_support"] == "0.95"
    assert row["bPTP_ML_Partitioned_Species"] == "2"
    assert row["mPTP_Partitioned_Species"] == "3"