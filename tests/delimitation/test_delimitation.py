"""bPTP/mPTP の partition parser と、taxonomy CSV への反映を確認する。"""

import csv
from pathlib import Path

import pytest

from conftest import FIXTURES


@pytest.mark.parametrize("filename,species,support", [
    ("bptp.PTPhSupportPartition.txt", {"q1": "1", "q2": "1", "q3": "2"},
     {"q1": "0.95", "q2": "0.95", "q3": "0.80"}),
    ("bptp.PTPMLPartition.txt", {"q1": "1", "q2": "2", "q3": "2"},
     {"q1": "0.90", "q2": "0.75", "q3": "0.75"}),
    ("mptp_species_delimitation.txt", {"q1": "3", "q2": "3", "q3": "4"},
     {"q1": None, "q2": None, "q3": None}),
])
def test_partition_parser_matches_fixture(phylo_module, filename, species, support):
    assert phylo_module["extract_species_data"](FIXTURES / "delimitation" / filename) == (species, support)


def test_delimitation_wrappers_discovery_and_csv_update(phylo_module, monkeypatch):
    timestamp = phylo_module["timestamp"]
    tree_file = str(FIXTURES / "phylogeny/tree.nwk")
    bptp_dir = Path(phylo_module["bptp_base_dir"]) / f"{timestamp}_bPTP_analysis"
    mptp_dir = Path(phylo_module["mptp_base_dir"]) / f"{timestamp}_mPTP_analysis"
    bptp_prefix = bptp_dir / f"{timestamp}_bPTP_species_delimitation"
    mptp_prefix = mptp_dir / f"{timestamp}_mPTP_species_delimitation"
    calls = []

    def fake_run(command, check):
        assert check
        calls.append(command)
        if "mptp" in command:
            assert command == ["xvfb-run", "-a", "mptp", "-tree_file", tree_file,
                               "-output_file", str(mptp_prefix), "-ml", "-single"]
            prefix = Path(command[command.index("-output_file") + 1])
            assert prefix.parent.is_dir()  # Created by the production wrapper.
            Path(str(prefix) + ".txt").write_bytes(
                (FIXTURES / "delimitation/mptp_species_delimitation.txt").read_bytes())
            Path(str(prefix) + ".svg").write_text('<svg xmlns="http://www.w3.org/2000/svg"/>')
        else:
            assert command == ["xvfb-run", "-a", "python3", "/app/PTP/bin/bPTP.py",
                               "-t", tree_file, "-o", str(bptp_prefix),
                               "-s", "42", "-i", "1000", "-n", "10", "-b", "0.1"]
            prefix = Path(command[command.index("-o") + 1])
            assert prefix.parent.is_dir()
            for suffix in [".PTPMLPartition.txt", ".PTPhSupportPartition.txt"]:
                Path(str(prefix) + suffix).write_bytes(
                    (FIXTURES / "delimitation" / ("bptp" + suffix)).read_bytes())
            # Opaque external artifacts; real-tool tests validate actual production.
            for suffix in [".PTPPartitions.txt", ".PTPPartitonSummary.txt", ".PTPllh.txt"]:
                Path(str(prefix) + suffix).write_text("fixture external artifact\n")

    monkeypatch.setattr(phylo_module["subprocess"], "run", fake_run)
    phylo_module["run_bptp"](tree_file, 1000, 10, 0.1, 42, phylo_module["bptp_base_dir"])
    phylo_module["run_mptp"](tree_file, phylo_module["mptp_base_dir"])
    assert len(calls) == 2
    assert {p.name for p in bptp_dir.iterdir()} == {
        bptp_prefix.name + suffix for suffix in [
            ".PTPMLPartition.txt", ".PTPhSupportPartition.txt", ".PTPPartitions.txt",
            ".PTPPartitonSummary.txt", ".PTPllh.txt"]}
    assert {p.name for p in mptp_dir.iterdir()} == {mptp_prefix.name + ".txt", mptp_prefix.name + ".svg"}
    assert all(p.stat().st_size for directory in [bptp_dir, mptp_dir] for p in directory.iterdir())

    taxonomy = Path(phylo_module["taxonomy_dir"]) / "taxonomic_data.csv"
    taxonomy.write_text("qseqid,accessionID,taxonomic_name\nq1,ACC1,One\nq2,ACC2,Two\nq3,ACC3,Three\nq4,ACC4,Unassigned\n")
    phylo_module["process_ptp_outputs"]()
    with taxonomy.open(newline="") as handle:
        reader = csv.DictReader(handle)
        assert reader.fieldnames == [
            "qseqid", "accessionID", "taxonomic_name", "bPTP_Bayesian_Partitioned_Species",
            "PTPhSupport_support", "bPTP_ML_Partitioned_Species", "PTPML_support", "mPTP_Partitioned_Species"]
        rows = list(reader)
    assert len(rows) == 4
    for row, expected in zip(rows, [
        ("q1", "ACC1", "One", "1", "0.95", "1", "0.90", "3"),
        ("q2", "ACC2", "Two", "1", "0.95", "2", "0.75", "3"),
        ("q3", "ACC3", "Three", "2", "0.80", "2", "0.75", "4"),
        ("q4", "ACC4", "Unassigned", "", "", "", "", ""),
    ]):
        assert tuple(row.values()) == expected
