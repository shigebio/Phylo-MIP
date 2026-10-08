"""実際の bPTP/mPTP が出力する成果物の存在・parse 可否・species 対応を確認する。"""

from pathlib import Path

import pytest
from Bio import Phylo

from conftest import FIXTURES

pytestmark = pytest.mark.integration


def test_real_bptp_artifacts_and_parser(phylo_module, require_tools):
    require_tools("xvfb-run", "python3")
    if not Path("/app/PTP/bin/bPTP.py").is_file():
        pytest.skip("bPTP wrapper requires /app/PTP/bin/bPTP.py; run in the project image")
    tree_path = FIXTURES / "phylogeny/tree.nwk"
    phylo_module["run_bptp"](str(tree_path), 1000, 10, 0.1, 42, phylo_module["bptp_base_dir"])
    timestamp = phylo_module["timestamp"]
    directory = Path(phylo_module["bptp_base_dir"]) / f"{timestamp}_bPTP_analysis"
    prefix = directory / f"{timestamp}_bPTP_species_delimitation"
    for suffix in [".PTPMLPartition.txt", ".PTPhSupportPartition.txt", ".PTPPartitions.txt",
                   ".PTPPartitonSummary.txt", ".PTPllh.txt",
                   ".PTPMLPartition.txt.png", ".PTPMLPartition.txt.svg",
                   ".PTPhSupportPartition.txt.png", ".PTPhSupportPartition.txt.svg", ".llh.pdf"]:
        assert Path(str(prefix) + suffix).stat().st_size > 0
    for suffix in [".PTPMLPartition.txt", ".PTPhSupportPartition.txt"]:
        species, support = phylo_module["extract_species_data"](str(prefix) + suffix)
        assert set(species) == {"q1", "q2", "q3", "q4"}
        assert set(support) == set(species)
        assert all(value.isdigit() for value in species.values())
        assert all(value is None or 0 <= float(value) <= 1 for value in support.values())
    expected = {leaf.name for leaf in Phylo.read(tree_path, "newick").get_terminals()}
    for suffix in [".PTPMLPartition.txt.ml.tre", ".PTPhSupportPartition.txt.sh.tre"]:
        terminals = Phylo.read(str(prefix) + suffix, "newick").get_terminals()
        assert {leaf.name for leaf in terminals} == expected
        assert len(terminals) == len(expected)


def test_real_mptp_artifacts_and_parser(phylo_module, require_tools):
    require_tools("xvfb-run", "mptp")
    phylo_module["run_mptp"](str(FIXTURES / "phylogeny/tree.nwk"), phylo_module["mptp_base_dir"])
    timestamp = phylo_module["timestamp"]
    directory = Path(phylo_module["mptp_base_dir"]) / f"{timestamp}_mPTP_analysis"
    prefix = directory / f"{timestamp}_mPTP_species_delimitation"
    assert Path(str(prefix) + ".txt").stat().st_size > 0
    assert Path(str(prefix) + ".svg").stat().st_size > 0
    species, support = phylo_module["extract_species_data"](str(prefix) + ".txt")
    assert set(species) == {"q1", "q2", "q3", "q4"}
    assert all(value.isdigit() for value in species.values())
    assert support == dict.fromkeys(species)
