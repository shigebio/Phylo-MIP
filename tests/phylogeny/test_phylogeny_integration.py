"""実際の FastTree が入力 alignment の terminal taxa を保持することを確認する。"""

from pathlib import Path

import pytest
from Bio import Phylo, SeqIO

from conftest import FIXTURES

pytestmark = pytest.mark.integration


@pytest.mark.parametrize("method", ["NJ", "ML"])
def test_real_fasttree_preserves_taxa(phylo_module, require_tools, method):
    require_tools("fasttree", "mafft")
    aligned = phylo_module["run_mafft"](FIXTURES / "alignment/sequences.fasta", "aligned.fasta")
    result = phylo_module["run_fasttree"](aligned, "tree.nwk", method=method, bootstrap=0)
    assert Path(result).name == f"{phylo_module['timestamp']}_{method}_tree.nwk"
    terminals = Phylo.read(result, "newick").get_terminals()
    expected = {r.id for r in SeqIO.parse(aligned, "fasta")}
    assert {leaf.name for leaf in terminals} == expected
    assert len(terminals) == len(expected)
