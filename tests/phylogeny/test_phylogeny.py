"""VSEARCH→MAFFT→FastTree の mocked pipeline と Newick/Nexus 出力を確認する。"""

import shlex
from pathlib import Path

import pytest
from Bio import Phylo

from conftest import load_phylo_module, write_blast_csv


@pytest.mark.parametrize("method", ["NJ", "ML"])
def test_tree_pipeline_writes_newick_and_nexus(tmp_path, monkeypatch, method):
    input_path = tmp_path / "input.csv"
    write_blast_csv(input_path, [{"qseqid": "q1", "sallacc": "ACC1", "pident": 99, "qseq": "ATGC"}])
    calls = []

    def fake_run(command, **kwargs):
        calls.append(command)
        assert kwargs["check"] is True
        if not isinstance(command, str):
            # Delimitation has separate wrapper and parser tests.
            assert command[:2] == ["xvfb-run", "-a"]
            return
        tokens = shlex.split(command)
        if tokens[0] == "vsearch":
            Path(tokens[tokens.index("--centroids") + 1]).write_text(">q1\nATGC\n>q3\nGGTA\n")
            Path(tokens[tokens.index("--mothur_shared_out") + 1]).write_text(
                "label\tGroup\tnumOtus\tOtu1\tOtu2\nunique\tall\t2\t1\t1\n")
        elif tokens[0] == "mafft":
            Path(tokens[-1]).write_bytes(Path(tokens[2]).read_bytes())
        else:
            assert tokens[:2] == ["fasttree", "-nt"]
            assert ("-nj" in tokens) == (method == "NJ")
            assert Path(tokens[-3]).is_file()
            Path(tokens[-1]).write_text("(q1:0.1,q3:0.2);\n")

    import subprocess
    monkeypatch.setattr(subprocess, "run", fake_run)
    module = load_phylo_module(input_path, monkeypatch, "--onlyp", "--tree", "--method", method, "--seed", "42")
    timestamp = module["timestamp"]
    # The CLI and wrapper both add a prefix; this existing behavior is documented.
    newick = Path(module["phylogeny_dir"]) / f"{timestamp}_{method}_{timestamp}_{method}_phylogenetic_tree.nwk"
    nexus = Path(module["phylogeny_dir"]) / f"{timestamp}_{method}_phylogenetic_tree.nex"
    for path, fmt in [(newick, "newick"), (nexus, "nexus")]:
        tree = Phylo.read(path, fmt)
        assert {leaf.name for leaf in tree.get_terminals()} == {"q1", "q3"}
        assert len(tree.get_terminals()) == 2
    assert len(calls) == 5
