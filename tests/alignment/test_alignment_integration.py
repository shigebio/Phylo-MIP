from pathlib import Path

import pandas as pd
import pytest
from Bio import SeqIO

from conftest import FIXTURES

pytestmark = pytest.mark.integration


def test_real_vsearch_collapses_only_identical_sequences(phylo_module, require_tools):
    require_tools("vsearch")
    source = FIXTURES / "alignment/sequences.fasta"
    alignment = Path(phylo_module["alignment_dir"])
    result = phylo_module["run_vsearch"](source, "clustered.fasta", str(alignment))
    inputs = list(SeqIO.parse(source, "fasta"))
    outputs = list(SeqIO.parse(result, "fasta"))
    assert len(inputs) == 4 and len(outputs) == 3
    assert {str(r.seq) for r in outputs} == {str(r.seq) for r in inputs}
    assert len({r.id for r in outputs}) == 3
    assert {r.id for r in outputs} <= {r.id for r in inputs}
    tsv = pd.read_csv(alignment / f"{phylo_module['timestamp']}_haplotype_clusters.tsv", sep="\t")
    converted = pd.read_csv(alignment / f"{phylo_module['timestamp']}_haplotype_clusters.csv")
    pd.testing.assert_frame_equal(tsv, converted)
    assert list(tsv.columns[:3]) == ["label", "Group", "numOtus"]
    assert len(tsv) == 1
    assert int(tsv.iloc[0, 2]) == 3
    assert sorted(tsv.iloc[0, 3:].astype(int).tolist()) == [1, 1, 2]


def test_real_mafft_preserves_otus_and_ungapped_sequences(phylo_module, require_tools):
    require_tools("mafft")
    source = FIXTURES / "alignment/sequences.fasta"
    result = phylo_module["run_mafft"](source, "aligned.fasta")
    inputs = {r.id: str(r.seq).upper() for r in SeqIO.parse(source, "fasta")}
    outputs = list(SeqIO.parse(result, "fasta"))
    assert len(outputs) == len(inputs)
    assert {r.id: str(r.seq).replace("-", "").upper() for r in outputs} == inputs
    assert len({len(r.seq) for r in outputs}) == 1
    assert Path(result).name == "aligned.fasta"
