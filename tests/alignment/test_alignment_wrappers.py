"""VSEARCH/MAFFT wrapper のコマンド、出力ファイル名、変換結果を mock で確認する。"""

import shlex
from pathlib import Path

import pandas as pd
from Bio import SeqIO

from conftest import FIXTURES


def test_vsearch_command_and_cluster_table_conversion(phylo_module, tmp_path, monkeypatch):
    input_fasta = tmp_path / "input.fasta"
    input_fasta.write_text(">q1\nATGC\n>q2\nATGC\n>q3\nGGTA\n")
    alignment_dir = Path(phylo_module["alignment_dir"])
    timestamp = phylo_module["timestamp"]
    output_name = f"{timestamp}_clustered_sequences.fasta"
    calls = []

    def fake_run(command, shell, check):
        assert shell and check
        tokens = shlex.split(command)
        assert tokens == ["vsearch", "--cluster_fast", str(input_fasta), "-id", "1",
                          "--centroids", str(alignment_dir / output_name), "--mothur_shared_out",
                          str(alignment_dir / f"{timestamp}_haplotype_clusters.tsv")]
        calls.append(command)
        Path(tokens[tokens.index("--centroids") + 1]).write_bytes(
            (FIXTURES / "alignment/clustered.fasta").read_bytes())
        Path(tokens[tokens.index("--mothur_shared_out") + 1]).write_bytes(
            (FIXTURES / "alignment/haplotype_clusters.tsv").read_bytes())

    monkeypatch.setattr(phylo_module["subprocess"], "run", fake_run)
    result = phylo_module["run_vsearch"](input_fasta, output_name, str(alignment_dir))
    assert len(calls) == 1
    assert Path(result) == alignment_dir / output_name
    assert {r.id: str(r.seq) for r in SeqIO.parse(result, "fasta")} == {"q1": "ATGC", "q3": "GGTA"}
    tsv = pd.read_csv(alignment_dir / f"{timestamp}_haplotype_clusters.tsv", sep="\t")
    converted = pd.read_csv(alignment_dir / f"{timestamp}_haplotype_clusters.csv")
    pd.testing.assert_frame_equal(tsv, converted)
    assert converted.to_dict("list") == {
        "label": ["unique"], "Group": ["all"], "numOtus": [2], "Otu1": [2], "Otu2": [1]}


def test_mafft_command_and_alignment_output(phylo_module, tmp_path, monkeypatch):
    input_fasta = tmp_path / "clustered.fasta"
    input_fasta.write_text(">q1\nATGC\n>q3\nATGCA\n")
    output_name = f"{phylo_module['timestamp']}_aligned_sequences.fasta"
    expected_path = Path(phylo_module["alignment_dir"]) / output_name
    calls = []

    def fake_run(command, shell, check):
        assert shell and check
        tokens = shlex.split(command)
        assert tokens == ["mafft", "--auto", str(input_fasta), ">", str(expected_path)]
        calls.append(command)
        Path(tokens[-1]).write_text(">q1\nATGC-\n>q3\nATGCA\n")

    monkeypatch.setattr(phylo_module["subprocess"], "run", fake_run)
    result = phylo_module["run_mafft"](input_fasta, output_name)
    assert len(calls) == 1
    assert Path(result) == expected_path
    records = list(SeqIO.parse(result, "fasta"))
    assert {r.id: str(r.seq).replace("-", "") for r in records} == {"q1": "ATGC", "q3": "ATGCA"}
    assert len(records) == 2
    assert {len(r.seq) for r in records} == {5}
