from pathlib import Path

from Bio import SeqIO

from conftest import load_phylo_module, write_blast_csv


def load_module(tmp_path, monkeypatch):
    input_path = tmp_path / "input.csv"
    write_blast_csv(input_path, [{"qseqid": "q1", "sallacc": "ACC1", "pident": 99.0, "qseq": "ATGC"}])
    return load_phylo_module(input_path, monkeypatch, "--onlyp")


def test_vsearch_outputs_collapsed_fasta_and_consistent_cluster_tables(tmp_path, monkeypatch):
    # VSEARCH成果物のcollapse結果とCSV/TSV整合を確認する / Verify VSEARCH collapse output and CSV/TSV consistency.
    module = load_module(tmp_path, monkeypatch)
    input_fasta = tmp_path / "input.fasta"
    input_fasta.write_text(">q1\nATGC\n>q2\nATGC\n>q3\nGGTA\n", encoding="utf-8")

    def fake_run(command, shell, check):
        centroid_path = Path(module["alignment_dir"]) / "clustered.fasta"
        tsv_path = Path(module["alignment_dir"]) / f"{module['timestamp']}_haplotype_clusters.tsv"
        centroid_path.write_text(">q1\nATGC\n>q3\nGGTA\n", encoding="utf-8")
        tsv_path.write_text("qseqid\thaplotype\nq1\t1\nq2\t1\nq3\t2\n", encoding="utf-8")

    monkeypatch.setattr(module["subprocess"], "run", fake_run)
    result = module["run_vsearch"](input_fasta, "clustered.fasta", module["alignment_dir"])

    assert Path(result).name == "clustered.fasta"
    assert [record.id for record in SeqIO.parse(result, "fasta")] == ["q1", "q3"]
    tsv = Path(module["alignment_dir"]) / f"{module['timestamp']}_haplotype_clusters.tsv"
    csv_path = Path(module["alignment_dir"]) / f"{module['timestamp']}_haplotype_clusters.csv"
    assert tsv.exists() and tsv.stat().st_size > 0
    assert csv_path.exists() and csv_path.stat().st_size > 0
    assert csv_path.read_text(encoding="utf-8").splitlines() == [
        "qseqid,haplotype", "q1,1", "q2,1", "q3,2"
    ]


def test_mafft_outputs_parseable_equal_length_alignment(tmp_path, monkeypatch):
    # MAFFT aligned FASTAのparse、OTU保持、alignment長を確認する / Verify parseability, OTU retention, and alignment length.
    module = load_module(tmp_path, monkeypatch)
    input_fasta = tmp_path / "clustered.fasta"
    input_fasta.write_text(">q1\nATGC\n>q3\nGGTA\n", encoding="utf-8")

    def fake_run(command, shell, check):
        aligned_path = Path(module["alignment_dir"]) / "aligned.fasta"
        aligned_path.write_text(">q1\nATGC-\n>q3\nGGTAA\n", encoding="utf-8")

    monkeypatch.setattr(module["subprocess"], "run", fake_run)
    result = module["run_mafft"](input_fasta, "aligned.fasta")
    records = list(SeqIO.parse(result, "fasta"))

    assert Path(result).name == "aligned.fasta"
    assert [record.id for record in records] == ["q1", "q3"]
    assert len({len(record.seq) for record in records}) == 1


def test_fasttree_newick_and_nexus_preserve_terminal_taxa(tmp_path, monkeypatch):
    # FastTree Newick/Nexusのparseとtaxa集合を確認する / Verify FastTree Newick/Nexus parsing and terminal taxa.
    module = load_module(tmp_path, monkeypatch)
    aligned = tmp_path / "aligned.fasta"
    aligned.write_text(">q1\nATGC\n>q3\nGGTA\n", encoding="utf-8")

    def fake_run(command, shell, check):
        tree_path = Path(module["phylogeny_dir"]) / f"{module['timestamp']}_NJ_tree.nwk"
        tree_path.write_text("(q1:0.1,q3:0.2);\n", encoding="utf-8")

    monkeypatch.setattr(module["subprocess"], "run", fake_run)
    newick_path = module["run_fasttree"](aligned, "tree.nwk", method="NJ", bootstrap=0)

    from Bio import Phylo

    tree = Phylo.read(newick_path, "newick")
    nexus_path = tmp_path / "tree.nex"
    Phylo.write(tree, nexus_path, "nexus")
    nexus_tree = Phylo.read(nexus_path, "nexus")
    assert {terminal.name for terminal in tree.get_terminals()} == {"q1", "q3"}
    assert len(nexus_tree.get_terminals()) == 2