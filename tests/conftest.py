"""全テストで共有する fixture と、外部 API を遮断する補助関数。"""

import csv
import os
import shutil
import sys
from pathlib import Path

import pytest
import requests
from Bio import Entrez


PROJECT_ROOT = Path(__file__).parents[1]
PHYLO_SCRIPT = PROJECT_ROOT / "app" / "Phylo-MIP.py"
FIXTURES = PROJECT_ROOT / "tests" / "fixtures"
APP_ROOT = PROJECT_ROOT / "app"
sys.path.insert(0, str(APP_ROOT))

from phylomip.alignment import run_mafft as _run_mafft, run_vsearch as _run_vsearch
from phylomip.bptp import run_bptp as _run_bptp
from phylomip.delimitation import extract_species_data, process_ptp_outputs as _process_ptp_outputs
from phylomip.mptp import run_mptp as _run_mptp
from phylomip.output import OutputPaths, create_output_directories
from phylomip.phylogeny import run_fasttree as _run_fasttree
from phylomip.taxonomy import TaxonomyProcessor
from phylomip.pipeline import run as run_pipeline
from phylomip.cli import parse_args
import phylomip.taxonomy as taxonomy_module


@pytest.fixture
def require_tools():
    def require(*names):
        if os.environ.get("RUN_INTEGRATION") != "1":
            pytest.skip("Set RUN_INTEGRATION=1 to execute external analysis tools")
        missing = [name for name in names if shutil.which(name) is None]
        if missing:
            pytest.skip("Missing external tools: " + ", ".join(missing))
    return require


@pytest.fixture(autouse=True)
def block_live_taxonomy(monkeypatch):
    def unexpected_request(*args, **kwargs):
        pytest.fail("Live taxonomy requests are forbidden; supply a mock response")

    monkeypatch.setattr(requests.sessions.Session, "request", unexpected_request)
    monkeypatch.setattr(Entrez, "efetch", unexpected_request)


@pytest.fixture
def phylo_module(tmp_path, monkeypatch):
    input_path = tmp_path / "input.csv"
    write_blast_csv(input_path, [
        {"qseqid": "q1", "sallacc": "ACC1", "pident": 99.0, "qseq": "ATGC"}
    ])
    return load_phylo_module(input_path, monkeypatch, "--onlyp")


def load_phylo_module(input_csv, monkeypatch, *arguments):
    """Build a compatibility context while tests migrate to normal imports."""
    import pandas as pd
    paths = OutputPaths.from_input(input_csv)
    create_output_directories(paths)
    processor = TaxonomyProcessor(paths.taxonomy_dir)
    context = {
        "input_file_path": str(input_csv), "timestamp": paths.timestamp,
        "alignment_dir": paths.alignment_dir, "phylogeny_dir": paths.phylogeny_dir,
        "bptp_base_dir": paths.bptp_base_dir, "mptp_base_dir": paths.mptp_base_dir,
        "taxonomy_dir": paths.taxonomy_dir, "output_dir": paths.output_dir,
        "df": pd.read_csv(input_csv), "subprocess": __import__("subprocess"),
        "requests": taxonomy_module.requests, "Entrez": taxonomy_module.Entrez,
        "time": taxonomy_module.time,
    }
    def process_row(index, row):
        processor.fetch_ncbi_data = context["fetch_ncbi_data"]
        processor.get_gbif_taxonomic_info = context["get_gbif_taxonomic_info"]
        return processor.process_row(index, row)
    def process_with_progress(filter_class=None):
        processor.fetch_ncbi_data = context["fetch_ncbi_data"]
        processor.get_gbif_taxonomic_info = context["get_gbif_taxonomic_info"]
        return processor.process_with_progress(context["df"], filter_class)
    context.update({
        "fetch_ncbi_data": processor.fetch_ncbi_data,
        "get_gbif_taxonomic_info": processor.get_gbif_taxonomic_info,
        "process_row": process_row,
        "process_with_progress": process_with_progress,
        "save_csv": processor.save_csv, "save_fasta": processor.save_fasta,
        "run_vsearch": lambda input_fasta, output_centroids, output_dir: _run_vsearch(
            input_fasta, output_centroids, output_dir, paths.timestamp, context["subprocess"]),
        "run_mafft": lambda vsearch_output, output_aligned: _run_mafft(
            vsearch_output, output_aligned, paths.alignment_dir, context["subprocess"]),
        "run_fasttree": lambda input_aligned, output_tree, method="NJ", bootstrap=1000, gamma=False, outgroup=None: _run_fasttree(
            input_aligned, output_tree, paths.phylogeny_dir, paths.timestamp, method, bootstrap, gamma, outgroup, context["subprocess"]),
        "run_bptp": lambda tree_file, mcmc, thinning, burnin, seed, base_dir: _run_bptp(
            tree_file, mcmc, thinning, burnin, seed, base_dir, paths.timestamp, context["subprocess"]),
        "run_mptp": lambda tree_file, base_dir: _run_mptp(tree_file, base_dir, paths.timestamp, context["subprocess"]),
        "extract_species_data": extract_species_data,
        "process_ptp_outputs": lambda: _process_ptp_outputs(paths.bptp_base_dir, paths.mptp_base_dir,
                                                              paths.taxonomy_dir, paths.output_dir, paths.timestamp),
    })
    if "--tree" in arguments:
        pipeline_args = parse_args([str(input_csv), *arguments])
        pipeline_args._timestamp = paths.timestamp
        run_pipeline(pipeline_args)
    return context


def write_blast_csv(path, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["qseqid", "sallacc", "pident", "qseq"])
        writer.writeheader()
        writer.writerows(rows)
