"""全テストで共有する fixture と、外部 API を遮断する補助関数。"""

import csv
import importlib.util
import os
import shutil
import sys
import types
from pathlib import Path

import pytest
import requests
from Bio import Entrez


PROJECT_ROOT = Path(__file__).parents[1]
PHYLO_SCRIPT = PROJECT_ROOT / "app" / "Phylo-MIP.py"
FIXTURES = PROJECT_ROOT / "tests" / "fixtures"


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
    """Load the script after replacing its command-line environment."""
    ete3_stub = types.ModuleType("ete3")
    ete3_stub.NodeStyle = object
    monkeypatch.setitem(sys.modules, "ete3", ete3_stub)
    monkeypatch.setattr(sys, "argv", [str(PHYLO_SCRIPT), str(input_csv), *arguments])
    module_name = f"phylo_test_{id(input_csv)}"
    spec = importlib.util.spec_from_file_location(module_name, PHYLO_SCRIPT)
    module = importlib.util.module_from_spec(spec)
    monkeypatch.setitem(sys.modules, module_name, module)
    spec.loader.exec_module(module)
    return module.__dict__


def write_blast_csv(path, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["qseqid", "sallacc", "pident", "qseq"])
        writer.writeheader()
        writer.writerows(rows)
