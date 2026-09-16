import csv
import importlib.util
import sys
import types
from pathlib import Path


PROJECT_ROOT = Path(__file__).parents[1]
PHYLO_SCRIPT = PROJECT_ROOT / "app" / "Phylo-MIP.py"


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