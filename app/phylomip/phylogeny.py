"""Phylogenetic tree construction and format conversion."""

import os
import subprocess
from Bio import Phylo


def run_fasttree(input_aligned, output_tree, phylogeny_dir, timestamp, method="NJ", bootstrap=1000,
                 gamma=False, outgroup=None, subprocess_module=subprocess):
    tree_path = os.path.join(phylogeny_dir, f"{timestamp}_{method}_{output_tree}")
    command = "fasttree -nt"
    if method == "NJ":
        command += " -nj"
    elif method != "ML":
        raise ValueError(f"Invalid method '{method}'. Choose 'ML' or 'NJ'.")
    if bootstrap > 0:
        command += f" -boot {bootstrap}"
    if gamma:
        command += " -gamma"
    if outgroup:
        command += f" -outgroup {outgroup}"
    command += f" {input_aligned} > {tree_path}"
    print(f"Running FastTree: {command}")
    subprocess_module.run(command, shell=True, check=True)
    print(f"Phylogenetic tree saved to {tree_path}")
    return tree_path


def convert_newick_to_nexus(tree_path, nexus_path):
    tree = Phylo.read(tree_path, "newick")
    Phylo.write(tree, nexus_path, "nexus")
    print("Newick file from FastTree has been converted to NEXUS format.")
    return nexus_path

