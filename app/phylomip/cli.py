"""Command-line argument definition and validation."""

import argparse
import random


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Process Sequence file and generate phylogenetic trees")
    parser.add_argument("input_csv", help="Input CSV file path")
    parser.add_argument("--o", type=str, dest="output_base", help="Output base name")
    parser.add_argument("--top", type=int, default=1, choices=range(1, 10),
                        help="Number of top results to retain per qseqid (default: 1, range: 1-10)")
    parser.add_argument("--onlyp", action="store_true", help="Run only phylogenic analysis")
    parser.add_argument("--class", type=str, dest="class_name", nargs="+",
                        help="Select using class for phylogenetic analysis")
    parser.add_argument("--tree", action="store_true", help="Generate phylogenetic tree")
    tree_group = parser.add_argument_group("FastTree options", "Options for FastTree analysis")
    tree_group.add_argument("--method", default="NJ", choices=["NJ", "ML"],
                            help="Tree generation method: NJ or ML")
    tree_group.add_argument("--bootstrap", type=int, default=0, help="Number of bootstrap replicates")
    tree_group.add_argument("--gamma", action="store_true", help="Use gamma model")
    tree_group.add_argument("--outgroup", help="Outgroup for tree")
    parser.add_argument("--bptp", action="store_true", help="Enable bPTP options")
    bptp_group = parser.add_argument_group("bPTP options", "Options for bPTP analysis")
    bptp_group.add_argument("--mcmc", type=int, default=100000, help="Number of MCMC iterations (default: 100000)")
    bptp_group.add_argument("--thinning", type=int, default=100, help="Thinning value (default: 100)")
    bptp_group.add_argument("--burnin", type=float, default=0.1, help="Burn-in fraction (default: 0.1)")
    bptp_group.add_argument("--seed", type=int, default=random.randint(1, 10**6),
                            help="Random seed")
    return parser


def parse_args(argv=None):
    """Parse CLI options without executing analysis."""
    return build_parser().parse_args(argv)

