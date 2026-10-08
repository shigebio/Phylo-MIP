"""Execution orchestration for the Phylo-MIP analysis stages."""

import glob
import os
import re
import pandas as pd

from .alignment import run_mafft, run_vsearch
from .bptp import run_bptp
from .cli import parse_args
from .delimitation import process_ptp_outputs
from .mptp import run_mptp
from .output import (OutputPaths, copy_output_files_with_structure,
                     create_output_directories, verify_directory_structure)
from .phylogeny import convert_newick_to_nexus, run_fasttree
from .taxonomy import TaxonomyProcessor


INVALID_CHARS = r'[<>:"/\\|?*]'


def run(args):
    input_file_path = os.path.abspath(args.input_csv)
    paths = OutputPaths.from_input(input_file_path, getattr(args, "_timestamp", None))
    print(f"Output directory will be created at: {paths.output_dir}")
    create_output_directories(paths)
    verify_directory_structure(paths)
    try:
        df = pd.read_csv(input_file_path)
        print(f"File loaded successfully from: {input_file_path}")
    except Exception as exc:
        print(f"Error reading CSV file: {exc}")
        raise SystemExit(1)
    df = df.sort_values("pident", ascending=False).groupby("qseqid").head(args.top)
    taxonomy = TaxonomyProcessor(paths.taxonomy_dir)
    sanitized_input_name = re.sub(INVALID_CHARS, "_", os.path.basename(args.input_csv))
    output_base = args.output_base if args.output_base else f"{sanitized_input_name}_output"
    if args.onlyp:
        input_csv = os.path.join(input_file_path)
        filtered_filename = os.path.join(paths.taxonomy_dir, f"{output_base}_filtered.csv")
        if args.class_name:
            taxonomy.filter_by_class(input_csv, filtered_filename, args.class_name)
        fasta_filename = os.path.join(paths.taxonomy_dir, f"{output_base}.fasta")
        main_fasta = taxonomy.csv_to_fasta(input_csv, fasta_filename)
    else:
        main_fasta = taxonomy.process_with_progress(df, args.class_name)
    if args.tree:
        vsearch_output_name = f"{paths.timestamp}_clustered_sequences.fasta"
        aligned_fasta_name = f"{paths.timestamp}_aligned_sequences.fasta"
        tree_file_name = f"{paths.timestamp}_{args.method}_phylogenetic_tree.nwk"
        vsearch_output_path = run_vsearch(main_fasta, vsearch_output_name, paths.alignment_dir, paths.timestamp)
        aligned_fasta_path = run_mafft(vsearch_output_path, aligned_fasta_name, paths.alignment_dir)
        tree_file_path = run_fasttree(aligned_fasta_path, tree_file_name, paths.phylogeny_dir, paths.timestamp,
                                      method=args.method, bootstrap=args.bootstrap,
                                      gamma=args.gamma, outgroup=args.outgroup)
        nexus_path = os.path.join(paths.phylogeny_dir, f"{paths.timestamp}_{args.method}_phylogenetic_tree.nex")
        convert_newick_to_nexus(tree_file_path, nexus_path)
        print(f"MCMC: {args.mcmc}, Thinning: {args.thinning}, Burn-in: {args.burnin}, Seed: {args.seed}")
        run_bptp(tree_file_path, args.mcmc, args.thinning, args.burnin, args.seed, paths.bptp_base_dir, paths.timestamp)
        run_mptp(tree_file_path, paths.mptp_base_dir, paths.timestamp)
        process_ptp_outputs(paths.bptp_base_dir, paths.mptp_base_dir, paths.taxonomy_dir,
                            paths.output_dir, paths.timestamp)
    print("\n=== Final Output Summary ===")
    verify_directory_structure(paths)
    print("Script execution completed.")
    print("Copying output files with directory structure...")
    final_output_path = copy_output_files_with_structure(paths.output_dir, paths.input_dir)
    print(f"Process complete. Output files saved to: {final_output_path}")
    print("Checking output files location...")
    if os.path.exists(paths.output_dir):
        print(f"Output directory exists at: {paths.output_dir}")
        final_output_path = copy_output_files_with_structure(paths.output_dir, paths.input_dir)
        if final_output_path:
            print(f"Process complete. Output files saved to: {final_output_path}")
        else:
            print(f"Using original output directory: {paths.output_dir}")
    else:
        print(f"Error: Output directory not found at: {paths.output_dir}")
        alt_outputs = glob.glob(os.path.join(paths.input_dir, "phylomip_output_*"))
        if alt_outputs:
            print(f"Found alternative output directories: {alt_outputs}")
            print(f"Using alternative output: {alt_outputs[-1]}")
        else:
            print("No output directories found")
    print("Script execution completed.")
    return paths


def main(argv=None):
    return run(parse_args(argv))

