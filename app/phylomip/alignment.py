"""Sequence clustering and multiple sequence alignment."""

import os
import subprocess
import pandas as pd


def run_vsearch(input_fasta, output_centroids, output_dir, timestamp, subprocess_module=subprocess):
    os.makedirs(output_dir, exist_ok=True)
    haplotype_tsv = os.path.join(output_dir, f"{timestamp}_haplotype_clusters.tsv")
    output_centroids_path = os.path.join(output_dir, output_centroids)
    command = f"vsearch --cluster_fast {input_fasta} -id 1 --centroids {output_centroids_path} --mothur_shared_out {haplotype_tsv}"
    print(f"Running VSEARCH: {command}")
    subprocess_module.run(command, shell=True, check=True)
    haplotype_df = pd.read_csv(haplotype_tsv, sep="\t")
    csv_output_file = os.path.join(output_dir, f"{timestamp}_haplotype_clusters.csv")
    haplotype_df.to_csv(csv_output_file, index=False)
    print(f"Haplotype data saved to {csv_output_file}")
    return output_centroids_path


def run_mafft(vsearch_output, output_aligned, alignment_dir, subprocess_module=subprocess):
    aligned_path = os.path.join(alignment_dir, output_aligned)
    command = f"mafft --auto {vsearch_output} > {aligned_path}"
    print(f"Running MAFFT: {command}")
    subprocess_module.run(command, shell=True, check=True)
    return aligned_path

