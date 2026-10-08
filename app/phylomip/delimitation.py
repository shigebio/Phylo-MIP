"""Shared bPTP/mPTP partition parsing and taxonomy CSV updates."""

import glob
import os
import re
import pandas as pd


def extract_species_data(file_path):
    species_data, support_data = {}, {}
    current_species = current_support = None
    try:
        with open(file_path, "r") as handle:
            lines = handle.readlines()
        i = 0
        while i < len(lines):
            line = lines[i].strip()
            match = re.search(r"Species\s+(\d+)", line)
            if match:
                current_species = match.group(1)
                support_match = re.search(r"\(support\s+=\s+([\d\.]+)\)", line)
                current_support = support_match.group(1) if support_match else None
                i += 1
                if i < len(lines):
                    next_line = lines[i].strip()
                    if "," in next_line:
                        ids = [item.strip() for item in next_line.split(",") if item.strip()]
                        for seq_id in ids:
                            qseqid = seq_id.split("_")[0] if "_" in seq_id else seq_id
                            species_data[qseqid] = current_species; support_data[qseqid] = current_support
                    elif "_" in next_line and not re.search(r"Species\s+\d+", next_line):
                        qseqid = next_line.split("_")[0]
                        species_data[qseqid] = current_species; support_data[qseqid] = current_support
                        j = i + 1
                        while j < len(lines):
                            candidate = lines[j].strip()
                            if not re.search(r"Species\s+\d+", candidate) and candidate:
                                qseqid = candidate.split("_")[0]
                                species_data[qseqid] = current_species; support_data[qseqid] = current_support
                                j += 1
                            else:
                                break
                        i = j - 1
            i += 1
    except Exception as exc:
        print(f"Error processing file {file_path}: {exc}")
    return species_data, support_data


def update_csv_with_species_data(csv_path, bptp_bayes_data, bptp_bayes_support,
                                 bptp_ml_data, bptp_ml_support, mptp_data):
    try:
        df = pd.read_csv(csv_path)
        columns = {
            "bPTP_Bayesian_Partitioned_Species": None,
            "PTPhSupport_support": None,
            "bPTP_ML_Partitioned_Species": None,
            "PTPML_support": None,
            "mPTP_Partitioned_Species": None,
        }
        for name, default in columns.items():
            if name not in df.columns:
                df[name] = default
        for index, row in df.iterrows():
            qseqid = str(row["qseqid"])
            if qseqid in bptp_bayes_data:
                df.at[index, "bPTP_Bayesian_Partitioned_Species"] = bptp_bayes_data[qseqid]
                df.at[index, "PTPhSupport_support"] = bptp_bayes_support.get(qseqid)
            if qseqid in bptp_ml_data:
                df.at[index, "bPTP_ML_Partitioned_Species"] = bptp_ml_data[qseqid]
                df.at[index, "PTPML_support"] = bptp_ml_support.get(qseqid)
            if qseqid in mptp_data:
                df.at[index, "mPTP_Partitioned_Species"] = mptp_data[qseqid]
        df.to_csv(csv_path, index=False)
        print(f"Updated {csv_path} with PTP species information and support values")
    except Exception as exc:
        print(f"Error updating CSV file {csv_path}: {exc}")


def process_ptp_outputs(bptp_base_dir, mptp_base_dir, taxonomy_dir, output_dir, timestamp):
    print("Processing PTP output files...")
    bptp_dirs = glob.glob(os.path.join(bptp_base_dir, f"{timestamp}_bPTP_analysis*"))
    mptp_dirs = glob.glob(os.path.join(mptp_base_dir, f"{timestamp}_mPTP_analysis*"))
    csv_files = glob.glob(os.path.join(taxonomy_dir, "*taxonomic_data.csv"))
    if not csv_files:
        csv_files = glob.glob(os.path.join(taxonomy_dir, "*.csv"))
    if not csv_files:
        csv_files = glob.glob(os.path.join(output_dir, "*.csv"))
    bptp_bayes_data, bptp_bayes_support, bptp_ml_data, bptp_ml_support = {}, {}, {}, {}
    for directory in bptp_dirs:
        for path in glob.glob(os.path.join(directory, "*.PTPhSupportPartition.txt")):
            print(f"Processing bPTP Bayesian file: {path}")
            data, support = extract_species_data(path); bptp_bayes_data.update(data); bptp_bayes_support.update(support)
        for path in glob.glob(os.path.join(directory, "*.PTPMLPartition.txt")):
            print(f"Processing bPTP ML file: {path}")
            data, support = extract_species_data(path); bptp_ml_data.update(data); bptp_ml_support.update(support)
    mptp_data = {}
    for directory in mptp_dirs:
        for path in glob.glob(os.path.join(directory, "*.txt")):
            print(f"Processing mPTP file: {path}")
            data, _ = extract_species_data(path); mptp_data.update(data)
    for path in csv_files:
        print(f"Updating CSV file: {path}")
        update_csv_with_species_data(path, bptp_bayes_data, bptp_bayes_support,
                                      bptp_ml_data, bptp_ml_support, mptp_data)
    print("PTP output processing complete!")

