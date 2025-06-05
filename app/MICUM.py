import os
import shutil
import stat
from datetime import datetime
import pandas as pd
from Bio import Entrez, Phylo
import requests
import sys
import time
import random
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
import argparse
import re
import six
import csv
import glob
from ete3 import NodeStyle
sys.path.append('/usr/local/lib/python3.7/site-packages')

# If you want to receive notifications from NCBI, change it to your own email address
Entrez.email = "your_email@example.com"
default_seed = random.randint(1, 10**6)

# Parsing optional arguments
parser = argparse.ArgumentParser(description="Process Sequence file and generate phylogenetic trees")
parser.add_argument('input_csv', help="Input CSV file path")
parser.add_argument('--o',  type=str, dest="output_base", help="Output base name")
parser.add_argument('--top', type=int, default=1, choices=range(1, 10), help="Number of top results to retain per qseqid (default: 1, range: 1-10)")
parser.add_argument('--onlyp', action='store_true', help="Run only phylogenic analysis")
parser.add_argument('--class', type=str, dest="class_name", nargs='+', help="Select using class for phylogenetic analysis")
parser.add_argument('--tree', help="Generate phylogenetic tree", action='store_true')
tree_group = parser.add_argument_group('FastTree options', 'Options for FastTree analysis')
tree_group.add_argument('--method', default="NJ", choices=["NJ", "ML"], help="Tree generation method: NJ or ML")
tree_group.add_argument('--bootstrap', type=int, default=0, help="Number of bootstrap replicates")
tree_group.add_argument('--gamma', help="Use gamma model", action='store_true')
tree_group.add_argument('--outgroup', help="Outgroup for tree")
parser.add_argument('--bptp', help='Enable bPTP options', action='store_true')
bptp_group = parser.add_argument_group('bPTP options', 'Options for bPTP analysis')
bptp_group.add_argument('--mcmc', type=int, default=100000, help='Number of MCMC iterations (default: 100000)')
bptp_group.add_argument('--thinning', type=int, default=100, help='Thinning value (default: 100)')
bptp_group.add_argument('--burnin', type=float, default=0.1, help='Burn-in fraction (default: 0.1)')
bptp_group.add_argument('--seed', type=int, default=default_seed, help=f'Random seed (default: {default_seed})')

args = parser.parse_args()
input_file_path = os.path.abspath(args.input_csv)
input_dir = os.path.dirname(input_file_path)
timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
output_dir = os.path.join(input_dir, f"micum_output_{timestamp}")
print(f"Output directory will be created at: {output_dir}")
os.makedirs(output_dir, exist_ok=True)

taxonomy_dir = os.path.join(output_dir, "taxonomy")
alignment_dir = os.path.join(output_dir, "alignment")
phylogeny_dir = os.path.join(output_dir, "phylogeny")
bptp_base_dir = os.path.join(output_dir, "bptp")
mptp_base_dir = os.path.join(output_dir, "mptp")

def copy_with_directory_structure(src_dir, dest_dir):
    """
    A function to copy files while preserving the directory structure
    """
    try:
        # If the directory exists, delete it
        if os.path.exists(dest_dir):
            shutil.rmtree(dest_dir)

        shutil.copytree(src_dir, dest_dir)
        print(f"Successfully copied directory structure from {src_dir} to {dest_dir}")

    except Exception as e:
        print(f"Error copying directory structure: {e}")

def copy_output_files_with_structure(output_dir, original_dir):
    """
    A function that copies the output directory to the original directory while preserving the directory structure.
    """
    output_dirname = os.path.basename(output_dir)
    dest_path = os.path.join(original_dir, f"copied_{output_dirname}")

    try:
        if os.path.exists(dest_path):
            shutil.rmtree(dest_path)
        shutil.copytree(output_dir, dest_path)
        print(f"Successfully copied directory to: {dest_path}")
        return dest_path
    except Exception as e:
        print(f"Error copying output files: {e}")
        return None

# Create subdirectories
def create_output_directories():
    """
    Create the output directory and set the correct permissions
    """
    directories = [
        output_dir,
        taxonomy_dir,
        alignment_dir,
        phylogeny_dir,
        bptp_base_dir,
        mptp_base_dir
    ]

    for directory in directories:
        try:
            os.makedirs(directory, exist_ok=True)
            try:
                os.chmod(directory, stat.S_IRWXU | stat.S_IRWXG | stat.S_IRWXO)
            except (OSError, PermissionError) as e:
                print(f"Warning: Could not set permissions for {directory}: {e}")

            print(f"Created directory: {directory}")

            if not os.path.exists(directory):
                print(f"ERROR: Failed to create directory: {directory}")
            else:
                print(f"SUCCESS: Directory exists: {directory}")

        except Exception as e:
            print(f"Error creating directory {directory}: {e}")
            sys.exit(1)

create_output_directories()

def verify_directory_structure(): # debug code
    """
    Check the directory structure and file location
    """
    print("\n=== Directory Structure Verification ===")
    for root, dirs, files in os.walk(output_dir):
        level = root.replace(output_dir, '').count(os.sep)
        indent = '  ' * level
        print(f"{indent}{os.path.basename(root)}/")
        subindent = '  ' * (level + 1)
        for file in files:
            print(f"{subindent}{file}")
    print("=========================================\n")

verify_directory_structure()

def ensure_directory_exists(file_path):
    """
    Verify that a directory in a file path exists
    """
    directory = os.path.dirname(file_path)
    if not os.path.exists(directory):
        try:
            os.makedirs(directory, exist_ok=True)
            try:
                os.chmod(directory, stat.S_IRWXU | stat.S_IRWXG | stat.S_IRWXO)
            except (OSError, PermissionError):
                pass  # Attempt to create file even if permission setting fails
        except Exception as e:
            print(f"Warning: Could not create directory {directory}: {e}")

# For filename sanitization
invalid_chars = r'[<>:"/\\|?*]'

# Import CSV
try:
    df = pd.read_csv(input_file_path)
    print(f"File loaded successfully from: {input_file_path}")
except Exception as e:
    print(f"Error reading CSV file: {e}")
    sys.exit(1)

# If there are duplicate qseqids, only the number selected from the one with the highest pindent is left
df = df.sort_values('pident', ascending=False).groupby('qseqid').head(args.top)

fasta_list = []
csv_data = []
cache = {}

# Convert strings that are inappropriate for OTUs to hyphens
def sanitize_otu_name(name):
    name = re.sub(r'\.', '_', name)
    return re.sub(r'[^a-zA-Z0-9_>\-]', '_', name)

# Function to retrieve taxon data from NCBI
def fetch_ncbi_data(accessionID):
    if accessionID in cache:
        return cache[accessionID]

    for attempt in range(3):
        try:
            handle = Entrez.efetch(db="nucleotide", id=accessionID, rettype="gb", retmode="xml")
            records = Entrez.read(handle)
            handle.close()
            cache[accessionID] = records
            return records
        except Exception as e:
            print(f"NCBI API request failed: {e}. Retrying... ({attempt + 1}/3)")
            time.sleep(2 ** attempt)

    # On failure, instead of returning an error message, return None
    print(f"Failed to fetch data for accession ID: {accessionID}")
    return None

# Function to retrieve taxon data from GBIF
def get_gbif_taxonomic_info(species_name):
    url = f"https://api.gbif.org/v1/species/match?name={species_name}"
    for attempt in range(3):
        try:
            response = requests.get(url, timeout=30)
            if response.status_code == 200:
                data = response.json()
                if 'order' in data and data['order'] is not None:  # Ensuring taxonomic info exists
                    return {
                        'source': 'GBIF',
                        'species': data.get('species'),
                        'genus': data.get('genus'),
                        'family': data.get('family'),
                        'order': data.get('order'),
                        'class': data.get('class')
                    }
            print(f"GBIF API request incomplete or failed (attempt {attempt + 1}/3)")
        except requests.exceptions.RequestException as e:
            print(f"GBIF API request failed: {e}. Retrying...")
        time.sleep(1)

    return None

# Narrow the classes by referring to the --class argument
def filter_by_class(input_csv, output_csv, class_name):
    try:
        # Import CSV
        df = pd.read_csv(input_csv)

        # 'class' カラムが存在するか確認
        if 'class' not in df.columns:
            raise KeyError("'class' column not found in the input CSV.")

        # Check if 'class' column exists
        filtered_df = df[df['class'].isin(class_name)]

        # Check that the result is not empty
        if len(filtered_df) == 0:
            print("No matching rows found for the given class.")

        # Output the filtered results
        os.makedirs(os.path.dirname(output_csv), exist_ok=True)  # Create output directory
        filtered_df.to_csv(output_csv, index=False)

        print(f"Filtered data saved to {output_csv}.")
    except FileNotFoundError:
        print(f"Error: Input file {input_csv} not found.")
    except KeyError:
        print("Error: 'Class' column not found in the input CSV.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

def process_row(index, row):
    qseqid = row['qseqid']
    accessionID = row['sallacc']
    pident = row['pident']
    qseq = row['qseq'].replace("-", "N")  # Replace hyphens with N

    try:
        # Obtain species name from NCBI data by referencing accessionID (sallacc)
        records = fetch_ncbi_data(accessionID)
        if not records:
            # If you are unable to obtain information from NCBI
            taxonomic_name = "Uncertain_taxonomy"  # Indeterminate taxon
            fasta_entry = f">{sanitize_otu_name(f'{qseqid}_{accessionID}_{taxonomic_name}_{pident:.2f}')}\n{qseq}\n"
            csv_entry = [qseqid, accessionID, 'Unknown', 'Unknown', taxonomic_name, pident, qseq, 'NCBI Failed']
            return fasta_entry, csv_entry

        organism_name = records[0]['GBSeq_organism']

        # Obtain species name from GBIF data by referencing species name
        gbif_info = get_gbif_taxonomic_info(organism_name)

        if gbif_info:
            taxonomic_info = gbif_info
            source = 'GBIF'
        else:
            # If it fails to obtain from GBIF, obtain from NCBI
            taxonomy_list = records[0]['GBSeq_taxonomy'].split("; ")
            taxonomic_info = {
                'source': 'NCBI',
                'species': organism_name,
                'genus': taxonomy_list[-1] if len(taxonomy_list) > 0 else None,
                'family': taxonomy_list[-2] if len(taxonomy_list) > 1 else None,
                'order': taxonomy_list[-3] if len(taxonomy_list) > 2 else None,
                'class': taxonomy_list[-4] if len(taxonomy_list) > 3 else None
            }
            source = 'NCBI'

        taxonomic_name = "Low_Identity_Match"
        if pident >= 98.00:
            taxonomic_name = taxonomic_info.get('species', 'Uncertain_taxonomy')
        elif 95.00 <= pident < 98.00:
            taxonomic_name = taxonomic_info.get('genus', 'Uncertain_taxonomy')
        elif 90.00 <= pident < 95.00:
            taxonomic_name = taxonomic_info.get('family', 'Uncertain_taxonomy')
        elif 85.00 <= pident < 90.00:
            taxonomic_name = taxonomic_info.get('order', 'Uncertain_taxonomy')

        order = taxonomic_info.get('order', 'Unknown')
        class_name = taxonomic_info.get('class', 'Unknown')
        family_name = taxonomic_info.get('family', 'Unknown')

        # Generates FASTA and CSV output
        fasta_entry = f">{sanitize_otu_name(f'{qseqid}_{accessionID}_{taxonomic_name}_{pident:.2f}')}\n{qseq}\n"
        csv_entry = [qseqid, accessionID, class_name, order, family_name, taxonomic_name, pident, qseq, source]

        return fasta_entry, csv_entry

    except Exception as e:
        print(f"Error in qseqid {qseqid}: {e}")
        # Use "Uncertain_taxonomy" on error
        taxonomic_name = "Uncertain_taxonomy"
        fasta_entry = f">{sanitize_otu_name(f'{qseqid}_{accessionID}_{taxonomic_name}_{pident:.2f}')}\n{qseq}\n"
        csv_entry = [qseqid, accessionID, 'Unknown', 'Unknown', taxonomic_name, pident, qseq, 'Error']
        return fasta_entry, csv_entry

# Function to display the progress of an API request
def process_with_progress(filter_class=None):
    processed_rows = 0
    total_rows = len(df)
    fasta_list = []  # Unfiltered FASTA data
    csv_data = []    # Unfiltered CSV data
    filtered_fasta_list = []  # Filtered FASTA data
    filtered_csv_data = []    # Filtered CSV data

    with ThreadPoolExecutor(max_workers=5) as executor:
        futures = {executor.submit(process_row, index, row): index for index, row in df.iterrows()}
        for future in as_completed(futures):
            index = futures[future]
            try:
                fasta_entry, csv_entry = future.result()
                if fasta_entry and csv_entry:
                    fasta_list.append(fasta_entry)
                    csv_data.append(csv_entry)

                    # Filtering Process
                    if filter_class:
                        class_value = csv_entry[2]  # 'class' column value
                        if class_value is None:
                            class_value = ""
                        if class_value.strip().lower() in [item.lower() for item in filter_class]:
                            filtered_fasta_list.append(fasta_entry)
                            filtered_csv_data.append(csv_entry)
                    else:
                        # If no filter is used, all data is added.
                        filtered_fasta_list.append(fasta_entry)
                        filtered_csv_data.append(csv_entry)

                    processed_rows += 1
                    print(f"\rProcessed row {processed_rows} of {total_rows} ({processed_rows / total_rows:.2%})", end='')
            except Exception as e:
                print(f"\nError processing row {index + 1}: {e}")

    print()  # End progress tracking

    taxonomy_fasta = save_fasta(os.path.join(taxonomy_dir, 'taxonomic_sequences.fasta'), fasta_list)
    save_csv(os.path.join(taxonomy_dir, 'taxonomic_data.csv'), csv_data)

    print("Processing complete.")

    if filter_class:
        filtered_fasta = save_fasta(os.path.join(taxonomy_dir, 'filtered_taxonomic_sequences.fasta'), filtered_fasta_list)
        save_csv(os.path.join(taxonomy_dir, 'filtered_taxonomic_data.csv'), filtered_csv_data)
        return filtered_fasta

    return taxonomy_fasta

def save_fasta(file_path, fasta_entries):
    ensure_directory_exists(file_path)
    try:
        with open(file_path, 'w') as f:
            f.writelines(fasta_entries)
        print(f"SUCCESS: Saved FASTA to {file_path}")
        return file_path
    except Exception as e:
        print(f"ERROR: Failed to save FASTA to {file_path}: {e}")
        raise

def save_csv(file_path, csv_rows):
    ensure_directory_exists(file_path)
    try:
        with open(file_path, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['qseqid', 'accessionID', 'class', 'order', 'taxonomic_name', 'pident', 'qseq', 'source'])
            writer.writerows(csv_rows)
        print(f"SUCCESS: Saved CSV to {file_path}")
    except Exception as e:
        print(f"ERROR: Failed to save CSV to {file_path}: {e}")
        raise

def csv_to_fasta(input_csv, output_fasta): # IMO: It may be more streamlined to combine the conversion to FASTA with the process in process_row
    # Inport CSV
    df = pd.read_csv(input_csv)

    # Check if required columns exist
    if not {'qseqid', 'qseq'}.issubset(df.columns): # Check for the presence of qseq to ensure unique OTUs
        raise ValueError("Input CSV must contain 'qseqid' and 'qseq' columns.")

    with open(output_fasta, 'w') as fasta_file:
        for _, row in df.iterrows():
            # Generate OTUs
            otu_name_parts = [
                str(row[col]).replace(" ", "_").replace(",", "_").replace(".", "_").replace("-", "_")
                for col in df.columns if col != 'qseq' and pd.notna(row[col])
            ]
            otu_name = "_".join(otu_name_parts)

            # Obtain sequence
            sequence = row['qseq']

            # Write a FASTA entry
            fasta_file.write(f">{otu_name}\n{sequence}\n")

    return output_fasta

# VSEARCH
def run_vsearch(input_fasta, output_centroids, output_dir):
    # Specify the output directory
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Run VSEARCH
    haplotype_tsv = os.path.join(alignment_dir, f"{timestamp}_haplotype_clusters.tsv")
    output_centroids_path = os.path.join(alignment_dir, output_centroids)
    vsearch_cmd = f"vsearch --cluster_fast {input_fasta} -id 1 --centroids {output_centroids_path} --mothur_shared_out {haplotype_tsv}"
    print(f"Running VSEARCH: {vsearch_cmd}")
    subprocess.run(vsearch_cmd, shell=True, check=True)

    # Convert TSV to CSV
    haplotype_df = pd.read_csv(haplotype_tsv, sep="\t")
    csv_output_file = os.path.join(alignment_dir, f"{timestamp}_haplotype_clusters.csv")
    haplotype_df.to_csv(csv_output_file, index=False)
    print(f"Haplotype data saved to {csv_output_file}")

    # Return the actual path of the centroids file
    return output_centroids_path


# MAFFT
def run_mafft(vsearch_output, output_aligned):
    aligned_path = os.path.join(alignment_dir, output_aligned)
    mafft_cmd = f"mafft --auto {vsearch_output} > {aligned_path}"
    print(f"Running MAFFT: {mafft_cmd}")
    subprocess.run(mafft_cmd, shell=True, check=True)
    return aligned_path

# FastTree
def run_fasttree(input_aligned, output_tree, method="NJ", bootstrap=1000, gamma=False, outgroup=None):
    tree_path = os.path.join(phylogeny_dir, f"{timestamp}_{method}_{output_tree}")
    fasttree_cmd = "fasttree -nt"

    if method == "NJ":
        fasttree_cmd += " -nj"
    elif method == "ML":
        pass
    else:
        raise ValueError(f"Invalid method '{method}'. Choose 'ML' or 'NJ'.")

    if bootstrap > 0:
        fasttree_cmd += f" -boot {bootstrap}"

    if gamma:
        fasttree_cmd += " -gamma"

    if outgroup:
        fasttree_cmd += f" -outgroup {outgroup}"

    fasttree_cmd += f" {input_aligned} > {tree_path}"
    print(f"Running FastTree: {fasttree_cmd}")
    subprocess.run(fasttree_cmd, shell=True, check=True)
    print(f"Phylogenetic tree saved to {tree_path}")
    return tree_path

# bPTP
def run_bptp(tree_file, mcmc, thinning, burnin, seed, base_dir):
    try:
        bptp_path = os.path.join('/app/PTP/bin', 'bPTP.py')

        # Create timestamped bPTP directory
        bptp_output_dir = os.path.join(base_dir, f'{timestamp}_bPTP_analysis')
        os.makedirs(bptp_output_dir, exist_ok=True)

        # Set the output file name with descriptive prefix
        output_file = os.path.join(bptp_output_dir, f'{timestamp}_bPTP_species_delimitation')

        bptp_command = [
            'xvfb-run', '-a', 'python3', bptp_path,
            '-t', tree_file,
            '-o', output_file,
            '-s', str(seed),
            '-i', str(mcmc),
            '-n', str(thinning),
            '-b', str(burnin)
        ]

        print(f"Running bPTP with command: {' '.join(bptp_command)}")
        subprocess.run(bptp_command, check=True)
        print('bPTP analysis complete.')
    except subprocess.CalledProcessError as e:
        print(f"Error running bPTP.py: {e}")

def run_mptp(tree_file, base_dir):
    try:
        # Create timestamped mPTP directory
        mptp_output_dir = os.path.join(base_dir, f'{timestamp}_mPTP_analysis')
        os.makedirs(mptp_output_dir, exist_ok=True)

        # Set the output file name with descriptive prefix
        output_file = os.path.join(mptp_output_dir, f'{timestamp}_mPTP_species_delimitation')

        subprocess.run(
            ['xvfb-run', '-a', 'mptp', '-tree_file', tree_file,
             '-output_file', output_file, '-ml', '-single'],
            check=True
        )
        print('mPTP analysis complete.')
    except subprocess.CalledProcessError as e:
        print(f"Error running mPTP: {e}")

def extract_species_data(file_path):
    """
    Extract species grouping data and support values from PTP output files.

    Args:
        file_path: Path to the PTP output file

    Returns:
        Tuple of two dictionaries:
        - Dictionary mapping qseqid to species number
        - Dictionary mapping qseqid to support value (may be None for mPTP)
    """
    species_data = {}
    support_data = {}
    current_species = None
    current_support = None

    try:
        with open(file_path, 'r') as f:
            lines = f.readlines()

        i = 0
        while i < len(lines):
            line = lines[i].strip()

            # Search only for species number
            species_match = re.search(r'Species\s+(\d+)', line)
            if species_match:
                current_species = species_match.group(1)

                # Search for support values(for results of bPTP)
                support_match = re.search(r'\(support\s+=\s+([\d\.]+)\)', line)
                current_support = support_match.group(1) if support_match else None

                i += 1

                # Processing patterns with support values
                if i < len(lines):
                    next_line = lines[i].strip()

                    # A pattern where the IDs are separated by commas
                    if ',' in next_line:
                        for seq_id in next_line.split(','):
                            seq_id = seq_id.strip()
                            if seq_id:
                                qseqid = seq_id.split('_')[0] if '_' in seq_id else seq_id
                                species_data[qseqid] = current_species
                                support_data[qseqid] = current_support  # Store support value (may be None)
                    # A pattern with one ID per line
                    elif '_' in next_line and not re.search(r'Species\s+\d+', next_line):
                        seq_id = next_line.strip()
                        qseqid = seq_id.split('_')[0] if '_' in seq_id else seq_id
                        species_data[qseqid] = current_species
                        support_data[qseqid] = current_support  # Store support value (may be None)

                        # Check whether the next line is the same Species ID
                        j = i + 1
                        while j < len(lines):
                            next_line = lines[j].strip()
                            # If the next line is not a Species, it is judged to be a different ID of the same Species.
                            if not re.search(r'Species\s+\d+', next_line) and next_line:
                                seq_id = next_line.strip()
                                qseqid = seq_id.split('_')[0] if '_' in seq_id else seq_id
                                species_data[qseqid] = current_species
                                support_data[qseqid] = current_support  # Store support value (may be None)
                                j += 1
                            else:
                                break
                        i = j - 1  # Move to next Species row or last row
            i += 1

    except Exception as e:
        print(f"Error processing file {file_path}: {e}")

    return species_data, support_data

def update_csv_with_species_data(csv_path, bptp_bayes_data, bptp_bayes_support, bptp_ml_data, bptp_ml_support, mptp_data):
    """
    Update CSV files with species grouping information and support values from PTP analyses.

    Args:
        csv_path: Path to the CSV file to update
        bptp_bayes_data: Dictionary mapping qseqid to bPTP Bayesian species
        bptp_bayes_support: Dictionary mapping qseqid to bPTP Bayesian support values
        bptp_ml_data: Dictionary mapping qseqid to bPTP ML species
        bptp_ml_support: Dictionary mapping qseqid to bPTP ML support values
        mptp_data: Dictionary mapping qseqid to mPTP species
    """
    try:
        # Read the CSV file
        df = pd.read_csv(csv_path)

        # Add new columns if they don't exist
        if 'bPTP_Bayesian_Partitioned_Species' not in df.columns:
            df['bPTP_Bayesian_Partitioned_Species'] = None
        if 'PTPhSupport_support' not in df.columns:
            df['PTPhSupport_support'] = None
        if 'bPTP_ML_Partitioned_Species' not in df.columns:
            df['bPTP_ML_Partitioned_Species'] = None
        if 'PTPML_support' not in df.columns:
            df['PTPML_support'] = None
        if 'mPTP_Partitioned_Species' not in df.columns:
            df['mPTP_Partitioned_Species'] = None

        # Update the values
        for index, row in df.iterrows():
            qseqid = str(row['qseqid'])
            if qseqid in bptp_bayes_data:
                df.at[index, 'bPTP_Bayesian_Partitioned_Species'] = bptp_bayes_data[qseqid]
                df.at[index, 'PTPhSupport_support'] = bptp_bayes_support.get(qseqid)
            if qseqid in bptp_ml_data:
                df.at[index, 'bPTP_ML_Partitioned_Species'] = bptp_ml_data[qseqid]
                df.at[index, 'PTPML_support'] = bptp_ml_support.get(qseqid)
            if qseqid in mptp_data:
                df.at[index, 'mPTP_Partitioned_Species'] = mptp_data[qseqid]

        # Save the updated CSV
        df.to_csv(csv_path, index=False)
        print(f"Updated {csv_path} with PTP species information and support values")
    except Exception as e:
        print(f"Error updating CSV file {csv_path}: {e}")

def process_ptp_outputs():
    """
    Process PTP output files and update CSV files with species grouping information.
    """
    print("Processing PTP output files...")

    # Find bPTP and mPTP directories in their respective base directories
    bptp_dirs = glob.glob(os.path.join(bptp_base_dir, f'{timestamp}_bPTP_analysis*'))
    mptp_dirs = glob.glob(os.path.join(mptp_base_dir, f'{timestamp}_mPTP_analysis*'))

    # Find CSV files in taxonomy directory
    csv_files = glob.glob(os.path.join(taxonomy_dir, '*taxonomic_data.csv'))

    # If no filtered CSV files found, check for other relevant CSVs
    if not csv_files:
        csv_files = glob.glob(os.path.join(taxonomy_dir, '*.csv'))
    if not csv_files:
        csv_files = glob.glob(os.path.join(output_dir, '*.csv'))

    # Process bPTP files
    bptp_bayes_data = {}
    bptp_bayes_support = {}
    bptp_ml_data = {}
    bptp_ml_support = {}

    for bptp_dir in bptp_dirs:
        # Find the Bayesian support partition file
        bayes_files = glob.glob(os.path.join(bptp_dir, '*.PTPhSupportPartition.txt'))
        for bayes_file in bayes_files:
            print(f"Processing bPTP Bayesian file: {bayes_file}")
            bayes_species_data, bayes_support_data = extract_species_data(bayes_file)
            bptp_bayes_data.update(bayes_species_data)
            bptp_bayes_support.update(bayes_support_data)

        # Find the ML partition file
        ml_files = glob.glob(os.path.join(bptp_dir, '*.PTPMLPartition.txt'))
        for ml_file in ml_files:
            print(f"Processing bPTP ML file: {ml_file}")
            ml_species_data, ml_support_data = extract_species_data(ml_file)
            bptp_ml_data.update(ml_species_data)
            bptp_ml_support.update(ml_support_data)

    # Process mPTP files
    # Considering the possibility of extracting some value from the mPTP results. NOTE: Extensibility
    mptp_data = {}
    mptp_support = {}

    for mptp_dir in mptp_dirs:
        # Find all text files
        mptp_files = glob.glob(os.path.join(mptp_dir, '*.txt'))
        for mptp_file in mptp_files:
            print(f"Processing mPTP file: {mptp_file}")
            species_data, support_data = extract_species_data(mptp_file)
            mptp_data.update(species_data)
            mptp_support.update(support_data)

    # Update all CSV files
    for csv_file in csv_files:
        print(f"Updating CSV file: {csv_file}")
        update_csv_with_species_data(csv_file, bptp_bayes_data, bptp_bayes_support,
                                    bptp_ml_data, bptp_ml_support, mptp_data)

    print("PTP output processing complete!")

print("Creating output directories...")
create_output_directories()
print("Directories created successfully.")

sanitized_input_name = re.sub(invalid_chars, '_', os.path.basename(args.input_csv))
output_base = args.output_base if args.output_base else f"{sanitized_input_name}_output"

if args.onlyp: # Branch to execute only phylogenetic analysis and later
    input_csv = os.path.join(input_file_path)
    filtered_filename = os.path.join(taxonomy_dir, f"{output_base}_filtered.csv")

    if args.class_name:  # If --class option is specified, filter by the specified class.
        filter_by_class(input_csv, filtered_filename, args.class_name)

    fasta_filename = os.path.join(taxonomy_dir, f"{output_base}.fasta")
    main_fasta = csv_to_fasta(input_csv, fasta_filename)

else:
    # Obtaining classification infomation
    main_fasta = process_with_progress(args.class_name)

# When the --tree option is present: Create a phylogenetic tree using VSEARCH→MAFFT→FastTree
if args.tree:
    vsearch_output_name = f"{timestamp}_clustered_sequences.fasta"
    aligned_fasta_name = f"{timestamp}_aligned_sequences.fasta"
    tree_file_name = f"{timestamp}_{args.method}_phylogenetic_tree.nwk"

    # Run VSEARCH
    vsearch_output_path = run_vsearch(main_fasta, vsearch_output_name, alignment_dir)

    # Run MAFFT
    aligned_fasta_path = run_mafft(vsearch_output_path, aligned_fasta_name)

    # Run FastTree
    tree_file_path = run_fasttree(aligned_fasta_path, tree_file_name,
                                  method=args.method, bootstrap=args.bootstrap,
                                  gamma=args.gamma, outgroup=args.outgroup)

    # Export NEXUS files to phylogeny directory
    nwk = Phylo.read(tree_file_path, 'newick')
    nexus_path = os.path.join(phylogeny_dir, f"{timestamp}_{args.method}_phylogenetic_tree.nex")
    Phylo.write(nwk, nexus_path, 'nexus')
    print("Newick file from FastTree has been converted to NEXUS format.")

    print(f"MCMC: {args.mcmc}, Thinning: {args.thinning}, Burn-in: {args.burnin}, Seed: {args.seed}")
    # Run bPTP and mPTP with organized directories
    run_bptp(tree_file_path, args.mcmc, args.thinning, args.burnin, args.seed, bptp_base_dir)
    run_mptp(tree_file_path, mptp_base_dir)

    # Process PTP outputs and update CSV files
    process_ptp_outputs()

print("\n=== Final Output Summary ===")
verify_directory_structure()
print("Script execution completed.")

# Copy with directory structure intact
print("Copying output files with directory structure...")
final_output_path = copy_output_files_with_structure(output_dir, input_dir)
print(f"Process complete. Output files saved to: {final_output_path}")

print("Checking output files location...")
if os.path.exists(output_dir):
    print(f"Output directory exists at: {output_dir}")
    final_output_path = copy_output_files_with_structure(output_dir, input_dir)
    if final_output_path:
        print(f"Process complete. Output files saved to: {final_output_path}")
    else:
        print(f"Using original output directory: {output_dir}")
else:
    print(f"Error: Output directory not found at {output_dir}")
    # Search for an alternative output directory
    alt_outputs = glob.glob(os.path.join(input_dir, "micum_output_*"))
    if alt_outputs:
        print(f"Found alternative output directories: {alt_outputs}")
        final_output_path = alt_outputs[-1]  # Use the latest
        print(f"Using alternative output: {final_output_path}")
    else:
        print("No output directories found")

print("Script execution completed.")