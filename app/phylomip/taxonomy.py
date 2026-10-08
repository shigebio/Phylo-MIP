"""Taxonomy retrieval, rank selection, and taxonomy file generation."""

import csv
import os
import re
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

import pandas as pd
import requests
from Bio import Entrez

from .output import ensure_directory_exists

Entrez.email = "your_email@example.com"


def sanitize_otu_name(name):
    name = re.sub(r"\.", "_", name)
    return re.sub(r"[^a-zA-Z0-9_>\-]", "_", name)


class TaxonomyProcessor:
    def __init__(self, taxonomy_dir, cache=None):
        self.taxonomy_dir = str(taxonomy_dir)
        self.cache = cache if cache is not None else {}

    def fetch_ncbi_data(self, accessionID):
        if accessionID in self.cache:
            return self.cache[accessionID]
        for attempt in range(3):
            try:
                handle = Entrez.efetch(db="nucleotide", id=accessionID, rettype="gb", retmode="xml")
                records = Entrez.read(handle)
                handle.close()
                self.cache[accessionID] = records
                return records
            except Exception as exc:
                print(f"NCBI API request failed: {exc}. Retrying... ({attempt + 1}/3)")
                time.sleep(2 ** attempt)
        print(f"Failed to fetch data for accession ID: {accessionID}")
        return None

    def get_gbif_taxonomic_info(self, species_name):
        url = f"https://api.gbif.org/v1/species/match?name={species_name}"
        for attempt in range(3):
            try:
                response = requests.get(url, timeout=30)
                if response.status_code == 200:
                    data = response.json()
                    if "order" in data and data["order"] is not None:
                        return {"source": "GBIF", "species": data.get("species"),
                                "genus": data.get("genus"), "family": data.get("family"),
                                "order": data.get("order"), "class": data.get("class")}
                print(f"GBIF API request incomplete or failed (attempt {attempt + 1}/3)")
            except requests.exceptions.RequestException as exc:
                print(f"GBIF API request failed: {exc}. Retrying...")
            time.sleep(1)
        return None

    def filter_by_class(self, input_csv, output_csv, class_name):
        try:
            df = pd.read_csv(input_csv)
            if "class" not in df.columns:
                raise KeyError("'class' column not found in the input CSV.")
            filtered_df = df[df["class"].isin(class_name)]
            if len(filtered_df) == 0:
                print("No matching rows found for the given class.")
            os.makedirs(os.path.dirname(output_csv), exist_ok=True)
            filtered_df.to_csv(output_csv, index=False)
            print(f"Filtered data saved to {output_csv}.")
        except FileNotFoundError:
            print(f"Error: Input file {input_csv} not found.")
        except KeyError:
            print("Error: 'Class' column not found in the input CSV.")
        except Exception as exc:
            print(f"An unexpected error occurred: {exc}")

    def process_row(self, index, row):
        qseqid, accessionID, pident = row["qseqid"], row["sallacc"], row["pident"]
        qseq = row["qseq"].replace("-", "N")
        try:
            records = self.fetch_ncbi_data(accessionID)
            if not records:
                taxonomic_name = "Uncertain_taxonomy"
                fasta_entry = f">{sanitize_otu_name(f'{qseqid}_{accessionID}_{taxonomic_name}_{pident:.2f}')}\n{qseq}\n"
                return fasta_entry, [qseqid, accessionID, "Unknown", "Unknown", "Unknown", taxonomic_name, pident, qseq, "NCBI Failed"]
            organism_name = records[0]["GBSeq_organism"]
            gbif_info = self.get_gbif_taxonomic_info(organism_name)
            if gbif_info:
                taxonomic_info, source = gbif_info, "GBIF"
            else:
                taxonomy_list = records[0]["GBSeq_taxonomy"].split("; ")
                taxonomic_info = {"source": "NCBI", "species": organism_name,
                                  "genus": taxonomy_list[-1] if taxonomy_list else None,
                                  "family": taxonomy_list[-2] if len(taxonomy_list) > 1 else None,
                                  "order": taxonomy_list[-3] if len(taxonomy_list) > 2 else None,
                                  "class": taxonomy_list[-4] if len(taxonomy_list) > 3 else None}
                source = "NCBI"
            taxonomic_name = "Low_Identity_Match"
            if pident >= 98.00:
                taxonomic_name = taxonomic_info.get("species", "Uncertain_taxonomy")
            elif 95.00 <= pident < 98.00:
                taxonomic_name = taxonomic_info.get("genus", "Uncertain_taxonomy")
            elif 92.00 <= pident < 95.00:
                taxonomic_name = taxonomic_info.get("family", "Uncertain_taxonomy")
            elif 85.00 <= pident < 92.00:
                taxonomic_name = taxonomic_info.get("order", "Uncertain_taxonomy")
            class_name = taxonomic_info.get("class", "Unknown")
            order = taxonomic_info.get("order", "Unknown")
            family_name = taxonomic_info.get("family", "Unknown")
            fasta_entry = f">{sanitize_otu_name(f'{qseqid}_{accessionID}_{taxonomic_name}_{pident:.2f}')}\n{qseq}\n"
            return fasta_entry, [qseqid, accessionID, class_name, order, family_name, taxonomic_name, pident, qseq, source]
        except Exception as exc:
            print(f"Error in qseqid {qseqid}: {exc}")
            taxonomic_name = "Uncertain_taxonomy"
            fasta_entry = f">{sanitize_otu_name(f'{qseqid}_{accessionID}_{taxonomic_name}_{pident:.2f}')}\n{qseq}\n"
            return fasta_entry, [qseqid, accessionID, "Unknown", "Unknown", "Unknown", taxonomic_name, pident, qseq, "Error"]

    def process_with_progress(self, df, filter_class=None):
        processed_rows, total_rows = 0, len(df)
        fasta_list, csv_data, filtered_fasta_list, filtered_csv_data = [], [], [], []
        with ThreadPoolExecutor(max_workers=5) as executor:
            futures = {executor.submit(self.process_row, index, row): index for index, row in df.iterrows()}
            for future in as_completed(futures):
                index = futures[future]
                try:
                    fasta_entry, csv_entry = future.result()
                    if fasta_entry and csv_entry:
                        fasta_list.append(fasta_entry); csv_data.append(csv_entry)
                        if filter_class:
                            class_value = csv_entry[2] or ""
                            if class_value.strip().lower() in [item.lower() for item in filter_class]:
                                filtered_fasta_list.append(fasta_entry); filtered_csv_data.append(csv_entry)
                        else:
                            filtered_fasta_list.append(fasta_entry); filtered_csv_data.append(csv_entry)
                        processed_rows += 1
                        print(f"\rProcessed row {processed_rows} of {total_rows} ({processed_rows / total_rows:.2%})", end="")
                except Exception as exc:
                    print(f"\nError processing row {index + 1}: {exc}")
        print()
        taxonomy_fasta = self.save_fasta(os.path.join(self.taxonomy_dir, "taxonomic_sequences.fasta"), fasta_list)
        self.save_csv(os.path.join(self.taxonomy_dir, "taxonomic_data.csv"), csv_data)
        print("Processing complete.")
        if filter_class:
            filtered_fasta = self.save_fasta(os.path.join(self.taxonomy_dir, "class_filtered_taxonomic_sequences.fasta"), filtered_fasta_list)
            self.save_csv(os.path.join(self.taxonomy_dir, "class_filtered_taxonomic_data.csv"), filtered_csv_data)
            return filtered_fasta
        return taxonomy_fasta

    def save_fasta(self, file_path, fasta_entries):
        ensure_directory_exists(file_path)
        try:
            with open(file_path, "w") as handle:
                handle.writelines(fasta_entries)
            print(f"SUCCESS: Saved FASTA to {file_path}")
            return file_path
        except Exception as exc:
            print(f"ERROR: Failed to save FASTA to {file_path}: {exc}")
            raise

    def save_csv(self, file_path, csv_rows):
        ensure_directory_exists(file_path)
        try:
            with open(file_path, "w", newline="") as handle:
                writer = csv.writer(handle)
                writer.writerow(["qseqid", "accessionID", "class", "order", "family", "taxonomic_name", "pident", "qseq", "source"])
                writer.writerows(csv_rows)
            print(f"SUCCESS: Saved CSV to {file_path}")
        except Exception as exc:
            print(f"ERROR: Failed to save CSV to {file_path}: {exc}")
            raise

    def csv_to_fasta(self, input_csv, output_fasta):
        df = pd.read_csv(input_csv)
        if not {"qseqid", "qseq"}.issubset(df.columns):
            raise ValueError("Input CSV must contain 'qseqid' and 'qseq' columns.")
        with open(output_fasta, "w") as fasta_file:
            for _, row in df.iterrows():
                otu_name = "_".join(str(row[col]).replace(" ", "_").replace(",", "_").replace(".", "_").replace("-", "_")
                                     for col in df.columns if col != "qseq" and pd.notna(row[col]))
                fasta_file.write(f">{otu_name}\n{row['qseq']}\n")
        return output_fasta


def csv_to_fasta(input_csv, output_fasta):
    """Compatibility wrapper for callers that only need CSV conversion."""
    return TaxonomyProcessor(os.path.dirname(output_fasta)).csv_to_fasta(input_csv, output_fasta)

