"""Taxonomy retrieval, rank-aware resolution, and taxonomy file generation."""

import csv
import os
import re
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from typing import Optional

import pandas as pd
import requests
from Bio import Entrez

from .output import ensure_directory_exists


Entrez.email = "your_email@example.com"

TAXONOMY_RANKS = ("kingdom", "phylum", "class", "order", "family", "genus", "species")
LEGACY_COLUMNS = (
    "qseqid", "accessionID", "class", "order", "family", "taxonomic_name", "pident", "qseq", "source"
)
TAXONOMY_COLUMNS = (
    "ncbi_taxid",
    *(f"ncbi_{rank}" for rank in TAXONOMY_RANKS),
    *(f"gbif_{rank}" for rank in TAXONOMY_RANKS),
    *(f"resolved_{rank}" for rank in TAXONOMY_RANKS),
    "resolved_source", "taxonomy_conflict", "taxonomy_conflict_ranks",
    "gbif_status", "gbif_match_type", "gbif_confidence",
)
CSV_COLUMNS = LEGACY_COLUMNS + TAXONOMY_COLUMNS


def sanitize_otu_name(name):
    name = re.sub(r"\.", "_", name)
    return re.sub(r"[^a-zA-Z0-9_>\-]", "_", name)


def _empty_ranks():
    return {rank: None for rank in TAXONOMY_RANKS}


def _normalise_taxon_name(value):
    return str(value).strip().casefold() if value is not None else None


@dataclass
class TaxonomyRecord:
    """Canonical per-row taxonomy record used before legacy serialization."""

    qseqid: str
    accession_id: str
    pident: float
    qseq: str
    legacy_source: str
    ncbi_taxid: Optional[str] = None
    ncbi_ranks: dict = field(default_factory=_empty_ranks)
    gbif_ranks: dict = field(default_factory=_empty_ranks)
    resolved_ranks: dict = field(default_factory=_empty_ranks)
    resolved_source: str = "unresolved"
    taxonomic_name: str = "Uncertain_taxonomy"
    taxonomy_conflict: bool = False
    taxonomy_conflict_ranks: tuple = ()
    gbif_metadata: dict = field(default_factory=dict)

    def to_csv_row(self):
        """Serialize while retaining the existing nine leading columns."""
        row = [
            self.qseqid,
            self.accession_id,
            self.resolved_ranks.get("class"),
            self.resolved_ranks.get("order"),
            self.resolved_ranks.get("family"),
            self.taxonomic_name,
            self.pident,
            self.qseq,
            self.legacy_source,
            self.ncbi_taxid,
        ]
        row.extend(self.ncbi_ranks.get(rank) for rank in TAXONOMY_RANKS)
        row.extend(self.gbif_ranks.get(rank) for rank in TAXONOMY_RANKS)
        row.extend(self.resolved_ranks.get(rank) for rank in TAXONOMY_RANKS)
        row.extend([
            self.resolved_source,
            str(self.taxonomy_conflict).lower(),
            ";".join(self.taxonomy_conflict_ranks),
            self.gbif_metadata.get("gbif_status"),
            self.gbif_metadata.get("gbif_match_type"),
            self.gbif_metadata.get("gbif_confidence"),
        ])
        return row


class TaxonomyProcessor:
    def __init__(self, taxonomy_dir, cache=None, taxonomy_cache=None):
        self.taxonomy_dir = str(taxonomy_dir)
        self.cache = cache if cache is not None else {}
        self.taxonomy_cache = taxonomy_cache if taxonomy_cache is not None else {}

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

    @staticmethod
    def extract_ncbi_taxid(records):
        """Extract a TaxID from a nucleotide record's db_xref qualifiers."""
        def walk(value):
            if isinstance(value, dict):
                if value.get("GBQualifier_name") == "db_xref":
                    qualifier = str(value.get("GBQualifier_value", ""))
                    if qualifier.lower().startswith("taxon:"):
                        return qualifier.split(":", 1)[1]
                for child in value.values():
                    found = walk(child)
                    if found:
                        return found
            elif isinstance(value, (list, tuple)):
                for child in value:
                    found = walk(child)
                    if found:
                        return found
            return None

        return walk(records)

    @staticmethod
    def _parse_ncbi_taxonomy(records):
        """Convert an Entrez Taxonomy response into a rank-keyed mapping."""
        ranks = _empty_ranks()
        if not records:
            return ranks
        taxon = records[0] if isinstance(records, (list, tuple)) else records
        lineage = taxon.get("LineageEx", []) if isinstance(taxon, dict) else []
        entries = list(lineage)
        if isinstance(taxon, dict) and taxon.get("Rank") and taxon.get("ScientificName"):
            entries.append(taxon)
        for entry in entries:
            if not isinstance(entry, dict):
                continue
            rank = str(entry.get("Rank", "")).strip().lower()
            if rank in ranks and entry.get("ScientificName"):
                ranks[rank] = entry["ScientificName"]
        return ranks

    def fetch_ncbi_taxonomy(self, taxid):
        """Fetch rank-aware NCBI Taxonomy data for a TaxID with retry."""
        if not taxid:
            return None
        taxid = str(taxid)
        if taxid in self.taxonomy_cache:
            return self.taxonomy_cache[taxid]
        for attempt in range(3):
            handle = None
            try:
                handle = Entrez.efetch(db="taxonomy", id=taxid, retmode="xml")
                records = Entrez.read(handle)
                ranks = self._parse_ncbi_taxonomy(records)
                if any(value is not None for value in ranks.values()):
                    self.taxonomy_cache[taxid] = ranks
                    return ranks
                raise ValueError("NCBI Taxonomy response did not contain supported ranks")
            except Exception as exc:
                print(f"NCBI Taxonomy request failed: {exc}. Retrying... ({attempt + 1}/3)")
                time.sleep(2 ** attempt)
            finally:
                if handle is not None:
                    handle.close()
        return None

    def get_gbif_taxonomic_info(self, species_name):
        url = f"https://api.gbif.org/v1/species/match?name={species_name}"
        for attempt in range(3):
            try:
                response = requests.get(url, timeout=30)
                if response.status_code == 200:
                    data = response.json()
                    if "order" in data and data["order"] is not None:
                        result = {"source": "GBIF"}
                        for rank in ("species", "genus", "family", "order", "class"):
                            result[rank] = data.get(rank)
                        for rank in ("kingdom", "phylum"):
                            if rank in data:
                                result[rank] = data.get(rank)
                        optional = {
                            "gbif_status": data.get("status"),
                            "gbif_match_type": data.get("matchType"),
                            "gbif_confidence": data.get("confidence"),
                        }
                        result.update({key: value for key, value in optional.items() if value is not None})
                        return result
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

    @staticmethod
    def _legacy_lineage_ranks(organism_name, lineage):
        values = str(lineage or "").split("; ") if lineage else []
        ranks = _empty_ranks()
        ranks["species"] = organism_name
        ranks["genus"] = values[-1] if values else None
        ranks["family"] = values[-2] if len(values) > 1 else None
        ranks["order"] = values[-3] if len(values) > 2 else None
        ranks["class"] = values[-4] if len(values) > 3 else None
        return ranks

    @staticmethod
    def _ranked_info(info):
        ranks = _empty_ranks()
        if info:
            for rank in TAXONOMY_RANKS:
                if rank in info:
                    ranks[rank] = info.get(rank)
        return ranks

    @staticmethod
    def _conflicts(ncbi_ranks, gbif_ranks):
        return tuple(
            rank for rank in TAXONOMY_RANKS
            if ncbi_ranks.get(rank) and gbif_ranks.get(rank)
            and _normalise_taxon_name(ncbi_ranks[rank]) != _normalise_taxon_name(gbif_ranks[rank])
        )

    @staticmethod
    def _taxonomic_name(ranks, pident):
        selected_rank = None
        if pident >= 98.00:
            selected_rank = "species"
        elif 95.00 <= pident < 98.00:
            selected_rank = "genus"
        elif 92.00 <= pident < 95.00:
            selected_rank = "family"
        elif 85.00 <= pident < 92.00:
            selected_rank = "order"
        if selected_rank is None:
            return "Low_Identity_Match"
        return ranks.get(selected_rank) or "Uncertain_taxonomy"

    def process_row(self, index, row):
        qseqid, accessionID, pident = row["qseqid"], row["sallacc"], row["pident"]
        qseq = row["qseq"].replace("-", "N")
        try:
            records = self.fetch_ncbi_data(accessionID)
            if not records:
                record = TaxonomyRecord(qseqid, accessionID, pident, qseq, "NCBI Failed")
                fasta_entry = f">{sanitize_otu_name(f'{qseqid}_{accessionID}_{record.taxonomic_name}_{pident:.2f}')}\n{qseq}\n"
                return fasta_entry, record.to_csv_row()

            first_record = records[0]
            organism_name = first_record.get("GBSeq_organism")
            ncbi_taxid = self.extract_ncbi_taxid(records)
            ncbi_ranks = self.fetch_ncbi_taxonomy(ncbi_taxid) if ncbi_taxid else None
            ncbi_taxonomy_found = bool(ncbi_ranks and any(ncbi_ranks.values()))
            fallback_ranks = self._legacy_lineage_ranks(organism_name, first_record.get("GBSeq_taxonomy"))
            if not ncbi_taxonomy_found:
                ncbi_ranks = fallback_ranks

            gbif_info = self.get_gbif_taxonomic_info(organism_name) if organism_name else None
            gbif_ranks = self._ranked_info(gbif_info)
            conflict_ranks = self._conflicts(ncbi_ranks, gbif_ranks) if ncbi_taxonomy_found else ()
            if ncbi_taxonomy_found:
                resolved_ranks, resolved_source = ncbi_ranks, "NCBI"
            elif any(gbif_ranks.values()):
                resolved_ranks, resolved_source = gbif_ranks, "GBIF"
            elif any(fallback_ranks.values()):
                resolved_ranks, resolved_source = fallback_ranks, "NCBI_lineage_fallback"
            else:
                resolved_ranks, resolved_source = _empty_ranks(), "unresolved"

            legacy_source = "GBIF" if gbif_info else "NCBI"
            taxonomic_name = self._taxonomic_name(resolved_ranks, pident)
            record = TaxonomyRecord(
                qseqid, accessionID, pident, qseq, legacy_source,
                ncbi_taxid=ncbi_taxid,
                ncbi_ranks=ncbi_ranks,
                gbif_ranks=gbif_ranks,
                resolved_ranks=resolved_ranks,
                resolved_source=resolved_source,
                taxonomic_name=taxonomic_name,
                taxonomy_conflict=bool(conflict_ranks),
                taxonomy_conflict_ranks=conflict_ranks,
                gbif_metadata={key: gbif_info.get(key) for key in ("gbif_status", "gbif_match_type", "gbif_confidence")}
                if gbif_info else {},
            )
            fasta_entry = f">{sanitize_otu_name(f'{qseqid}_{accessionID}_{taxonomic_name}_{pident:.2f}')}\n{qseq}\n"
            return fasta_entry, record.to_csv_row()
        except Exception as exc:
            print(f"Error in qseqid {qseqid}: {exc}")
            record = TaxonomyRecord(qseqid, accessionID, pident, qseq, "Error")
            fasta_entry = f">{sanitize_otu_name(f'{qseqid}_{accessionID}_{record.taxonomic_name}_{pident:.2f}')}\n{qseq}\n"
            return fasta_entry, record.to_csv_row()

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
                        fasta_list.append(fasta_entry)
                        csv_data.append(csv_entry)
                        if filter_class:
                            class_value = csv_entry[2] or ""
                            if class_value.strip().lower() in [item.lower() for item in filter_class]:
                                filtered_fasta_list.append(fasta_entry)
                                filtered_csv_data.append(csv_entry)
                        else:
                            filtered_fasta_list.append(fasta_entry)
                            filtered_csv_data.append(csv_entry)
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
                writer.writerow(CSV_COLUMNS)
                for row in csv_rows:
                    values = list(row)
                    if len(values) < len(CSV_COLUMNS):
                        values.extend([None] * (len(CSV_COLUMNS) - len(values)))
                    writer.writerow(["" if value is None else value for value in values[:len(CSV_COLUMNS)]])
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
