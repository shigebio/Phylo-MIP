#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import csv
import os
import sys
import datetime
import re

def parse_arguments():
    parser = argparse.ArgumentParser(description="Insert data from Phylo-MIP output file into Qiime output file (TSV or CSV)")

    parser.add_argument("-q", "--file_qiime", required=True, help="Path to the Qiime output file (TSV or CSV)")
    parser.add_argument("-p", "--file_pm", required=True, help="Path to the Phylo-MIP output file (TSV or CSV)")
    parser.add_argument("-o", "--output", default="", help="Output filename (default: timestamp prefix)")
    parser.add_argument("-f", "--format", choices=["tsv", "csv"], required=True, help="Output format (tsv or csv)")
    parser.add_argument("-d", "--debug", action="store_true", help="Enable debug mode to show detailed matching information")

    return parser.parse_args()

def detect_delimiter(file_path):
    # Detect delimiter in file (TSV or CSV)
    with open(file_path, 'r', encoding='utf-8') as file:
        first_line = file.readline().strip()
        if '\t' in first_line:
            return '\t'
        else:
            return ','

def normalize_id(id_value):
    """
    Normalize ID by removing any whitespace and converting to string.
    This helps with matching IDs between files that might have slight formatting differences.
    """
    if id_value is None:
        return ""
    return str(id_value).strip()

def read_file_into_dict(file_path, key_column, debug=False):
    # Key columns specified by reading the file
    delimiter = detect_delimiter(file_path)
    result_dict = {}
    headers = []
    all_rows = []

    with open(file_path, 'r', encoding='utf-8') as file:
        reader = csv.reader(file, delimiter=delimiter)
        headers = next(reader)

        # Find indexes on key columns
        key_idx = -1
        for i, header in enumerate(headers):
            if header == key_column:
                key_idx = i
                break

        if key_idx == -1:
            print(f"Error: Column '{key_column}' not found in '{file_path}'")
            sys.exit(1)

        # Read all rows first (for debugging)
        all_rows = list(reader)

        if debug:
            print(f"Read {len(all_rows)} rows from {file_path}")
            print(f"Headers: {headers}")
            if len(all_rows) > 0:
                print(f"Sample row: {all_rows[0]}")

        # Process rows and build dictionary
        for row in all_rows:
            if row and len(row) > key_idx:
                # Use normalized ID as key
                normalized_key = normalize_id(row[key_idx])
                if normalized_key:  # Skip empty keys
                    result_dict[normalized_key] = row

        if debug:
            print(f"Created dictionary with {len(result_dict)} entries")

    return headers, result_dict

def merge_files(qiime_file_path, fm_file_path, output_filename, output_format, debug=False):
    # Create a new file by combining the Qiime output file and the Phylo-MIP output file
    # Detect delimiter in qiime output file
    delimiter_q = detect_delimiter(qiime_file_path)

    # Detect delimiter in Phylo-MIP output file
    delimiter_m = detect_delimiter(fm_file_path)

    # Set the delimiter for the output file
    output_delimiter = '\t' if output_format == 'tsv' else ','

    with open(qiime_file_path, 'r', encoding='utf-8') as file_qiime:
        reader_qiime = csv.reader(file_qiime, delimiter=delimiter_q)
        headers_qiime = next(reader_qiime)

    # Find the column index for "#OTU ID" column
    otu_id_idx = -1
    for i, header in enumerate(headers_qiime):
        if header == '#OTU ID':
            otu_id_idx = i
            break

    if otu_id_idx == -1:
        print("Error: '#OTU ID' column not found in Qiime output file")
        sys.exit(1)

    # Read data from Phylo-MIP output file
    headers_pm, data_pm = read_file_into_dict(fm_file_path, 'qseqid', debug)

    # Debug: Analyze ID overlap
    if debug:
        # Collect QIIME IDs for comparison
        qiime_ids = []
        with open(qiime_file_path, 'r', encoding='utf-8') as file_qiime:
            reader_qiime = csv.reader(file_qiime, delimiter=delimiter_q)
            next(reader_qiime)  # Skip header
            for row in reader_qiime:
                if len(row) > otu_id_idx:
                    qiime_ids.append(normalize_id(row[otu_id_idx]))

        _ids = list(data_pm.keys())

        print(f"\nID MATCHING ANALYSIS:")
        print(f"QIIME file has {len(qiime_ids)} OTU IDs")
        print(f"Phylo-MIP file has {len(pm_ids)} qseqids")

        # Check overlap
        common_ids = set(qiime_ids) & set(pm_ids)
        print(f"Number of common IDs: {len(common_ids)}")

        if len(qiime_ids) > 0:
            match_percentage = len(common_ids)/len(qiime_ids)*100
            print(f"Percentage of QIIME IDs matched: {match_percentage:.2f}%")

        # Print some examples
        if len(qiime_ids) > 0:
            print(f"First 3 QIIME IDs: {qiime_ids[:3]}")
        if len(pm_ids) > 0:
            print(f"First 3 Phylo-MIP IDs: {pm_ids[:3]}")

        # Check if IDs need formatting
        if len(common_ids) == 0 and len(qiime_ids) > 0 and len(pm_ids) > 0:
            print("\nWARNING: No matching IDs found! ID format may be different.")
            print(f"QIIME ID example format: '{qiime_ids[0]}'")
            print(f"Phylo-MIP ID example format: '{pm_ids[0]}'")
            print("Consider checking if the ID formats match or need preprocessing.")

    # Get the directory of the qiime file for output
    output_dir = os.path.dirname(qiime_file_path)
    if not output_dir:  # If dirname returns empty string (for relative paths)
        output_dir = os.getcwd()

    # Generate filename with timestamp if not provided
    if not output_filename:
        timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        output_filename = f"{timestamp}_merged"

    # Set the output file extension and path
    output_ext = '.tsv' if output_format == 'tsv' else '.csv'
    if not output_filename.endswith(output_ext):
        output_filename += output_ext

    output_path = os.path.join(output_dir, output_filename)

    # Debug counters
    match_count = 0
    total_rows = 0

    # Create output file
    with open(qiime_file_path, 'r', encoding='utf-8') as file_qiime, open(output_path, 'w', encoding='utf-8', newline='') as output_file:
        reader_qiime = csv.reader(file_qiime, delimiter=delimiter_q)
        writer = csv.writer(output_file, delimiter=output_delimiter)

        headers_qiime = next(reader_qiime)

        # Output headers: Qiime output file headers up to the #OTU ID + all headers from Phylo-MIP output file + remaining headers from Qiime output file
        new_headers = headers_qiime[:otu_id_idx+1] + headers_pm + headers_qiime[otu_id_idx+1:]
        writer.writerow(new_headers)

        # Process each line
        for row_qiime in reader_qiime:
            total_rows += 1

            if len(row_qiime) <= otu_id_idx:
                # Skip if line is too short
                continue

            otu_id = normalize_id(row_qiime[otu_id_idx])

            # Get data matching qseqid from Phylo-MIP output file
            matching_row_pm = data_.get(otu_id, None)

            if matching_row_pm is not None:
                match_count += 1
                if debug and match_count <= 3:
                    print(f"Match found for ID: {otu_id}")
            else:
                matching_row_pm = [""] * len(headers_pm)
                if debug and total_rows <= 10:
                    print(f"No match found for ID: {otu_id}")

            # New line: values ​​up to #OTU ID column in Qiime output file + all values ​​in Phylo-MIP output file + remaining values ​​in Qiime output file
            new_row = row_qiime[:otu_id_idx+1] + matching_row_pm + row_qiime[otu_id_idx+1:]
            writer.writerow(new_row)

    if debug:
        print(f"\nProcessed {total_rows} rows from QIIME file")
        print(f"Found matches for {match_count} rows ({match_count/total_rows*100:.2f}% match rate)")

    print(f"Merge complete: result saved to '{output_path}'")

def main():
    args = parse_arguments()

    # Use the input file paths directly
    qiime_file_path = args.file_qiime
    fm_file_path = args.file_pm
    debug_mode = args.debug

    # Print for debugging
    print(f"Looking for Qiime file at: {qiime_file_path}")
    print(f"Looking for Phylo-MIP file at: {fm_file_path}")
    if debug_mode:
        print("Debug mode enabled - extra information will be displayed")

    # Input file existence check
    if not os.path.exists(qiime_file_path):
        print(f"Error: File '{qiime_file_path}' not found")
        sys.exit(1)

    if not os.path.exists(fm_file_path):
        print(f"Error: File '{fm_file_path}' not found")
        sys.exit(1)

    # Perform a file merge
    merge_files(qiime_file_path, fm_file_path, args.output, args.format, debug_mode)

if __name__ == "__main__":
    main()