#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import csv
import os
import sys
import datetime

def parse_arguments():
    parser = argparse.ArgumentParser(description="Insert data from MICUM output file into Qiime output file (TSV or CSV)")

    parser.add_argument("-q", "--file_qiime", required=True, help="Path to the Qiime output file (TSV or CSV)")
    parser.add_argument("-m", "--file_micum", required=True, help="Path to the MICUM output file (TSV or CSV)")
    parser.add_argument("-o", "--output", default="", help="Output filename (default: timestamp prefix)")
    parser.add_argument("-f", "--format", choices=["tsv", "csv"], required=True, help="Output format (tsv or csv)")

    return parser.parse_args()

def detect_delimiter(file_path):
    # Detect delimiter in file (TSV or CSV)
    with open(file_path, 'r', encoding='utf-8') as file:
        first_line = file.readline().strip()
        if '\t' in first_line:
            return '\t'
        else:
            return ','

def read_file_into_dict(file_path, key_column):
    # Key columns specified by reading the file
    delimiter = detect_delimiter(file_path)
    result_dict = {}
    headers = []

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

        # Read the remaining lines
        for row in reader:
            if row and len(row) > key_idx:
                result_dict[row[key_idx]] = row

    return headers, result_dict

def merge_files(qiime_file_path, micum_file_path, output_filename, output_format):
    # Create a new file by combining the Qiime output file and the MICUM output file
    # Detect delimiter in qiime output file
    delimiter_q = detect_delimiter(qiime_file_path)

    # Detect delimiter in MICUM output file
    delimiter_m = detect_delimiter(micum_file_path)

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

    # Read data from MICUM output file
    headers_micum, data_micum = read_file_into_dict(micum_file_path, 'qseqid')

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

    # Create output file
    with open(qiime_file_path, 'r', encoding='utf-8') as file_qiime, open(output_path, 'w', encoding='utf-8', newline='') as output_file:
        reader_qiime = csv.reader(file_qiime, delimiter=delimiter_q)
        writer = csv.writer(output_file, delimiter=output_delimiter)

        headers_qiime = next(reader_qiime)

        # Output headers: Qiime output file headers up to the #OTU ID + all headers from MICUM output file + remaining headers from Qiime output file
        new_headers = headers_qiime[:otu_id_idx+1] + headers_micum + headers_qiime[otu_id_idx+1:]
        writer.writerow(new_headers)

        # Process each line
        for row_qiime in reader_qiime:
            if len(row_qiime) <= otu_id_idx:
                # Skip if line is too short
                continue

            otu_id = row_qiime[otu_id_idx]

            # Get data matching qseqid from MICUM output file
            matching_row_micum = data_micum.get(otu_id, [""] * len(headers_micum))

            # New line: values ​​up to #OTU ID column in Qiime output file + all values ​​in MICUM output file + remaining values ​​in Qiime output file
            new_row = row_qiime[:otu_id_idx+1] + matching_row_micum + row_qiime[otu_id_idx+1:]
            writer.writerow(new_row)

    print(f"Merge complete: result saved to '{output_path}'")

def main():
    args = parse_arguments()

    # Use the input file paths directly
    qiime_file_path = args.file_qiime
    micum_file_path = args.file_micum

    # Print for debugging
    print(f"Looking for Qiime file at: {qiime_file_path}")
    print(f"Looking for MICUM file at: {micum_file_path}")

    # Input file existence check
    if not os.path.exists(qiime_file_path):
        print(f"Error: File '{qiime_file_path}' not found")
        sys.exit(1)

    if not os.path.exists(micum_file_path):
        print(f"Error: File '{micum_file_path}' not found")
        sys.exit(1)

    # Perform a file merge
    merge_files(qiime_file_path, micum_file_path, args.output, args.format)

if __name__ == "__main__":
    main()