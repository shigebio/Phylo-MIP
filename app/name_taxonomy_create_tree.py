import pandas as pd
from Bio import Entrez, Phylo
import requests
import sys
import time
from datetime import datetime
import random
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
import argparse
import re
import os
import six
import csv
from datetime import datetime
sys.path.append('/usr/local/lib/python3.7/site-packages')
from ete3 import NodeStyle

# NCBI API setup:NCBIに怒られないか心配な人は自分のメールアドレスに変更
Entrez.email = "your_email@example.com"
default_seed = random.randint(1, 10**6)

# オプション引数をパース
parser = argparse.ArgumentParser(description="Process Sequence file and generate phylogenetic trees")
parser.add_argument('input_csv', help="Input CSV file (relative path to ../input)")
parser.add_argument('output_base', help="Output base name (relative path to ../output)")
parser.add_argument('--top', type=int, default=1, choices=range(1, 6), help="Number of top results to retain per qseqid (default: 1, range: 1-5)")
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
bptp_group.add_argument('--seed', type=int, default=default_seed, help='Random seed (default: {default_seed})')

args = parser.parse_args()

# CSV読み込み
df = pd.read_csv(os.path.join('..', 'input', args.input_csv))

# 重複するqseqidの場合はpindentが高い方から選択された数だけ残す
df = df.sort_values('pident', ascending=False).groupby('qseqid').head(args.top)

fasta_list = []
csv_data = []
cache = {}

# 英数字、アンダースコア、ハイフン以外の文字、スペースをハイフンに変換
def sanitize_otu_name(name):
    name = re.sub(r'\.', '_', name)
    return re.sub(r'[^a-zA-Z0-9_>\-]', '_', name)

# NCBIから分類群データ取得する用関数
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

    # 最終的に失敗した場合、エラーメッセージを返す代わりにNoneを返す
    print(f"Failed to fetch data for accession ID: {accessionID}")
    return None

# GBIFから分類群データ取得する用関数
def get_gbif_taxonomic_info(species_name):
    url = f"https://api.gbif.org/v1/species/match?name={species_name}"
    for attempt in range(3):
        try:
            response = requests.get(url, timeout=600)
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

# --classの引数で絞り込む
def filter_by_class(input_csv, output_csv, class_name):
    try:
        # CSVを読み込む
        df = pd.read_csv(input_csv)

        # 'class' カラムが存在するか確認
        if 'class' not in df.columns:
            raise KeyError("'class' column not found in the input CSV.")

        # クラスが一致する行をフィルタリング
        filtered_df = df[df['class'].isin(class_name)]

        # 結果が空でないことを確認
        if len(filtered_df) == 0:
            print("No matching rows found for the given class.")

        # フィルタリングされた結果を出力
        os.makedirs(os.path.dirname(output_csv), exist_ok=True)  # 出力ディレクトリを作成
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
    qseq = row['qseq'].replace("-", "N")  # '-' を 'N' に置換

    try:
        # accessionID(sallacc)からNCBIで種名取得
        records = fetch_ncbi_data(accessionID)
        if not records:
            # NCBIから情報を取得できなかった場合
            taxonomic_name = "Uncertain_taxonomy"  # 不確定な分類群
            fasta_entry = f">{sanitize_otu_name(f'{qseqid}_{accessionID}_{taxonomic_name}_{pident:.2f}')}\n{qseq}\n"
            csv_entry = [qseqid, accessionID, 'Unknown', 'Unknown', taxonomic_name, pident, qseq, 'NCBI Failed']
            return fasta_entry, csv_entry

        organism_name = records[0]['GBSeq_organism']

        # Step 2: Use species name from NCBI to query GBIF
        gbif_info = get_gbif_taxonomic_info(organism_name)

        if gbif_info:
            taxonomic_info = gbif_info
            source = 'GBIF'
        else:
            # GBIFから取得失敗→NCBIから取得
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

        # FASTAエントリーとCSVエントリーを作成
        fasta_entry = f">{sanitize_otu_name(f'{qseqid}_{accessionID}_{taxonomic_name}_{pident:.2f}')}\n{qseq}\n"
        csv_entry = [qseqid, accessionID, class_name, order, taxonomic_name, pident, qseq, source]

        print(f"Processing row {index + 1} of {len(df)}")

        return fasta_entry, csv_entry

    except Exception as e:
        print(f"Error processing row {index + 1}: {e}")
        # エラー時には不確定な分類群名を使用
        taxonomic_name = "Uncertain_taxonomy"
        fasta_entry = f">{sanitize_otu_name(f'{qseqid}_{accessionID}_{taxonomic_name}_{pident:.2f}')}\n{qseq}\n"
        csv_entry = [qseqid, accessionID, 'Unknown', 'Unknown', taxonomic_name, pident, qseq, 'Error']
        return fasta_entry, csv_entry

# APIリクエストの進捗表示用
def process_with_progress(filter_class=None):
    processed_rows = 0
    total_rows = len(df)
    fasta_list = []  # フィルター前のFASTAデータ
    csv_data = []    # フィルター前のCSVデータ
    filtered_fasta_list = []  # フィルター後のFASTAデータ
    filtered_csv_data = []    # フィルター後のCSVデータ

    with ThreadPoolExecutor(max_workers=5) as executor:
        futures = {executor.submit(process_row, index, row): index for index, row in df.iterrows()}
        for future in as_completed(futures):
            index = futures[future]
            try:
                fasta_entry, csv_entry = future.result()
                if fasta_entry and csv_entry:
                    fasta_list.append(fasta_entry)
                    csv_data.append(csv_entry)

                    # フィルタリング処理
                    if filter_class:
                        class_value = csv_entry[2]  # 'class' カラムの値
                        if class_value is None:
                            class_value = ""
                        if class_value.strip().lower() in [item.lower() for item in filter_class]:
                            filtered_fasta_list.append(fasta_entry)
                            filtered_csv_data.append(csv_entry)
                    else:
                        # フィルターなしの場合は全データを追加
                        filtered_fasta_list.append(fasta_entry)
                        filtered_csv_data.append(csv_entry)

                    processed_rows += 1
                    print(f"\rProcessed row {processed_rows} of {total_rows} ({processed_rows / total_rows:.2%})", end='')
            except Exception as e:
                print(f"\nError processing row {index + 1}: {e}")

    print()  # End progress tracking

    output_fasta = save_fasta('../output/pre_filtered_output.fasta', fasta_list)
    save_csv('../output/pre_filtered_output.csv', csv_data)

    print("Processing complete.")

    if filter_class:
        output_fasta = save_fasta('../output/filtered_output.fasta', filtered_fasta_list)
        save_csv('../output/filtered_output.csv', filtered_csv_data)

    return output_fasta

def save_fasta(file_path, fasta_entries):
    with open(file_path, 'w') as f:
        f.writelines(fasta_entries)
    return file_path

def save_csv(file_path, csv_rows):
    with open(file_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['qseqid', 'accessionID', 'class', 'order', 'taxonomic_name', 'pident', 'qseq', 'source']) # ヘッダー行
        writer.writerows(csv_rows)

def csv_to_fasta(input_csv, output_fasta): # FASTAへの変換をprocess_row内の処理とまとめたほうがスリムになるかも
    # CSVを読み込む
    df = pd.read_csv(input_csv)

    # 必須列が存在するか確認
    if not {'qseqid', 'qseq'}.issubset(df.columns): # 最低限OTU名が一意になるようにqseqの存在確認
        raise ValueError("Input CSV must contain 'qseqid' and 'qseq' columns.")

    with open(output_fasta, 'w') as fasta_file:
        for _, row in df.iterrows():
            # OTU名を生成
            otu_name_parts = [
                str(row[col]).replace(" ", "_").replace(",", "_").replace(".", "_").replace("-", "_")
                for col in df.columns if col != 'qseq' and pd.notna(row[col])
            ]
            otu_name = "_".join(otu_name_parts)

            # 配列情報を取得
            sequence = row['qseq']

            # FASTAエントリを書き込む
            fasta_file.write(f">{otu_name}\n{sequence}\n")

    return output_fasta

# VSEARCH
def run_vsearch(input_fasta, output_centroids):
    # 出力ディレクトリを指定
    output_dir = "../output/"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # VSEARCHコマンド
    haplotype_tsv = os.path.join(output_dir, "haplotype_list.tsv")
    output_centroids_path = os.path.join(output_dir, output_centroids)  # Centroidsの出力パスを設定
    vsearch_cmd = f"vsearch --cluster_fast {input_fasta} -id 1 --centroids {output_centroids_path} --mothur_shared_out {haplotype_tsv}"
    print(f"Running VSEARCH: {vsearch_cmd}")
    subprocess.run(vsearch_cmd, shell=True, check=True)

    # TSVファイルをCSVに変換
    haplotype_df = pd.read_csv(haplotype_tsv, sep="\t")
    csv_output_file = os.path.join(output_dir, "haplotype_list.csv")
    haplotype_df.to_csv(csv_output_file, index=False)
    print(f"Haplotype data saved to {csv_output_file}")


# MAFFT
def run_mafft(vsearch_output, output_aligned):
    mafft_cmd = f"mafft --auto {vsearch_output} > ../output/{output_aligned}"
    print(f"Running MAFFT: {mafft_cmd}")
    subprocess.run(mafft_cmd, shell=True, check=True)

# FastTree
def run_fasttree(input_aligned, output_tree, method="NJ", bootstrap=250, gamma=False, outgroup=None):
    fasttree_cmd = "fasttree -nt"

    if method == "NJ":
        fasttree_cmd += " -nj"
    elif method == "ML":
        # ML法はデフォルト
        pass
    else:
        raise ValueError(f"Invalid method '{method}'. Choose 'ML' or 'NJ'.")

    if bootstrap > 0:
        fasttree_cmd += f" -boot {bootstrap}"

    if gamma:
        fasttree_cmd += " -gamma"

    if outgroup:
        fasttree_cmd += f" -outgroup {outgroup}"

    fasttree_cmd += f" {input_aligned} > ../output/{output_tree}"
    print(f"Running FastTree: {fasttree_cmd}")
    subprocess.run(fasttree_cmd, shell=True, check=True)
    print(f"Phylogenetic tree saved to {output_tree}")

# bPTP
def run_bptp(tree_file, mcmc, thinning, burnin, seed):
    try:
        # bPTP.pyのパスを構築
        bptp_path = os.path.join('/app/PTP/bin', 'bPTP.py')

        # 出力ディレクトリの作成
        now = datetime.now().strftime('%Y%m%d_%H%M%S')  # 年月日と分秒を取得
        output_dir = os.path.join('../output', f'bPTP_{now}')
        os.makedirs(output_dir, exist_ok=True)  # ディレクトリを作成（存在する場合はスキップ）

        # 出力ファイル名を設定
        output_file = os.path.join(output_dir, 'output_base_tree_bptp_output.txt')

        # bPTP.pyを実行するコマンドを構築
        bptp_command = [
            'xvfb-run', '-a', 'python3', bptp_path,
            '-t', tree_file,
            '-o', output_file,
            '-s', str(seed),
            '-i', str(mcmc),
            '-n', str(thinning),
            '-b', str(burnin)
        ]

        # コマンドを実行
        print(f"Running bPTP with command: {' '.join(bptp_command)}")
        subprocess.run(bptp_command, check=True)
        print('bPTP analysis complete.')
    except subprocess.CalledProcessError as e:
        print(f"Error running bPTP.py: {e}")

# PTP
# 保留
# def run_ptp(tree_file):
#     try:
#         # PTP.pyのパスを構築
#         ptp_path = os.path.join('/app/PTP/bin', 'PTP.py')

#         # 出力ディレクトリの作成
#         now = datetime.now().strftime('%Y%m%d_%H%M%S')  # 年月日と分秒を取得
#         output_dir = os.path.join('../output', f'PTP_{now}')
#         os.makedirs(output_dir, exist_ok=True)  # ディレクトリを作成（存在する場合はスキップ）

#         # 出力ファイル名を設定
#         output_file = os.path.join(output_dir, 'output_base_tree_ptp_output.txt')  # 出力ファイルのパス
#         seed = '1234'  # 任意のシード値を設定

#         # PTP.pyを実行
#         subprocess.run(
#             ['xvfb-run', '-a', 'python3', ptp_path, '-t', tree_file, '-o', output_file, '-s', seed],
#             check=True
#         )
#         print('PTP analysis complete.')
#     except subprocess.CalledProcessError as e:
#         print(f"Error running PTP.py: {e}")

def run_mptp(tree_file):
    try:
        # 出力ディレクトリの作成
        now = datetime.now().strftime('%Y%m%d_%H%M%S')  # 年月日と分秒を取得
        output_dir = os.path.join('../output', f'mPTP_{now}')
        os.makedirs(output_dir, exist_ok=True)  # ディレクトリを作成（存在する場合はスキップ）

        # 出力ファイル名を設定
        output_file = os.path.join(output_dir, 'output_base_tree_mptp_output.txt')  # 出力ファイルのパス

        # mPTPを実行
        subprocess.run(
            ['xvfb-run', '-a', 'mptp', '-tree_file', tree_file,
             '-output_file', output_file, '-ml', '-single'],
            check=True
        )
        print('mPTP analysis complete.')
    except subprocess.CalledProcessError as e:
        print(f"Error running mPTP: {e}")

if args.onlyp: # 系統解析以降だけ実行するための分岐
    input_csv = os.path.join('..', 'input', f"{args.input_csv}")
    filtered_filename = os.path.join('..', 'output', f"{args.output_base}_filtered.csv")  # ../output 下に保存

    if args.class_name:  # --classがあるとき、指定されたclassで絞込
        filter_by_class(input_csv, filtered_filename, args.class_name)

    fasta_filename = os.path.join('..', 'output', f"{args.output_base}.fasta")
    output_fasta = csv_to_fasta(input_csv, fasta_filename)

else:
    # 分類データの取得
    output_fasta = process_with_progress(args.class_name)

# --treeオプションがあるとき：VSEARCH→MAFFT→FastTreeで系統樹作成
if args.tree:
    vsearch_output = os.path.join("../output/", "output_vsearch.fasta")
    aligned_fasta = f"{args.output_base}_aligned.fasta"
    tree_file = f"{args.output_base}_tree.nwk"

    # Run VSEARCH to cluster sequences
    run_vsearch(output_fasta, vsearch_output)

    # Run MAFFT on the VSEARCH output
    run_mafft(vsearch_output, aligned_fasta)

    # Run FastTree
    mafft_output = os.path.join("../output/", aligned_fasta)
    run_fasttree(mafft_output, tree_file, method=args.method, bootstrap=args.bootstrap, gamma=args.gamma, outgroup=args.outgroup)

    # NEXUSファイルの書き出し
    nwk = Phylo.read(f"../output/{tree_file}", 'newick')
    Phylo.write(nwk, f"../output/{args.output_base}_tree.nex", 'nexus')
    print("Newick file from FastTree has been converted to NEXUS format.")

    nexus_output = os.path.join("../output/", tree_file)

    # Run bPTP
    mcmc = args.mcmc
    thinning = args.thinning
    burnin = args.burnin
    seed = args.seed

    print(f"MCMC: {mcmc}, Thinning: {thinning}, Burn-in: {burnin}, Seed: {seed}")
    run_bptp(nexus_output, mcmc, thinning, burnin, seed)

    # Run PTP
    # 保留
    # run_ptp(nexus_output)

    # Run mPTP
    run_mptp(nexus_output)
