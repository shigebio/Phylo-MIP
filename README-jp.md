# Phylo-MIP(Phylogeny-based Molecular Identification Pipeline)
## What`s app
- 系統解析～種の判定を簡易的かつほぼ自動的に行うことを目標に作成したパイプラインツールです
- 仮想環境の構成にDockerを使用しており、使用者のOSや環境に依存しません
### 使用用途
環境DNAを使用した系統解析のために作られました

### できること
- `localBLAST`or`BLAST+`で出力したファイルをもとに種名を含むOTUを持ったFASTAファイルの作成(NCBIのデータで構成されたDBに限られます)
  - 配列のアライメント([MAFFT](https://mafft.cbrc.jp/alignment/server/index.html))
  - ハプロタイプの検出([VSEARCH](https://github.com/torognes/vsearch))
  - 系統樹の作成([FASTTree](https://morgannprice.github.io/fasttree/))
  - 種の区分決定解析([bPTP](https://species.h-its.org/), [mPTP](https://github.com/Pas-Kapli/mptp))
- Qiime2 の出力ファイルと Phylo-MIP パイプラインの出力ファイルのマージ: <b>merge_data.py</b>
- MAcOSまたはLinuxOSで実行できることを確認済みですが、WindowsOSでもWSLを利用して実行可能です。

## Installation
1. Dockerの導入
     - https://docs.docker.jp/engine/getstarted/step_one.html
1. Phylo-MIP or クローン
    <details><summary>DLリンク</summary>

    ![image](https://github.com/user-attachments/assets/ad10015a-dbe1-4498-a751-ae2e0c42a47a)

    - `Download ZIP`からDL
    </details>

    Git Clone
    ```
    git clone https://github.com/shigebio/Phylo-MIP
    ```
1. DL or Cloneしてきたファイルのディレクトリへ移動
    ```
    cd /{path to Phylo-MIP}/Phylo-MIP
    ```

1. 仮想環境の構築
    ```
    chmod +x setup.sh
    chmod +x entrypoint.sh
    ./setup.sh
    ```

---
## How to Use Phylo-MIP
**実行の前に**
Phylo-MIP.pyはNCBIおよびGBIFのAPIを利用しています。[NCBI](https://blast.ncbi.nlm.nih.gov/doc/blast-help/developerinfo.html#developerinfo)と[GBIF](https://techdocs.gbif.org/en/openapi/v1/species)のAPIガイドラインを参照してください。
可能な限り、初回実効の前に[こちらのコード](https://github.com/shigebio/Phylo-MIP/blob/main/app/Phylo-MIP.py#L22)のメールアドレスを自身のものに書き換えい。

NCBIのガイドラインより一部抜粋
>- Do not contact the server more often than once every 10 seconds.
>- Do not poll for any single RID more often than once a minute.
>- Use the URL parameter email and tool, so that the NCBI can contact you if there is a problem.
>- Run scripts weekends or between 9 pm and 5 am Eastern time on weekdays if more than 50 searches will be submitted.
---

1. コマンドの実行

    **基本のコマンド**
      ```
      phylo-mip {入力CSVファイルのパス} --tree {オプション}
      ```
    **例**
      ```
      phylo-mip ./paht/your_input.csv --tree --method ML --bootstrap 250
      ```

   インプット用CSVの作成方法は[こちら](https://github.com/shigebio/Phylo-MIP/blob/main/README-Preparing_the_input_files.jp.md)を参照

      <details><summary>オプション</summary>

      - `--top` : `qseqid`が同じものを`pident`の上位から1~10まで指定できます
        - デフォルト：`1`
      - `--o` : outputファイルの出力名。無い場合はデフォルト名が適用されます
      - `--class` : 特定の分類群(現段階では綱まで)のみに絞って解析を行うことができます
        - NCBIやGBIFから取得した分類群情報が書き込まれていない場合は、検索対象にならないので注意が必要です
        - 手動で`class`カラムを作成・入力することでも使用可能です
      - `--tree`
          ```
          ## 子オプション
          # --method : 系統樹作成手法の選択ができます。`ML`でML法になります。
          # --bootstrap : ブートストラップの反復回数の指定。デフォルト：`250`。
          # --gamma : ガンマ分布を適用するかどうか。デフォルト：`False`。
          # --outgroup {OTU名} :  外群の指定
      - `--onlyp`： 系統解析以降を実行(アライメント→同一ハプロタイプ除去→系統樹作成→種決定解析)
          - https://github.com/shigebio/MICUM/pull/7
        - `--bptp` : bPTP解析のオプション
          ```
          ## Sub options
          # --mcmc : MCMC chainの回数を設定できます。デフォルト：100000
          # --thinning : サンプリング数を設定できます。サンプリングは、指定された数ごとに各 MCMC チェーンに対して実行されます。デフォルト: `100` = 100回ごとに1回サンプリング
          # --burnin : burn-in率 (0.1~1.0)。デフォルト: `0.1`
          # --seed : 割り当てるシード。同じシード値を指定すると、同じ入力に対して毎回同じ結果が得られます。デフォルトはランダム シードです。
      </details>

   **FASTAファイルとCSVファイル出力だけしたい場合**
      ```
      phylo-mip {入力CSVファイル名} {出力ファイル名}
      ```
    **例**
      ```
      phylo-mip ./path/your_data.csv output_mame
      ```

## 出力
```
=== Directory Structure Verification ===
micum_output_{実行時間}/
  taxonomy/ # 分類情報を付与したCSV、FASTAファイル
    taxonomic_data.csv
    taxonomic_sequences.fasta
  phylogeny/ # fattreeによる系統樹
    {実行時間}_ML_phylogenetic_tree.nex # ファイルの頭には実行時のprefixが付きます
    {実行時間}_ML_{実行時間}_ML_phylogenetic_tree.nwk
  alignment/ # mafft、VSEARCH実施後のファイル
    {実行時間}_haplotype_clusters.tsv # VSEARCHによるハプロタイプ同定結果TSVファイル
    {実行時間}_haplotype_clusters.csv　# VSEARCHによるハプロタイプ同定結果CSVファイル
    {実行時間}_clustered_sequences.fasta # VSEARCHによる同一ハプロタイプの除去後のファイル
    {実行時間}_aligned_sequences.fasta # MAFFTによるアライメント後のファイル
  mptp/ # mptp出力結果
    {実行時間}_mPTP_analysis/
      {実行時間}_mPTP_species_delimitation.txt # 割り当てられたSpeciesと各Speciesの尤度を含むTXTファイル
      {実行時間}_mPTP_species_delimitation.svg # 樹形ファイル(SVG)
  bptp/ # bptp出力結果
    {実行時間}_bPTP_analysis/
      {実行時間}_bPTP_species_delimitation.PTPhSupportPartition.txt.svg # 簡易ヒューリスティック検索(simple heuristic search)により構築された樹形ファイル(SVG)
      {実行時間}_bPTP_species_delimitation.llh.pdf # MCMC chainのtraceロググラフ(PDF)
      {実行時間}_bPTP_species_delimitation.PTPhSupportPartition.txt.sh.tre # 簡易ヒューリスティック検索(simple heuristic search)により構築された樹形ファイル(TREE)
      {実行時間}_bPTP_species_delimitation.PTPMLPartition.txt # ML法により割り当てられたSpeciesと各Speciesの尤度を含むTXTファイル
      {実行時間}_bPTP_species_delimitation.PTPMLPartition.txt.ml.tre # ML法により構築された樹形ファイル(TREE)
      {実行時間}_bPTP_species_delimitation.PTPMLPartition.txt.png # ML法により構築された樹形ファイル(PNG)
      {実行時間}_bPTP_species_delimitation.PTPllh.txt # MCMC chainのtraceログ
      {実行時間}_bPTP_species_delimitation.PTPhSupportPartition.txt.png # 簡易ヒューリスティック検索(simple heuristic search)により構築された樹形ファイル(PNG)
      {実行時間}_bPTP_species_delimitation.PTPPartitions.txt
      {実行時間}_bPTP_species_delimitation.PTPhSupportPartition.txt # 簡易ヒューリスティック検索(simple heuristic search)により割り当てられたSpeciesと各Speciesの尤度を含むTXTファイル
      {実行時間}_bPTP_species_delimitation.PTPMLPartition.txt.svg # ML法により構築された樹形ファイル(SVG)
      {実行時間}_bPTP_species_delimitation.PTPPartitonSummary.txt # 各手法で分割されたSpeiceis
```

OTUにつけられる配列ごとの分類学的ステータスは、pident (localBLAST 検索での検索配列の一致率) の値によって下記ルールで割り当てられます:
```
pident >= 98.00 : 種名
95.00 <= pident < 98.00 : 属名
90.00 <= pident < 95.00 : 科名
85.00 <= pident < 90.00 : 目名
```
https://github.com/shigebio/Phylo-MIP/blob/main/app/Phylo-MIP.py#L291-L298

<b>＊注意＊
     Entrez API(NCBI)から取得したデータはNCCBIのDBの構造上、分類群情報のカラムがずれて取得される可能性がありますので、一度ファイルの分類群情報をチェックしてください
</b>

---
## How to Use merge_data.py
**コマンドの実行**
    ```
    merge_data -q {Qiimeの出力ファイルのパス} -m {Phylo-MIP pipelineの出力ファイルのパス} -f {任意の主力形式: csv/tsv} -o  {出力ファイル名}
    ```
出力ファイル名のデフォルト: 実行時間のprefix
結合後のファイルは実行時のディレクトリ下に出力されます。


# 引用
個人的な使用の場合はお気兼ねなくご使用ください。出版物の発行に際して、このパイプラインを研究活動などにご使用になられた場合は以下を引用いただけますと幸いです:

**Phylo-MIP: Phylogeny-based Molecular Identification Pipeline for DNA metabarcoding, and assessment of insect communities in subalpine river ecosystems
Takumi Yshida, Shonosuke Shigeta, Yuta Hasebe, Masaki Takenaka
bioRxiv 2025.11.10.687572; doi: https://doi.org/10.1101/2025.11.10.687572**


# 使用時に感じた問題点
→[new issue作成](https://github.com/shigebio/Phylo-MIP/issues)して記載いただけると🙏
