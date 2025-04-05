# MICUM(Moleculer Identification pipeline Computational Unit Manager)
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
- Qiime2 の出力ファイルと MICUM パイプラインの出力ファイルのマージ: <b>merge_data.py</b>

## Installation
1. Dockerの導入
     - https://docs.docker.jp/engine/getstarted/step_one.html
1. MICUMをDL or クローン
    <details><summary>DLリンク</summary>

    ![image](https://github.com/user-attachments/assets/ad10015a-dbe1-4498-a751-ae2e0c42a47a)

    - `Download ZIP`からDL
    </details>

    Git Clone
    ```
    git clone https://github.com/shigebio/MICUM
    ```
1. DL or Cloneしてきたファイルのディレクトリへ移動
    ```
    cd /{path to MICUM}/MICUM
    ```

1. 仮想環境の構築
    ```
    chmod +x setup.sh
    chmod +x entrypoint.sh
    ./setup.sh
    ```

---
## How to Use MICUM
**実行の前に**
MICUM.pyはNCBIおよびGBIFのAPIを利用しています。[NCBI](https://blast.ncbi.nlm.nih.gov/doc/blast-help/developerinfo.html#developerinfo)と[GBIF](https://techdocs.gbif.org/en/openapi/v1/species)のAPIガイドラインを参照してください。
可能な限り、初回実効の前に[こちらのコード](https://github.com/shigebio/MICUM/blob/main/app/MICUM.py#L21)のメールアドレスを自身のものに書き換えい。

NCBIのガイドラインより一部抜粋
>- Do not contact the server more often than once every 10 seconds.
>- Do not poll for any single RID more often than once a minute.
>- Use the URL parameter email and tool, so that the NCBI can contact you if there is a problem.
>- Run scripts weekends or between 9 pm and 5 am Eastern time on weekdays if more than 50 searches will be submitted.
---
1. 事前に用意したCSVファイルを`input`フォルダ下に移動(インプット用CSVの作成方法は[こちら](https://github.com/shigebio/MICUM/blob/main/README-Preparing_the_input_files.jp.md)を参照)
   - 直接指定可能
1. 仮想環境の構築
  Docker Desktopを起動
    ```
    # Docker Engineを使用する場合
    service docker start
    ```

1. コマンドの実行
    **基本のコマンド**
      ```
      micum {入力CSVファイルのパス} --tree {オプション}
      ```
    **例**
      ```
      micum ./paht/your_input.csv --tree --method ML --bootstrap 250
      ```

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
      python3 MICUM.py {入力CSVファイル名} {出力ファイル名}
      ```
    **例**
      ```
      micum ./path/your_data.csv output_mame
      ```

## 出力
### 種名の割りてられたOTUで構成されたFASTAファイル/CSVファイル
- `--class`オプションなし
  - `pre_filtered_{出力ファイル名}.csv`
  - `pre_filtered_{出力ファイル名}.festa`
- `--class`オプションあり
  - `filtered_{出力ファイル名}.csv`
  - `filtered_{出力ファイル名}.festa`

OTUにつけられる配列ごとの分類学的ステータスは、pident (localBLAST 検索での検索配列の一致率) の値によって下記ルールで割り当てられます:
```
pident >= 98.00 : 種名
95.00 <= pident < 98.00 : 属名
90.00 <= pident < 95.00 : 科名
85.00 <= pident < 90.00 : 目名
```
https://github.com/shigebio/MICUM/blob/main/app/MICUM.py#L184-L192

<b>＊注意＊
     Entrez API(NCBI)から取得したデータはNCCBIのDBの構造上、分類群情報のカラムがずれて取得される可能性がありますので、一度ファイルの分類群情報をチェックしてください
</b>

### MAFFTによるアライメント後のファイル
  `{input/出力ファイル名}_aligned.fasta`

### VSEARCHによる同一ハプロタイプの除去後のファイル
  `{出力ファイル名}_vsearch.fasta`
  このファイルが系統樹作成や種区分の決定解析に使用されています

### bPTP解析の結果
  <details><summary>主要なファイル</summary>

    - `output_base_tree_bptp_{出力ファイル名}.txt.PTPhSupportPartition.txt`
      - 簡易ヒューリスティック検索(simple heuristic search)による解析結果。テキスト形式。
    - `output_base_tree_bptp_{出力ファイル名}.txt.PTPhSupportPartition.txt.png`
      - 簡易ヒューリスティック検索(simple heuristic search)による解析結果。画像(png)形式。
    - `output_base_tree_bptp_{出力ファイル名}.txt.PTPhSupportPartition.txt.svg`
      - 簡易ヒューリスティック検索(simple heuristic search)による解析結果。画像(svg)形式。
    - `output_base_tree_bptp_output.txt.PTPMLPartition.txt`
      - ML法による解析結果。テキスト形式。
    - `output_base_tree_bptp_output.txt.PTPMLPartition.txt.png`
      - ML法による解析結果。画像(png)形式。
    - `output_base_tree_bptp_output.txt.PTPMLPartition.txt.svg`
      - ML法による解析結果。画像(svg)形式。

  </details>

### mPTP解析の結果
  <details><summary>主要なファイル</summary>

    - `output_base_tree_mptp_{出力ファイル名}.txt.txt`
      - ML法による解析結果。テキスト形式。
    - `output_base_tree_mptp_{出力ファイル名}.txt.svg`
      - ML法によるの解析結果。画像(svg)形式。
  </details>

---
## How to Use merge_data.py
**コマンドの実行**
    ```
    merge_data -q {Qiimeの出力ファイルのパス} -m {MICUM pipelineの出力ファイルのパス} -f {任意の主力形式: csv/tsv} -o  {出力ファイル名}
    ```
出力ファイル名のデフォルト: 実行時間のprefix
結合後のファイルは実行時のディレクトリ下に出力されます。

# 使用時に感じた問題点
→[new issue作成](https://github.com/shigebio/MICUM/issues)して記載いただけると🙏
