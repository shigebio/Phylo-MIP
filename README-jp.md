# MICUM(Moleculer Identification pipeline Computational Unit Manager)
## What`s app
- 系統解析～種の判定を簡易的かつほぼ自動的に行うことを目標に作成したパイプラインツールです
- 仮想環境の構成にDockerを使用しており、使用者のOSや環境に依存しません
### 使用用途
- 環境DNAを使用した系統解析のために作られました

### できること
- `localBLAST`or`BLAST+`で出力したファイルをもとに種名を含むOTUを持ったFASTAファイルの作成(NCBIのデータで構成されたDBに限られます)
  - 配列のアライメント([MAFFT](https://mafft.cbrc.jp/alignment/server/index.html))
  - ハプロタイプの検出([VSEARCH](https://github.com/torognes/vsearch))
  - 系統樹の作成([FASTTree](https://morgannprice.github.io/fasttree/))
  - 種の区分決定解析([bPTP](https://species.h-its.org/), [mPTP](https://github.com/Pas-Kapli/mptp))
- Qiime2 の出力ファイルと MICUM パイプラインの出力ファイルのマージ: <b>merge_data.py</b>

## Installation
1. Dockerの導入 ※このあたりも責務の範疇でないので後で消してqiitaなどに載せようかと考えています
     - https://docs.docker.jp/engine/getstarted/step_one.html
2. GitHubリポジトリをDL or クローン
     - DL
       <details><summary>DLリンク</summary>

       ![image](https://github.com/user-attachments/assets/ad10015a-dbe1-4498-a751-ae2e0c42a47a)

       - `Download ZIP`からDL
      </details>

     - クローン
      ```
      git clone https://github.com/shigebio/MICUM
      ```
3. DL or クローンしたファイルに移動
    ```
    cd /{path to MICUM}/MICUM
    ```
4. 仮想環境の構築
環境によっては`docker`コマンドの前に`sudo`が必要になります
    1. Docker Desktopを起動
         ```
          # dockerの起動確認
          docker version

          # Docker Engineを利用している場合
          servise docker start

          # Docker Composeが無いと表示された場合(Linux)
          sudo apt update
          sudo apt install docker-compose
         ```

   1. 仮想環境の構築・起動
        ```
        docker-compose up -d
        ```
    ---
     - 仮想環境の構築には時間がかかります。
     - イメージをDocker hubから取得することも可能です(上記手順を行った場合は不要です<b>※最新化されていません</b>)。
       - [shigebio/name_taxonomy_create_tree](https://hub.docker.com/r/shigebio/name_taxonomy_create_tree)
     - DL or クローンしてきたファイルの`app`フォルダ下に`input`フォルダ、`output`フォルダが作成されていることを確認してください。

### プログラムの更新方法

  ```
  # update.shのあるパスへ移動
  cd /path/to/update.sh

  # 更新
  bash update.sh
  ```
---
## How to Use MICUM.py
1. 事前に用意したCSVファイルを`input`フォルダ下に移動(インプット用CSVの作成方法は[こちら](#input_section)を参照)
   - 直接指定可能
2. 仮想環境の構築
環境によっては`docker`コマンドの前に`sudo`が必要になります
  1. Docker Desktopを起動
       ```
        # Docker Engineを使用する場合
        service docker start
       ```

 1. 仮想環境の起動
      ```
      docker-compose up -d
      ```
 2. 仮想環境に入る
      ```
      docker exec -it micum /bin/bash
      ```
      - 以下のようになれば、仮想環境に入ることができています。
      ```
      root@absd1234:/app#
      ```
1. コマンドの実行
   1. 基本のコマンド
      ```
      python3 MICUM.py {入力CSVファイル名} --tree {オプション}
      ```
      例
      ```
      python3 MICUM.py your_input.csv --tree --method ML --bootstrap 250
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
        <br>
   2. FASTAファイルとCSVファイル出力だけしたい場合
      ```
      python3 MICUM.py {入力CSVファイル名} {出力ファイル名}
      ```
      例
      ```
      python3 MICUM.py your_data.csv output
      ```

2. コンテナの停止
   ```
   sudo docker-compose down
   ```
     - 常にコンテナを動かしているとメモリ消費しそうなので、停止させておくとよさそうです。

<div id="input_section"></div>
<details><summary>入力ファイルの準備</summary>

- `localBLAST`or`BLAST+`で出力した検索結果のCSVファイルを入力ファイルとして使用します。
   ## **localBLAST実行例**
   ### NCBIのデータからDBの作成
   1. メタゲノム解析ツール(e.g. [qiime2](https://qiime2.org/), [USEARCH](https://www.drive5.com/usearch/)...)を利用して、localBLASTに使用可能な配列データを作成してください。
   2. DBとして使用したい配列ファイルを検索し、send toボタンからFASTA形式で出力
      - 例：動物門のmtDNA 16S rRNA領域を検索(https://www.ncbi.nlm.nih.gov/nuccore/advanced)
          ```
          ((((Animalia) AND 16S) NOT whole genome) NOT chromosome) NOT complete genome
          ```

      <details><summary>send toリンクの場所</summary>

        ![image](https://github.com/user-attachments/assets/7424ed3d-86ba-4afd-96b1-dec875544b98)
      </details>

   3. `localBLAST`or`BLAST+`でDB作成
      ```
      makeblastdb -in {DBとして扱いたいFASTAファイル} -dbtype nucl -out {任意のDB名}.nc -hash_index -parse_seqids
      ```
       - 例：動物門のmtDNA 16S rRNA領域をDBにする
        ```
        makeblastdb -in animalia_16S.fasta -dbtype nucl -out animalia_16S_db.nc -hash_index -parse_seqids
        ```
   4. 作成したDBに対してBLAST検索をかけたい配列のFASTAファイルでBLAST検索
         ```
         blastn -db {検索対象DB名}.nc -query {BLAST検索をしたい配列を含むFASTAファイル名} -out {出力ファイル名}.csv -outfmt "10 qseqid sallacc pident qseq" -max_target_seqs 10 -evalue 1e-40 && sed -i '1i qseqid,sallacc,pident,qseq' {出力ファイル名}.csv
         ```
         - ※`sed`コマンドはLinuxコマンドなので、他OSの場合は注意
         - 例：`2.`で作成したDBに対して、BLAST検索をかけたいFASTAファイル`query_sequence.fasta`でBLAST検索したい場合
         ```
         blastn -db animalia_16S_db.nc -query query_sequence.fasta -out output_quried.csv -outfmt "10 qseqid sallacc pident qseq" -max_target_seqs 10 -evalue 1e-40 && sed -i '1i qseqid,sallacc,pident,qseq' output_quried.csv
         ```
         - `-outfmt "10 xx yy"`と`sed -i '1i xx,yy'`の項目と順番は揃えてください
         - `&& sed`コマンドが通らない場合は、下記入力ファイル例と同じ形式でカラム名を手動でつければ大丈夫です
         - `-outfmt "10 xx yy zz"`の引数は任意で増やすことができますが、最低限、`qseqid` `sallacc` `pident` `qseq`があれば動きます

   ### 入力ファイル例
   CSV形式です
   | qseqid | sallacc | pident | qseq |
   | ---- | ---- | ---- | ---- |
   | 9534cfe94fa593ed71 | AB1234 | 98.805 | GATCGAT・・・ |
   | 9534cfe94fa593ed72 | AB2345 | 96.016 | GATCGAT・・・ |
   | 9534cfe94fa593ed73 | AB3456 | 96.032 | GATCGAT・・・ |
   | 9534cfe94fa593ed74 | AB4567 | 96.032 | GATCGAT・・・ |

  <details><summary>各項目について</summary>

   - `qseqid`
     - BLAST検索実行時につく通し番号で、クエリ上の通し番号でサンプルごとに割り当てられます
     - この値をもとにBLAST結果を`pident`の高い順から選択してデータセットに含めるオプションがあります
   - `sallacc`
     - NCBIのデータをDBに使用する場合、`Accession ID`に相当します
     - 種名の検索に使用します
   - `pident`
     - BLAST検索実行時に出力されるサンプル配列とリファレンス配列の一致率
   - `qseq`
     - BLAST検索に使用したデータの塩基配列
  </details>

- `localBLAST`,`BLAST+`などで検索した後のファイルを使用する想定ですが、上記形式と一致していれば、手動で作成しても問題ないです
</details>

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
- `{input/出力ファイル名}_aligned.fasta`

### VSEARCHによる同一ハプロタイプの除去後のファイル
- `{出力ファイル名}_vsearch.fasta`
- このファイルが系統樹作成や種区分の決定解析に使用されています

### bPTP解析の結果
- `bPTP_{作成日時}`フォルダが作成されます
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
- `mPTP_{作成日時}`フォルダが作成されます
  <details><summary>主要なファイル</summary>

    - `output_base_tree_mptp_{出力ファイル名}.txt.txt`
      - ML法による解析結果。テキスト形式。
    - `output_base_tree_mptp_{出力ファイル名}.txt.svg`
      - ML法によるの解析結果。画像(svg)形式。
  </details>

---
## How to Use merge_data.py
**実行の前に**
MICUM.pyはNCBIおよびGBIFのAPIを利用しています。[NCBI](https://blast.ncbi.nlm.nih.gov/doc/blast-help/developerinfo.html#developerinfo)と[GBIF](https://techdocs.gbif.org/en/openapi/v1/species)のAPIガイドラインを参照してください。
可能な限り、初回実効の前に[こちらのコード](https://github.com/shigebio/MICUM/blob/main/app/MICUM.py#L21)のメールアドレスを自身のものに書き換えてください。

NCBIのガイドラインより一部抜粋
>- Do not contact the server more often than once every 10 seconds.
>- Do not poll for any single RID more often than once a minute.
>- Use the URL parameter email and tool, so that the NCBI can contact you if there is a problem.
>- Run scripts weekends or between 9 pm and 5 am Eastern time on weekdays if more than 50 searches will be submitted.
---
1. MICUM.py を実行したときと同じ方法で Docker コンテナを起動し、仮想環境 (MICUM.py を実行したときと同じ仮想環境) に入る。
   - [MICUM.py の使用方法](https://github.com/shigebio/MICUM?tab=readme-ov-file#how-to-use-micum.py)の手順 2 ～ 4 に従ってください。※MICUM.py を実行するときにこの手順を実行した場合は、この手順を実行する必要はありません。
2. QiimeとMICUMパイプラインの出力ファイルをinputフォルダに移動する
3. コマンドの実行
   - 基本のコマンド
    ```
    python3 merge_data.py -q {Qiime output file name} -m {MICUM pipeline output file name} -f {The file format you want to output: csv/tsv}
    ```
    <details><summary>オプション</summary>

      - `-o` : 出力ファイル名
        - Default: 実行時間のprefix
    </details>

4. 結合後のファイルはoutputファイルに出力されます。

# 使用時に感じた問題点
→[new issue作成](https://github.com/shigebio/MICUM/issues)して記載いただけると🙏
