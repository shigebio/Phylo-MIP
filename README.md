# MICUM(Moleculer Identification pipeline Computational Unit Manager：仮称)
## これはなに
- 系統解析～種の判定を簡易的かつほぼ自動的に行うことを目標に作成したパイプラインツールです
- Dockerを使用しており、使用者の環境に依存しません
- 使用用途
  - 環境DNA

## 入力ファイルの準備
### 入力ファイル例
CSV形式です
| qseqid | sallacc | pident | qseq |
| ---- | ---- | ---- | ---- |
| 9534cfe94fa593ed71 | AB1234 | 98.805 | GATCGAT・・・ |
| 9534cfe94fa593ed72 | AB2345 | 96.016 | GATCGAT・・・ |
| 9534cfe94fa593ed73 | AB3456 | 96.032 | GATCGAT・・・ |
| 9534cfe94fa593ed74 | AB4567 | 96.032 | GATCGAT・・・ |
#### 各項目について
- `qseqid`
  - `localBLAST`実行時につく通し番号で、クエリ上の通し番号でサンプルごとに割り当てられます
  - ここの値をもとにBLAST結果を`pident`の高い順から選択してデータセットに含めるオプションがあります
- `sallacc`
  - NCBIの`Accession ID`
  - 種名の検索に使います
- `pident`
  - `localBLAST`実行時に出力されるサンプル配列とリファレンス配列の一致率
- `qseq`
  - サンプル配列の塩基配列情報

- `localBLAST`,`BLAST+`などで行う想定ですが、最終的にinputが要求する形式と一致していれば、手動で作成しても問題ないです

### inputファイル作成
- このプロジェクトの説明ではないので、全体のフローはあとでPDFなどで別添えにする予定
#### **localBLAST実行例**
- 配列情報の取得
    1. GenBankで下記のように検索し、send toからFASTA形式で出力
       1. `Animalia 16S NOT "whole genome" NOT "chromosome" NOT "complete genome" 423750'`
    2. localBLAST、またはBLAST+でDB作成
        ```
        makeblastdb -in {DLしてきたFASTAファイル名} -dbtype nucl -out {出力したいファイル名}.nc -hash_index -parse_seqids
    3. 作成したDBに対して検索をかけたい配列のFASTAファイルでBLAST検索
        ```
        blastn -db {2.で作成したDB名}.nc -query {検索かけたいFASTAファイル名} -out {出力したいファイル名}.csv -outfmt "10 qseqid sseqid sallacc length pident mismatch gapopen qstart qend sstart send evalue bitscore qseq" -max_target_seqs 10 -evalue 1e-40 && sed -i '1i qseqid,sseqid,sallacc,length,pident,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore,qseq' {出力したいファイル名}.csv
    - `-outfmt "10 xx yy"`と`sed -i '1i xx, yy'`の項目と順番は揃えてください
    - 最低限、先述の4項目があれば動きます
## プログラムを使う
1. Dockerの導入
     - https://docs.docker.jp/engine/getstarted/step_one.html
2. リポジトリをDL or クローン
     - DL or クローンする場所は、デスクトップなどファイルシステムからアクセスしやすい場所がおすすめです
     - クローンする場合
       `git clone https://github.com/shigebio/MICUM`
3. クローンしたファイルに移動
  `cd /{path to MICUM}/MICUM`
       - クローンした場所にファイルが作られるので、`cd MICUM`でも行けると思います
1. Dockerの起動
   ```
   sudo service docker start
1. Docker環境のビルド
    ```
    $ sudo docker-compose build
    # 後述するコンテナの立ち上げと同時に行う場合
    $ sudo docker-compose up --build -d
    # おそらくこれだとbuild一瞬？
    $ sudo docker pull shigebio/name_taxonomy_create_tree
    ```
  - イメージのDocker hubレポジトリ：[shigebio/name_taxonomy_create_tree](https://hub.docker.com/r/shigebio/name_taxonomy_create_tree)
2. コンテナの起動
  `sudo docker-compose up -d`
     - ビルドした環境を起動しています
     - `sudo docker-compose up --build -d`やった場合は、やらなくていいです
1. コンテナ内に移動
  `sudo docker-compose exec -it app /bin/bash`
     - `root@ca174dd9ea32:/app# `コンソールがこんなかんじになれば大丈夫です
1. 事前に用意したCSVファイルを`imput`下に移動
2. コマンド実行
     - FASTAファイルとCSVファイル出力だけしたい場合
       - `python name_taxonomy_create_tree.py {入力するCSVファイル名} {出力したいファイル名}`
       - 例：`python name_taxonomy_create_tree.py your_data.csv output`
       - オプション
         - `--top`：`qseqid`が同じものを`pident`の上位から1~5まで指定できます
     - 系統樹作成～種同定(PTP)まで
       - `python3 name_taxonomy_create_tree.py {入力するCSVファイル名} {出力したいファイル名} --tree {各種オプション}`
       - 例：`python3 name_taxonomy_create_tree.py your_input.csv output --tree --method NJ --bootstrap 250`
       - オプション
           - `--method`：系統樹作成手法の選択。デフォルトでは`NJ`になっています。`ML`と`NJ`が選べます
           - `--bootstrap`：ブートストラップ値の指定。デフォルトでは`250`になっています。
           - `--gamma`：ガンマ分布を適用するかどうか。デフォルトでは`False`になっています。
           - `--outgroup`：外群の指定。デフォルトでは`None`になっています。
1. コンテナの停止
   - `sudo docker-compose down`
     - ずっとコンテナ動かしているとメモリ消費しそうなので、停止させておくとよさそうです
     - upで再起動→手順2から再開できます
# 出力
- 後で書きます
