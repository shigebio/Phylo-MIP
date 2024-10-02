# パイプライン名
---
## これはなに
- localBLASTのblastnの出力結果から得られた値をもとにOTUを書き換えるプログラムです
- name_taxonomy.py以外のファイルは入出力例なので消しても動作に問題はありません
---
## 入力ファイルの準備
- ここはlocalBLAST,BLAST+などで行います
- **CSVで出力(`-out_fmt 10`)してください**
- **`qseqid`、`sallacc`、`pident`、`qseq`が必須です
- 例：`st1_kamikochi_max10_BLAST10.csv`
  - このファイルの出どころは、以下手順
    1. GenBankで下記のように検索し、send toからFASTA形式で出力
       1. `Animalia 16S NOT "whole genome" NOT "chromosome" NOT "complete genome" 423750'`
    2. localBLAST、またはBLAST+でDB作成
       1. ```
          $makeblastdb -in {DLしてきたFASTAファイル名} -dbtype nucl -out {出力したいファイル名}.nc -hash_index -parse_seqids
       2. 例：
      	```
     		$makeblastdb -in 16S_genebank_all_240521_plus_Kamikochi.fasta -dbtype nucl -out all_16S.nc -hash_index -parse_seqids
    3. 作成したDBに対して検索をかけたい配列のFASTAファイルでBLAST検索
       1. ```
          $blastn -db {2.で作成したDB名}.nc -query {検索かけたいFASTAファイル名} -out {出力したいファイル名。拡張子は.csv} -outfmt "10 qseqid sseqid sallacc length pident mismatch gapopen qstart qend sstart send evalue bitscore qseq" -max_target_seqs 10 -evalue 1e-40 && sed -i '1i qseqid,sseqid,sallacc,length,pident,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore,qseq' {↑で入力した{出力したいファイル名}}
       2. 例：
      	```
     		$blastn -db genebank_animallia_mito.nc -query st1Representative_seq.fasta -out st1_kamikochi_max10_BLAST10.csv -outfmt "10 qseqid sseqid sallacc length pident mismatch gapopen qstart qend sstart send evalue bitscore qseq" -max_target_seqs 10 -evalue 1e-40 && sed -i '1i qseqid,sseqid,sallacc,length,pident,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore,qseq' st1_kamikochi_max10_BLAST10.csv
  - (余談)今回のプログラムでは、qseqid、sallacc、pident、qseqのみしか使わないので、他の値も確認したいとかでなければ、blastnコマンドは以下でもいいかもです(結局他の値も取った方が他の用途にも使えそうなので、3-1のコマンドでいいと思います...速さも変わらないと思うので...)
    ```
  	$blastn -db genebank_animallia_mito.nc -query st1Representative_seq.fasta -out st1_kamikochi_max10_BLAST10.csv -outfmt "10 qseqid sallacc pident qseq" -max_target_seqs 10 -evalue 1e-40 && sed -i '1i qseqid,sallacc,pident,qseq' st1_kamikochi_max10_BLAST10.csv
---
## プログラムを使う
- 適当な場所にファイルを展開
- `name_taxonomy_create_tree`に移動
  - コンソールのpathが`~/fuga/hoge/name_taxonomy_create_tree$`みたいになっていればOKです
- 事前に用意したCSVファイルを`name_taxonomy_create_tree`下に移動
- FastTreeとMAFFTが必要なのでいれる。↓で入ると思う
  - sudo apt install fasttee
  - sudo apt install mafft
- コマンド入力
  - FASTAファイルとCSVファイル出力だけしたい場合
    - `$python name_taxonomy_create_tree.py {入力するCSVファイル名} {出力したいFASTA、CSVファイル名}`
    - 例：`python name_taxonomy_create_tree.py st1_kamikochi_max10_BLAST10.csv output1`
  - Newickファイル(樹形データが欲しい場合)
    - `$python name_taxonomy_create_tree.py {入力するCSVファイル名} {出力したいFASTA、CSVファイル名} --tree {各種オプション}`
    - 例：`python name_taxonomy_create_tree.py st1_kamikochi_max10_BLAST10.csv output1 --tree --method NJ --bootstrap 250`
      - `--method`：系統樹作成手法の選択。デフォルトでは`ML`になっています。`ML`と`NJ`が選べます
      - `--bootstrap`：ブートストラップ値の指定。デフォルトでは`0`になっています。
      - `--gamma`：ガンマ分布を適用するかどうか。デフォルトでは`False`になっています。
      - `--outgroup`：外群の指定。デフォルトでは`None`になっています。使用したことがないのでわかりませんが、加工後のOTU名を参照することになりそうなので、このオプション使う場合は、FASTAのOTU名を参考することになるかもしれません。
  - 出力されるFASTAファイルとCSVファイル、Newickファイルの拡張子以前は同じ名前になります
    - MAFFTでアライメント後のFASTAファイル名は`_aligned`がつきます
# 出力
- ※`qseqid`が同じlocalBLAST結果は、`qseqid`が高いものだけ使われます。
- FASTA
  - 例：output1.fasta
  - OTU名は`{qseqid}_{sallacc}_{taxonomic_name}_{pident}`→`{localoBLAST時のqueryid}_{accession ID}_{分類群名}_{localoBLASTの一致率}`
  - 配列は、それぞれ検索時の配列です
  - `qseqid`は、変換前の配列との整合性チェックとOTUの一意性確保のためにつけてます
  - 分類群名は一致率`pident`を参照して以下のように場合分け
    1. 98.00 <= pident ⇒ 種小名取得(亜種名があれば含める)
    2. 98.00>pident>=95.00 ⇒ 属名取得
    3. 95.00>pident>=90.00 ⇒ 科名取得
    4. 90.00>pident>=85.00 ⇒ 目名取得
    5. pident<85.00 ⇒ 取得しない
	- 属以上は分類群に欠損データあったら、うまく機能しなさそうなので、運用してみて不都合出たら教えてください...
- CSV
  - 例：output1.csv
  - カラム名は、`qseqid, sallacc, taxonomic_name, pident, qseq`になっています
  - エクセル管理用
- MAFFT後のFASTA
  - ↑で出力されるFASTAと内容はほぼ同じです
    - アライメントされているかされていないかの違い
- FastTreeのNewick
  - 指定したオプション下で作成された樹形ファイルが出ます