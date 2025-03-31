# 入力ファイルの準備

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
