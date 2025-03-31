# Preparing the input files

   - Use the CSV file of the search results output by `localBLAST` or `BLAST+` as the input file.
   ## **Example of localBLAST execution**
   #### Creating a database from NCBI data
   1. Prepare the processed sequence data using metagenomic analysis tools(e.g. [qiime2](https://qiime2.org/), [USEARCH](https://www.drive5.com/usearch/)...).
   2. Search for the sequence file you want to use as a database and click the "send to" button to output it in FASTA format.
      - Example: Search for mtDNA 16S rRNA regions of animal phyla(https://www.ncbi.nlm.nih.gov/nuccore/advanced)
          ```
          ((((Animalia) AND 16S) NOT whole genome) NOT chromosome) NOT complete genome
          ```

        <details><summary>The location of the "send to" link</summary>

          ![image](https://github.com/user-attachments/assets/7424ed3d-86ba-4afd-96b1-dec875544b98)

         </details>

   3. Create a database using `localBLAST` or `BLAST+`
      ```
      makeblastdb -in {FASTA file to be treated as DB} -dbtype nucl -out {Any DB name}.nc -hash_index -parse_seqids
      ```
       - Example: Creating a database of mtDNA 16S rRNA regions of animal phyla
        ```
        makeblastdb -in animalia_16S.fasta -dbtype nucl -out animalia_16S_db.nc -hash_index -parse_seqids
        ```
   4. Run a BLAST search against the created DB using the FASTA file of the sequence you want to perform a BLAST search on.
         ```
         blastn -db {Search target database name}.nc -query {The name of the FASTA file containing the sequence you want to perform a BLAST search on.} -out {Output File Name}.csv -outfmt "10 qseqid sallacc pident qseq" -max_target_seqs 10 -evalue 1e-40 && sed -i '1i qseqid,sallacc,pident,qseq' {Output File Name}.csv
         ```
         - ※The `sed` command is a Linux command, so be careful if you are using another OS.
         - Example: To perform a BLAST search on the database created in `2.`, use the FASTA file `query_sequence.fasta`.
         ```
         blastn -db animalia_16S_db.nc -query query_sequence.fasta -out output_quried.csv -outfmt "10 qseqid sallacc pident qseq" -max_target_seqs 10 -evalue 1e-40 && sed -i '1i qseqid,sallacc,pident,qseq' output_quried.csv
         ```
         - Please make sure the order of the items in `-outfmt "10 xx yy"` and `sed -i '1i xx,yy'` is the same.
         - If the `&& sed command` does not work, you can manually name the columns in the same format as the input file example below.
         - You can increase the arguments of `-outfmt "10 xx yy zz"` as needed, but at the very least, it will work with `qseqid`, `sallacc`, `pident` and `qseq`.

   ### Example input file
   CSV format
   | qseqid | sallacc | pident | qseq |
   | ---- | ---- | ---- | ---- |
   | 9534cfe94fa593ed71 | AB1234 | 98.805 | GATCGAT・・・ |
   | 9534cfe94fa593ed72 | AB2345 | 96.016 | GATCGAT・・・ |
   | 9534cfe94fa593ed73 | AB3456 | 96.032 | GATCGAT・・・ |
   | 9534cfe94fa593ed74 | AB4567 | 96.032 | GATCGAT・・・ |

    <details><summary>About each column</summary>

      - `qseqid`
        - This is a sequence number given when performing a BLAST search. It is assigned to each sample by the sequence number in the query.
        - There is an option to select BLAST results with high `pident` to include in the dataset based on this value.
      - `sallacc`
        - If you use NCBI data in your database, this corresponds to `Accession ID`.
        - Used to search for species names
      - `pident`
        - Match rate between sample sequence and reference sequence output when performing BLAST search
      - `qseq`
        - Nucleotide sequences used in BLAST searches
    </details>

   - It is assumed that the files obtained after searching with `localBLAST`, `BLAST+`, etc. will be used, but as long as they match the above format, there is no problem in creating them manually.
