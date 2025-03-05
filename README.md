# MICUM(Moleculer Identification pipeline Computational Unit Manager)
- [README(日本語)](README-jp.md)
## What`s app
- This app is a pipeline tool that aims to perform phylogenetic analysis and species determination simply and almost automatically.
- It uses Docker to configure the virtual environment, so it is largely independent of the user's OS or environment.
### Concept
- Designed for phylogenetic analysis using environmental DNA and more ...

### What we can do
- Create a FASTA file with OTUs including species names based on the file output by `localBLAST` or `BLAST+` (only for DBs composed of NCBI data)
- Sequence alignment([MAFFT](https://mafft.cbrc.jp/alignment/server/index.html))
- Haplotype detection([VSEARCH](https://github.com/torognes/vsearch))
- Creating a phylogenetic tree([FASTTree](https://morgannprice.github.io/fasttree/))
- Species classification analysis([bPTP](https://species.h-its.org/), [mPTP](https://github.com/Pas-Kapli/mptp))

## Installation
1. Installing Docker
     - https://docs.docker.jp/engine/getstarted/step_one.html
2. Download or clone the GitHub repository
     - DL
       <details><summary>DL link</summary>

       ![image](https://github.com/user-attachments/assets/ad10015a-dbe1-4498-a751-ae2e0c42a47a)

       - Download from the `Download ZIP` button
      </details>

     - Clone
      ```
      git clone https://github.com/shigebio/MICUM
      ```
3. Go to DL or cloned file
    ```
    cd /{path to MICUM}/MICUM
    ```
4. Creating a virtual environment
Depending on your environment, you may need to use `sudo` before the `docker` command.
    1. Launch Docker Desktop or Docker Engine
         ```
          # If you use Docker Engine
          service docker start
          # Check that docker is running
          docker version

          # If Docker Compose is not available (Linux)
          sudo apt update
          sudo apt install docker-compose
         ```

   1. Creating and starting a virtual environment
        ```
        docker-compose up -d
        ```
     - Building a virtual environment takes time.
     - You can also get the image from Docker hub (this is not necessary if you have followed the steps above <b>※not up to date</b>).
       - [shigebio/name_taxonomy_create_tree](https://hub.docker.com/r/shigebio/name_taxonomy_create_tree)
     - Please make sure that the `input` and `output` folders have been created under the `app` folder of the downloaded or cloned file.

### How to update the version

  ```
  # Go to the path where update.sh is located
  cd /path/to/update.sh

  # Update
  bash update.sh
  ```

## How to Use
1. Move the CSV file you prepared in advance to the `input` folder(See [here](#input_section) for how to prepare your CSV file.).
   - Direct specification is also possible
1. Launch Docker Desktop or Docker Engine
      ```
      # If you use Docker Engine
      service docker start
      ```

1. Starting a virtual environment
    ```
    docker-compose up -d
    ```
2. Enter the virtual environment
    ```
    docker exec -it micum /bin/bash
    ```
    - If it looks like this, you are in a virtual environment.
    ```
    root@absd1234:/app#
    ```
3. Executing commands
   1. Basic commands
      - `python3 MICUM.py {Input CSV file name} --tree {Options}`
        - Example: `python3 MICUM.py your_input.csv --tree --method -ml --bootstrap 250`
           <details><summary>Options</summary>

            - `--top` : You can specify 1 to 10 of the top `pident` with the same `qseqid`.
              - Default: `1`
            - `--o` : Output name of the output file. If not present, the default name is applied.
            - `--class` : Analysis can be limited to specific taxa (currently only classes).
              - Please note that if taxonomic group information obtained from NCBI or GBIF is not entered, it will not be searchable.
              - You can also manually create a `class` column and enter any value.
            - `--tree`
               ```
               ## Sub options
               # --method : You can select the method for constructing a phylogenetic tree. Only `-ml` can be used, which will use the ML method.
               # --bootstrap : Number of bootstrap iterations. Default: `250`
               # --gamma : Whether to apply gamma distribution. Default: `False`
               # --outgroup {OTU name} : Specify the outgroup
            - `--onlyp`： Perform phylogenetic analysis and beyond (alignment → removal of identical haplotypes → creation of phylogenetic tree → species determination analysis)
               - https://github.com/shigebio/MICUM/pull/7
             - `--bptp` : bPTP Analysis Options
               ```
               ## Sub options
               # --mcmc : You can set the number of MCMC iterations. Only `-ml` can be used, which will use the ML method.
               # --thinning : You can set the number of samplings. Sampling is performed for each MCMC chain every specified number. Default: `100`
               # --burnin : Burn-in ratio (0.1~1.0). Default: `0.1`
               # --seed : The seed you want to assign. The same seed value will give you the same results every time for the same input. The default is a random seed.
           </details>
   2. If you only want to output FASTA and CSV files
      - `python3 MICUM.py {input CSV file name} {output file name}`
        - Example: `python3 MICUM.py your_data.csv output`

4. Stopping a container
   - `sudo docker-compose down`
     - If you keep running the container, it will consume memory, so it is better to stop it.
     - To reboot, see `4. Building a virtual environment > Each OS > 2. Starting the virtual environment`
    
<div id="input_section"></div>
<details><summary>Preparing the input files</summary>

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
</details>

## Outputs
### FASTA file/CSV file of OTUs with assigned species names
- Without the `--class` option
  - `pre_filtered_{output file name}.csv`
  - `pre_filtered_{output file name}.festa`
- With the `--class` option
  - `filtered_{output file name}.csv`
  - `filtered_{output file name}.festa`

<b>＊CAUTION＊
     Due to the structure of the NCBI database, the taxon information columns may be out of sync when retrieved from the Entrez API (NCBI). Please check the taxon information in the file.
</b>

### MAFFT aligned file
- `{input/output file name}_aligned.fasta`

### File after removal of identical haplotypes by VSEARCH
- `{output file name}_vsearch.fasta`
- This file is used for phylogenetic tree construction and species classification analysis.

### Results of bPTP analysis
- The `bPTP_{folder created date}` folder will be created.
  <details><summary>A file containing the main results</summary>

    - `output_base_tree_bptp_{output file name}.txt.PTPhSupportPartition.txt`
      - Analysis results by simple heuristic search. Text format.
    - `output_base_tree_bptp_{output file name}.txt.PTPhSupportPartition.txt.png`
      - Analysis results by simple heuristic search. Image (png) format.
    - `output_base_tree_bptp_{output file name}.txt.PTPhSupportPartition.txt.svg`
      - Analysis results by simple heuristic search. Image (svg) format.
    - `output_base_tree_bptp_output.txt.PTPMLPartition.txt`
      - Analysis results using the ML method. Text format.
    - `output_base_tree_bptp_output.txt.PTPMLPartition.txt.png`
      - Analysis results using the ML method. Image (png) format.
    - `output_base_tree_bptp_output.txt.PTPMLPartition.txt.svg`
      - Analysis results using the ML method. Image (svg) format.

  </details>

### Results of mPTP analysis
- The `mPTP_{folder created date}` folder will be created.
  <details><summary>A file containing the main results</summary>

    - `output_base_tree_mptp_{output file name}.txt.txt`
      - Analysis results using the ML method. Text format.
    - `output_base_tree_mptp_{output file name}.txt.svg`
      - Analysis results using the ML method. Image (svg) format.
  </details>

# The problems you have experienced
→Please [create a new issue](https://github.com/shigebio/MICUM/issues) and include the details🙏
