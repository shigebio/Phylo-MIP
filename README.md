# MICUM(Moleculer Identification pipeline Computational Unit Manager)
- [README(日本語)](README-jp.md)
## What`s app
- This app is a pipeline tool that aims to perform phylogenetic analysis and species determination simply and almost automatically.
- It uses Docker to configure the virtual environment, so it is largely independent of the user's OS or environment.
### Concept
- Designed for phylogenetic analysis using environmental DNA and more ...

### What we can do
- Create a FASTA file with OTUs including species names based on the file output by `localBLAST` or `BLAST+` (only for DBs composed of NCBI data) : <b>MICUM.py</b>
  - Sequence alignment([MAFFT](https://mafft.cbrc.jp/alignment/server/index.html))
  - Haplotype detection([VSEARCH](https://github.com/torognes/vsearch))
  - Creating a phylogenetic tree([FASTTree](https://morgannprice.github.io/fasttree/))
  - Species classification analysis([bPTP](https://species.h-its.org/), [mPTP](https://github.com/Pas-Kapli/mptp))
- Merging the output file of Qiime2 and the output file of MICUM pipline : <b>merge_data.py</b>

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
---
## How to Use MICUM.py
**Before you run**
MICUM.py uses NCBI and GBIF APIs. Please refer to the [NCBI](https://blast.ncbi.nlm.nih.gov/doc/blast-help/developerinfo.html#developerinfo) and [GBIF](https://techdocs.gbif.org/en/openapi/v1/species) API guidelines.
Before running for the first time, as much as possible, please change the example email address in [this code](https://github.com/shigebio/MICUM/blob/main/app/MICUM.py#L21) to your own.

Part of the NCBI guidelines
>- Do not contact the server more often than once every 10 seconds.
>- Do not poll for any single RID more often than once a minute.
>- Use the URL parameter email and tool, so that the NCBI can contact you if there is a problem.
>- Run scripts weekends or between 9 pm and 5 am Eastern time on weekdays if more than 50 searches will be submitted.
---
1. Move the CSV file you prepared in advance to the `input` folder(See [here](#input_section) for how to prepare your CSV file.).
   - Direct specification is also possible
2. Launch Docker Desktop or Docker Engine
      ```
      # If you use Docker Engine
      service docker start
      ```

3. Starting a virtual environment
    ```
    docker-compose up -d
    ```
4. Enter the virtual environment
    ```
    docker exec -it micum /bin/bash
    ```
    - If it looks like this, you are in a virtual environment.
    ```
    root@absd1234:/app#
    ```
5. Executing commands
   1. Basic commands
      ```
      python3 MICUM.py {Input CSV file name} --tree {Options}
      ```
      Example
      ```
      python3 MICUM.py your_input.csv --tree --method ML --bootstrap 250
      ```

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
            # --method : You can select the method for constructing a phylogenetic tree. Only `ML` can be used, which will use the ML method.
            # --bootstrap : Number of bootstrap iterations. Default: `250`
            # --gamma : Whether to apply gamma distribution. Default: `False`
            # --outgroup {OTU name} : Specify the outgroup
        - `--onlyp`： Perform phylogenetic analysis and beyond (alignment → removal of identical haplotypes → creation of phylogenetic tree → species determination analysis)
            - https://github.com/shigebio/MICUM/pull/7
          - `--bptp` : bPTP Analysis Options
            ```
            ## Sub options
            # --mcmc : You can set the number of MCMC iterations.
            # --thinning : You can set the number of samplings. Sampling is performed for each MCMC chain every specified number. Default: `100`
            # --burnin : Burn-in ratio (0.1~1.0). Default: `0.1`
            # --seed : The seed you want to assign. The same seed value will give you the same results every time for the same input. The default is a random seed.
        </details>
        <br>
   2. If you only want to output FASTA and CSV files
      ```
      python3 MICUM.py {input CSV file name} {output file name}
      ```
      Example
      ```
      python3 MICUM.py your_data.csv output
      ```

6. Stopping a container
   ```
   sudo docker-compose down
   ```
     - If you keep running the container, it will consume memory, so it is better to stop it.
     - To reboot, see `4. Building a virtual environment > Each OS > 2. Starting the virtual environment`

## Outputs
### FASTA file/CSV file of OTUs with assigned species names
- Without the `--class` option
  - `pre_filtered_{output file name}.csv`
  - `pre_filtered_{output file name}.festa`
- With the `--class` option
  - `filtered_{output file name}.csv`
  - `filtered_{output file name}.festa`

The taxonomic status of sequences assigned to an OTU is assigned according to the value of pident (match rate of search sequence in localBLAST search) as follows:
```
pident >= 98.00 : Species name
95.00 <= pident < 98.00 : Genus name
90.00 <= pident < 95.00 : Family name
85.00 <= pident < 90.00 : Order name
```
https://github.com/shigebio/MICUM/blob/main/app/MICUM.py#L184-L192

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

---
## How to Use merge_data.py
1. Start the Docker container in the same way as when you run MICUM.py, and enter the virtual environment(same virtual environment as when you run MICUM.py).
   - Please follow steps 2 to 4 in [How to Use MICUM.py](https://github.com/shigebio/MICUM?tab=readme-ov-file#how-to-use-micum.py). ※if you did so when running MICUM.py, you do not need to do this.
2. Move the output files of Qiime and MICUM pipeline to the input folder.
3. Executing commands
   - Basic commands
    ```
    python3 merge_data.py -q {Qiime output file name} -m {MICUM pipeline output file name} -f {The file format you want to output: csv/tsv}
    ```
    <details><summary>More Options</summary>

      - `-o` : output file name.
        - Default: time plefix on executed
    </details>

4. The merged file will be output to output file.

# The problems you have experienced
→Please [create a new issue](https://github.com/shigebio/MICUM/issues) and include the details🙏
