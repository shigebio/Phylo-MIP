# Phylo-MIP(Phylogeny-based Molecular Identification Pipeline)
- [README(日本語)](README-jp.md)
## What`s app
- This app is a pipeline tool that aims to perform phylogenetic analysis and species determination simply and almost automatically.
- It uses Docker to configure the virtual environment, so it is largely independent of the user's OS or environment.
### Concept
Designed for phylogenetic analysis using environmental DNA and more ...

### What we can do
- Create a FASTA file with OTUs including species names based on the file output by `localBLAST` or `BLAST+` (only for DBs composed of NCBI data) : <b>Phylo-MIP.py</b>
  - Sequence alignment([MAFFT](https://mafft.cbrc.jp/alignment/server/index.html))
  - Haplotype detection([VSEARCH](https://github.com/torognes/vsearch))
  - Creating a phylogenetic tree([FASTTree](https://morgannprice.github.io/fasttree/))
  - Species classification analysis([bPTP](https://species.h-its.org/), [mPTP](https://github.com/Pas-Kapli/mptp))
- Merging the output file of Qiime2 and the output file of Phylo-MIP pipline : <b>merge_data.py</b>
- Confirmed that it can be run on MAcOS or LinuxOS, but it can also be run on WindowsOS using WSL.

## Installation
1. Installing Docker
     - https://docs.docker.jp/engine/getstarted/step_one.html
2. Download or clone Phylo-MIP
  Phylo-MIP is available as a direct download or by cloning the repository.
    <details><summary>DL link</summary>

    ![image](https://github.com/user-attachments/assets/ad10015a-dbe1-4498-a751-ae2e0c42a47a)

    - Download from the `Download ZIP` button
    </details>

    Clone
    ```
    git clone https://github.com/shigebio/Phylo-MIP
    ```

3. Go to DL or cloned file
    ```
    cd /{path to Phylo-MIP}/Phylo-MIP
    ```
4. Creating a virtual environment
    ```
    chmod +x setup.sh
    chmod +x entrypoint.sh
    ./setup.sh
    ```
---
## How to Use Phylo-MIP
**Before you run**
Phylo-MIP.py uses NCBI and GBIF APIs. Please refer to the [NCBI](https://blast.ncbi.nlm.nih.gov/doc/blast-help/developerinfo.html#developerinfo) and [GBIF](https://techdocs.gbif.org/en/openapi/v1/species) API guidelines.
Before running for the first time, as much as possible, please change the example email address in [this code](https://github.com/shigebio/Phylo-MIP/blob/main/app/Phylo-MIP.py#L22) to your own.

Part of the NCBI guidelines
>- Do not contact the server more often than once every 10 seconds.
>- Do not poll for any single RID more often than once a minute.
>- Use the URL parameter email and tool, so that the NCBI can contact you if there is a problem.
>- Run scripts weekends or between 9 pm and 5 am Eastern time on weekdays if more than 50 searches will be submitted.
---

1. Executing commands

    **Basic commands**
      ```
      phylo-mip {Path to Input CSV} --tree {Options}
      ```
    **Example**
      ```
      phylo-mip ./path/to/your_input.csv --tree --method ML --bootstrap 250
      ```

   See [here](https://github.com/shigebio/Phylo-MIP/blob/main/README-Preparing_the_input_files.md) for acceptable input file formats

      <details><summary>See Options</summary>

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

   **If you only want to output FASTA and CSV files**
      ```
      phylo-mip {input CSV file name} {output file name}
      ```
      **Example**
      ```
      phylo-mip your_data.csv output
      ```

## Outputs
```
=== Directory Structure Verification ===
micum_output_{timestamp}/
  taxonomy/ # CSV, FASTA files with taxonomic information
    taxonomic_data.csv # Final output file integrating PTP analysis results and taxon information
    taxonomic_sequences.fasta
  phylogeny/ # Phylogenetic tree by fattree
    {timestamp}_ML_phylogenetic_tree.nex # The file will be prefixed with the execution time
    {timestamp}_ML_{execution time}_ML_phylogenetic_tree.nwk
  alignment/ # File after mafft and VSEARCH
    {timestamp}_haplotype_clusters.tsv # TSV file of haplotype identification results by VSEARCH
    {timestamp}_haplotype_clusters.csv　# CSV file of haplotype identification results by VSEARCH
    {timestamp}_clustered_sequences.fasta # File after removal of identical haplotypes by VSEARCH
    {timestamp}_aligned_sequences.fasta # File after alignment by MAFFT
  mptp/ # mptp output result
    {timestamp}_mPTP_analysis/
    {timestamp}_mPTP_species_delimitation.txt # TXT file containing assigned species and likelihood of each species
    {timestamp}_mPTP_species_delimitation.svg # Tree file (SVG)
  bptp/ # bptp output result
    {timestamp}_bPTP_analysis/
    {timestamp}_bPTP_species_delimitation.PTPhSupportPartition.txt.svg # Tree file (SVG) constructed by simple heuristic search
    {timestamp}_bPTP_species_delimitation.llh.pdf # Trace log graph of MCMC chain (PDF)
    {timestamp}_bPTP_species_delimitation.PTPhSupportPartition.txt.sh.tre # Tree file (TREE) constructed by simple heuristic search
    {timestamp}_bPTP_species_delimitation.PTPMLPartition.txt # TXT file containing species assigned by ML method and likelihood of each species
    {timestamp}_bPTP_species_delimitation.PTPMLPartition.txt.ml.tre # Tree file (TREE) constructed by ML method
    {timestamp}_bPTP_species_delimitation.PTPMLPartition.txt.png # Tree file (PNG) constructed by ML method
    {timestamp}_bPTP_species_delimitation.PTPllh.txt # MCMC chain trace log
    {timestamp}_bPTP_species_delimitation.PTPhSupportPartition.txt.png # Tree file (PNG) constructed by simple heuristic search
    {timestamp}_bPTP_species_delimitation.PTPPartitions.txt
    {execution time}_bPTP_species_delimitation.PTPhSupportPartition.txt # TXT file containing species assigned by simple heuristic search and likelihood of each species
    {execution time}_bPTP_species_delimitation.PTPMLPartition.txt.svg # Tree file (SVG) constructed by ML method
    {execution time}_bPTP_species_delimitation.PTPPartitonSummary.txt # Species divided by each method
```

The taxonomic status of sequences assigned to an OTU is assigned according to the value of pident (match rate of search sequence in localBLAST search) as follows:
```
pident >= 98.00 : Species name
95.00 <= pident < 98.00 : Genus name
92.00 <= pident < 95.00 : Family name
85.00 <= pident < 92.00 : Order name
```
code: https://github.com/shigebio/Phylo-MIP/blob/main/app/Phylo-MIP.py#L291-L298

<b>＊CAUTION＊
     Due to the structure of the NCBI database, the taxon information columns may be out of sync when retrieved from the Entrez API (NCBI). Please check the taxon information in the file.
</b>

---
## How to Use merge_data
**Basic commands**
  ```
  merge_data -q {Qiime output file path} -p {Phylo-MIP pipeline output file path} -f {The file format you want to output: csv/tsv} -o {output file name}
  ```

Output file name default: time plefix on executed.

The merged file will be output to woking directory.


# Citation
Please feel free to use it for personal use, but if you use this pipeline for research or other purposes in publications, we would appreciate it if you could cite the following:

**Phylo-MIP: Phylogeny-based Molecular Identification Pipeline for DNA metabarcoding, and assessment of insect communities in subalpine river ecosystems
Takumi Yshida, Shonosuke Shigeta, Yuta Hasebe, Masaki Takenaka
bioRxiv 2025.11.10.687572; doi: https://doi.org/10.1101/2025.11.10.687572**


# The problems you have experienced
→Please [create a new issue](https://github.com/shigebio/Phylo-MIP/issues) and include the details🙏
