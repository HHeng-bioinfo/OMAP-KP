# OMAP-KP
**OMAP-KP** is an R script that performs ordered mapping and assignment for plasmid identification in NGS data from KP.  
## Dependencies
+ To run OMAP-KP: 
  [**BLAST**](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html) [>= 2.13]
  [**R**](https://www.r-project.org/) [>= 3.6.1]
  [**data.table**](https://cran.r-project.org/web/packages/data.table/index.html) [>=1.13]
  [**IRanges**](https://bioconductor.org/packages/release/bioc/html/IRanges.html) [>=2.36]    

## Usage 
### Prepare your Input
1. BLAST reference
   The plasmids.repr.fas from KleTy can be found at https://github.com/zheminzhou/KleTy <br />
   The mobileOG-db_beatrix-1.6.faa can be found at https://mobileogdb.flsi.cloud.vt.edu/entries/database_download
    ```bash
    makeblastdb -in  ./dir/to/ref/plasmids.repr.fas  -dbtype nucl  -parse_seqids -out ./dir/to/ref/plasmids.repr.fas
    ```
    ```bash
    makeblastdb -in ./dir/to/ref/mobileog_target.faa -dbtype prot  -parse_seqids -out ./dir/to/ref/mobileog_target
    ```
2. Cluster and MGE annotation file
   Available here in the repository ("./KleTy/plasmids.repr.clu" and "./mobileOG-db/mge_target.txt")
### Run OMAP-KP
1. To run the script:
    ```bash
    Rscript OMAP_v2.R  ./dir/to/your/input_genome.fasta  ./dir/to/your/output.dir/  ./KleTy/plasmids.repr.fas ./mobileOG-db/mobileog_target ./KleTy/plasmids.repr.clu  ./mobileOG-db/mge_target.txt
    ```

## Output
1. [filename]_pls.tsv
    A file contains: filename, pls_id, cluster, scovs, and slen.
2. [filename]_raw.tsv
    A file contains the raw data of hc_dt.tsv file, including the names of contigs assigned to the plasmid.


## Supplementary files 
The validation dataset suggested a recall rate of 85.27% for plasmids exceeding 10,000 bp in 56 draft genomes. (Updated in 2025/01/17) <br />
https://github.com/HHeng-bioinfo/OMAP-KP_seq

