# OMAP-KP
**OMAP-KP** is an R script that performs ordered mapping and assignment for plasmid identification in NGS data from KP.  
## Dependencies
+ To prepare OMAP-KP's inputs:
  [**BLAST**](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html) [>= 2.13]
  [**PLASME**](https://github.com/HubertTang/PLASMe) [>= 1.1]

+ To run OMAP-KP:
  [**R**](https://www.r-project.org/) [>= 3.6.1]
  [**data.table**](https://cran.r-project.org/web/packages/data.table/index.html) [>=1.13]
  [**IRanges**](https://bioconductor.org/packages/release/bioc/html/IRanges.html) [>=2.36]    

## Sample Run 
The required demo files can be found in omap_test_data
```bash
Rscript OMAP-KP_v0.1.2.R omap_test_data/kp_demo_blast.out  omap_test_data/  omap_test_data/plasme_out_report.csv OMAP-KP_clusters_21.txt 
```

## Usage 
### Prepare your Input
1. PLSDB Blast input
    The plsdb.fna from PLSDB v2021_06_23 can be found in https://figshare.com/articles/dataset/PLSDB_2021_06_23_v2/19153574.
    ```bash
    makeblastdb -in  ./dir/to/ref/plsdb.fna  -dbtype nucl  -parse_seqids -out ./dir/to/ref/plsdb_ref
    ```
    ```bash
    blastn -query ./dir/to/your/genome.fna  -db ./dir/to/ref/plsdb_ref  -evalue 1e-10  -num_threads 4 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qcovhsp slen qlen' -dust no -soft_masking false  > ./dir/to/your/genome.blast.out
    ```
2. PLASME input
   Please install PLASMe and follow its instruction
    ```bash
    python PLASMe.py ./dir/to/your/genome.fna  ./dir/to/your/plasme_out 
    ```
4. Cluster annotation file
   Available here in the repository (OMAP-KP_clusters_21.txt)
### Run OMAP-KP
1. To run the script:
    ```bash
    Rscript OMAP-KP_v0.1.2.R ./dir/to/your/genome.blast.out  ./dir/to/your/output.dir/ ./dir/to/your/plasme_out_report.csv  OMAP-KP_clusters_21.txt 
    ```

## Output
1. [filename]_hc_dt.tsv
    A file contains: filename, pls_id, cluster, scovs, and slen.
2. [filename]_hc_raw.tsv
    A file contains the raw data of hc_dt.tsv file, including the names of contigs that are assigned to the plasmid.

