# OMAP-KP
**OMAP-KP** is a R script that perform ordered mapping and assignment for plasmids identification in NGS data from KP.
  
  
## Dependencies
To prepare OMAP-KP's inputs:
[**BLAST**](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html) [>= 2.13]
[**PLASME**](https://github.com/HubertTang/PLASMe) [>= 1.1]
  
To run OMAP-KP:
[**R**](https://www.r-project.org/) [>= 3.6.1]
[**data.table**](https://cran.r-project.org/web/packages/data.table/index.html) [>=1.13]
[**IRanges**](https://bioconductor.org/packages/release/bioc/html/IRanges.html) [>=2.36]
  
   
## Sample Run 
```Shell
Rscript listtable/OMAP-KP_v0.1.2.R ./test/genome.blast.out  ./test/ ./test/plasme_out_report.csv  OMAP-KP_clusters.txt 
```
  
   
## Usage 
### Prepare your the Input
1. PLSDB Blast input
    The plsdb.fna from PLSDB v2021_06_23 can be found in https://figshare.com/articles/dataset/PLSDB_2021_06_23_v2/19153574.
    ```Shell
    makeblastdb -in  ./dir/to/ref/plsdb.fna  -dbtype nucl  -parse_seqids -out ./dir/to/ref/plsdb_ref
    blastn -query ./dir/to/your/genome.fna  -db ./dir/to/ref/plsdb_ref  -evalue 1e-10  -num_threads 4 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qcovhsp slen qlen' -dust no -soft_masking false | awk -v filename="your_file_name" '{print filename"\t"$0}' > ./dir/to/your/genome.blast.out
    ```
2. PLASME input
    ```Shell
    python PLASMe.py ./dir/to/your/genome.fna  ./dir/to/your/plasme_out 
    ```
  
### Run OMAP-KP
1. To run the script:
    ```Shell
    Rscript listtable/Pls_01_extract_fus.R ./dir/to/your/genome.blast.out  ./dir/to/your/output.dir/ ./dir/to/your/plasme_out_report.csv  OMAP-KP_clusters.txt 
    ```
  
  
## Output
1. [filename]_hc_dt.tsv
    A file contains: 
2. [filename]_hc_raw.tsv
    A file contains: 

