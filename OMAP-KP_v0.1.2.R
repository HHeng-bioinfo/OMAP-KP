#!/usr/bin/env Rscript

## Info
# v0.01
# $1: blast file
# $2: out dir
# $3: PLASMe file
# $4: mob cluster file
# v0.1.2 2024/05/29



##### Sec 1 library and load #####
## Timer
#
start_time <- Sys.time()
print(paste0("Start at:", start_time))


## args in
# load
args = commandArgs(trailingOnly=TRUE)


## library
#
suppressPackageStartupMessages({
  library(data.table, quietly=TRUE, verbose=FALSE)
  library(IRanges, quietly=TRUE, verbose=FALSE)
})


## Load files
#
print("Read files...")

#
pls_blast <- fread(args[1])
outdir <- args[2]
ctg_included <- fread(args[3], 
                      col.names = c("contig", "length", "pls_id", "order", 
                                    "evidence", "score", "amb_region"))
mob_clu <- fread(args[4])

#
file_path <- args[1]
setDTthreads(4)
dir.create(outdir, recursive = T)

#
if (nrow(pls_blast) == 0) {stop("Empty blast file") 
} else {
  print(paste0("Read files in: ", file_path, " Done")) 
}
colnames(pls_blast) <- c("qseqid", "sseqid", "pident", "length", "mismatch", 
                         "gapopen", "qstart", "qend", "sstart", "send", 
                         "evalue", "bitscore", "qcovs", "qcovhsp", 
                         "slen", "qlen")

#
filename <- sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(file_path))
pls_blast[, filename := filename]
setcolorder(pls_blast, c("filename", "qseqid", "sseqid", "pident", "length", 
                         "mismatch", "gapopen", "qstart", "qend", "sstart", 
                         "send", "evalue", "bitscore", "qcovs", "qcovhsp", 
                         "slen", "qlen"))



##### Sec 2 Filter and Fun #####
## Function SCOVS
#
red_len <- function(start, end) {
  rd <- reduce(IRanges(start, end))
  return(sum(rd@width))
}


## merge clu
#
mob_clu <- mob_clu[, .(pls_id = sample_id, cluster = primary_cluster_id, pls_length = size)]
mob_clu[mob_clu == "-"] <- ""


## modi blast
#
pls_blast[, pls_id := gsub("\\|$", "", gsub("^[a-zA-Z]*\\|", "", sseqid))]

#
if (nrow(pls_blast[!(pls_id %in% mob_clu$pls_id)]) > 0) {
  print("opps")
  fwrite(merge.data.table(pls_blast, mob_clu, all.x = T)[is.na(cluster)],
         paste0(outdir, "/", "error.tsv"), sep = '\t')
  stop("Missing pls ID")
} 
pls_blast <- merge.data.table(pls_blast, mob_clu, by = "pls_id")

#
pls_blast[, sseqid := NULL]
pls_blast[, evalue := NULL]
pls_blast[, bitscore := NULL]
pls_blast[, qcovhsp := NULL]


## filter blast
#
pls_phit <- copy(pls_blast)
pls_phit <- pls_phit[pident >= 85][length >= 1000]


#
pls_blast <- pls_blast[qcovs >= 85][pident >= 85][length >= 1000]
pls_blast <- pls_blast[qseqid %in% ctg_included[(evidence == "BLASTn") | (score >= 0.85)]$contig]



##### Sec 3 HC plasmids #####
#
part_forward <- copy(pls_blast)
part_forward[sstart > send, c("sstart", "send") := .(send, sstart)]

#
part_map <- part_forward[, red_len(sstart, send), by = .(filename, pls_id, pls_length, cluster, slen)][, .(filename, pls_id, cluster, pls_length, scovs = V1/slen)]
part_map <- part_map[scovs >= 0.85][order(filename, -pls_length)]
part_map_set <- part_map[part_map[ , .I[which.max(scovs)], by = .(filename, cluster)]$V1][order(filename, -pls_length)]

##
#
indt <- copy(part_map_set)
target_out <- data.table(filename = NULL, pls_id = NULL, cluster = NULL, scovs = NULL, slen = NULL)
sub_out <- data.table(filename = NULL, pls_id = NULL, cluster = NULL, scovs = NULL, slen = NULL, pp = NULL, pc = NULL)
part_out_in <- merge(part_forward, part_map_set[, .(filename, pls_id)], by = c("filename", "pls_id"))

#
if (nrow(indt) == 0) {
  print("No Any HC plasmids")
} else {
  for (lineN in 1:nrow(indt)) {
    ##
    #
    filename_in <- indt[lineN, filename]
    cluster_in <- indt[lineN, cluster]
    pls_id_in <- indt[lineN, pls_id]
    print(paste0("Check: ", filename_in, " ", cluster_in))
    
    #
    target_per <- part_out_in[filename == filename_in][pls_id == pls_id_in][, red_len(sstart, send), by = .(filename, pls_id, cluster, slen)][, .(filename, pls_id, cluster, scovs = V1/slen, slen)]
    
    
    ##
    if (nrow(target_per) == 0) {
      print(paste0("Removed: ", filename_in, " ", pls_id_in))
    } else if (target_per$scovs >= 0.85) {
      ##
      #
      target_out <- rbind(target_out, target_per)
      
      #
      print(paste0("Keep: ", filename_in, " ", pls_id_in))
      
      
      ##
      #
      sub_part_in <- part_out_in[(qseqid %in% part_out_in[filename == filename_in][pls_id == pls_id_in, qseqid])][filename == filename_in][pls_id != pls_id_in]
      sub_part_out <- sub_part_in[, red_len(sstart, send), by = .(filename, pls_id, cluster, slen)][, .(filename, pls_id, cluster, scovs = V1/slen, slen)][scovs >= 0.85]
      
      #
      if (nrow(sub_part_out) != 0) {
        #
        part_out_blast <- part_out_in[filename == filename_in][pls_id == pls_id_in]
        part_out_blast_um <- unique(part_out_blast[, .(filename, pls_id, qseqid, qcovs, qlen)])[, .(filename, pls_id, qseqid, um_p = (1-qcovs/100)*qlen)]
        
        #
        sub_part_out_blast <- merge(sub_part_out[, .(filename, pls_id, scovs)], sub_part_in, by = c("filename", "pls_id"))
        sub_part_out_blast_um <- unique(sub_part_out_blast[, .(filename, pls_id, qseqid, qcovs, qlen)])[, .(filename, pls_id, qseqid, um_s = (1-qcovs/100)*qlen)]
        
        #
        for (subpls in unique(sub_part_out_blast_um$pls_id)) {
          print(paste0("check sub pls:", subpls))
          sub_um <- merge(sub_part_out_blast_um[filename == filename_in][pls_id == subpls], part_out_blast_um[, .(qseqid, um_p)], by = c("qseqid"))
          if (sub_um[, sum(um_s)] > sub_um[, sum(um_p)]) {
            print(paste0("Removed sub pls: ", filename_in, " ", subpls))
          } else {
            print(paste0("Keep sub pls: ", filename_in, " ", subpls))
            sub_out <- rbind(sub_out, sub_part_out[filename == filename_in][pls_id == subpls][, .(filename, pls_id, cluster, scovs, slen, pp = pls_id_in, pc = cluster_in)])
          }
        }
      }
      
      
      ##
      #
      part_out_in <- part_out_in[!(qseqid %in% part_out_in[filename == filename_in][pls_id == pls_id_in, qseqid])]
    } else {
      print(paste0("Removed: ", filename_in, " ", pls_id_in))
    }
  }
}


## write
# remove sub
if (nrow(sub_out) != 0) {
  target_out <- target_out[!cluster %in% sub_out[, .N, by = c("filename", "pc")][N > 1, pc]]
} else {
  print(paste0("No sub plasmids: ", filename))
}

# write HC
if (nrow(target_out) != 0) {
  #
  fwrite(target_out, paste0(outdir, "/", filename, "_hc_dt.tsv"), sep = '\t', eol = '\n')
  #
  fwrite(merge(target_out, pls_blast, by = c("filename", "pls_id", "cluster", "slen")), paste0(outdir, "/", filename, "_hc_raw.tsv"), sep = '\t', eol = '\n')
} else {
  print(paste0("No HC plasmids: ", filename))
}



##
#
print(paste0("Proc files: ", filename, " Done"))
end_time <- Sys.time()
print(paste0("End at:", end_time))

#
fly_time <- end_time - start_time
print(paste0("Used time: ", round(fly_time), "s"))


