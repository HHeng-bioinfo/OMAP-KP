#!/usr/bin/env Rscript

## Info
# Version 2.0.2
# conda activate R_4.3
# HHeng
# v2.0 2025/10/20





##### Sec 1 library and load #####
##### __Sec 1.1 library #####
## library
library(data.table)
library(IRanges)



##### __Sec 1.2 load #####
## Fun
# reduce length
red_len <- function(start, end) {
  rd <- reduce(IRanges(start, end))
  return(sum(rd@width))
}

# calculate overlap
calculate_overlap <- function(start1, end1, start2, end2) {
  overlap <- pmax(0, pmin(end1, end2) - pmax(start1, start2))
  overlap_percentage <- overlap / pmax(end1 - start1, end2 - start2)
  return(overlap_percentage)
}


## Input
#
args = commandArgs(trailingOnly=TRUE)

#
print(paste0("START input: ", args[1], " in ", Sys.time()))



# arg 1 strain
strain_path <- args[1]
strain_name <- sub("\\.[^.]+$", "", basename(strain_path))


# arg 2 output dir
out_dir <- args[2]
dir.create(out_dir, recursive = T)

# arg 3 pls db
pls_db <- args[3]

# arg 4 mge db
mge_db <- args[4]

# arg 5 pls cluster
print(args[5])
pls_clu <- fread(args[5])

# arg 6 target mge table
print(args[6])
mge_table <- fread(args[6])





##### Sec 2 proc input dt #####
##### __Sec 2.1 system blast #####
##
print(paste0("Start blast: ", Sys.time()))

## blast pls
cmd = paste0(" blastn -query ", strain_path, " -db ", pls_db, " -evalue 1e-10 ", 
             " -num_threads 4 ", 
             " -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qcovhsp slen qlen' ",
             " -dust no -soft_masking false -max_target_seqs 100000000",
             " -out ", out_dir, "/", strain_name, "_blast_pls.out")

#
print(cmd)
system(cmd)


## blast mge
cmd = paste0(" blastx -query ", strain_path, " -db ", mge_db, " -evalue 1e-10 ",
             " -num_threads 8 ",
             " -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qcovhsp slen qlen' ",
             " -seg no -soft_masking false -max_target_seqs 100000000",
             " -out ", out_dir, "/", strain_name, "_blast_mge.out")

#
print(cmd)
system(cmd)



##### __Sec 2.2 pls db res #####
## cluster define
#
pls_clu <- pls_clu[, .(sseqid = V1, cluster = V3, slen = V4)]


## pls blast
#
pls_ngs <- fread(paste0(out_dir, "/", strain_name, "_blast_pls.out"), colClasses=list(character=1:2))

#
if (nrow(pls_ngs) == 0) {
  print("Empty pls blast result")
} else {
  #
  blast_header <- c("qseqid", "sseqid", "pident", "length", "mismatch", 
                    "gapopen", "qstart", "qend", "sstart", "send", "evalue", 
                    "bitscore", "qcovs", "qcovhsp", "slen", "qlen")
  colnames(pls_ngs) <- blast_header
  
  #
  pls_ngs[, filename := strain_name]
  setcolorder(pls_ngs, c("filename", blast_header))
}


##### __Sec 2.3 mge db res #####
## mge blast
#
mge_ngs <- fread(paste0(out_dir, "/", strain_name, "_blast_mge.out"), colClasses=list(character=1:2))


#
if (nrow(mge_ngs) == 0) {
  print("Empty mge blast result")
  mge_ngs <- copy(pls_ngs[0])
} else {
  #
  colnames(mge_ngs) <- blast_header
  
  #
  mge_ngs[, filename := strain_name]
  setcolorder(mge_ngs, c("filename", blast_header))
  
  #
  mge_ngs <- mge_ngs[pident >= 0.85][length/slen >= 0.85]
  
  #
  mge_ngs <- mge_ngs[sseqid %in% mge_table$V1]
  
  #
  mge_ngs[sstart > send, c("sstart", "send") := .(send, sstart)]
  
  ## remove duplication
  #
  mge_rd <- copy(mge_ngs[0])
  
  #
  for (file_in in unique(mge_ngs$filename)) {
    for (qseq in mge_ngs[filename == file_in, unique(qseqid)]) {
      # print(qseq)
      dt_in <- mge_ngs[filename == file_in][qseqid == qseq]
      
      for (i in 1:nrow(dt_in)) {
        #
        ol_per <- calculate_overlap(dt_in[i, qstart], dt_in[i, qend], dt_in$qstart, dt_in$qend)
        
        #
        if (any(ol_per[-i] >= 0.85)) {
          out_dt_run <- dt_in[which(ol_per >= 0.85)][order(bitscore, decreasing = T)][1]
        } else {
          out_dt_run <- dt_in[i]
        }
        
        mge_rd <- rbind(mge_rd, out_dt_run)
      }
    }
  }
  
  #
  mge_ngs <- unique(mge_rd)
}



##### Sec 3 Call #####
##### __Sec 3.1 filter #####
##
print(paste0("Start call: ", Sys.time()))


##
#
part_blast <- pls_ngs[qcovs >= 85][pident >= 85][length >= 1000]
part_blast[sstart > send, c("sstart", "send") := .(send, sstart)]
part_blast[, strain := filename]
part_blast[, pls_id := sseqid]
print(paste0("Blast result contains:", nrow(part_blast)))

#
gc()

#
part_blast <- merge.data.table(part_blast, pls_clu[, .(sseqid, cluster)], by = "sseqid")

#
part_blast[, filename := NULL]
part_blast[, sseqid := NULL]
part_blast[, evalue := NULL]
part_blast[, bitscore := NULL]
part_blast[, qcovhsp := NULL]

#
gc()



##### __Sec 3.2 Loop #####
##
target_out <- data.table(strain = NULL, pls_id = NULL, cluster = NULL, scovs = NULL, slen = NULL)
sub_out <- data.table(strain = NULL, pls_id = NULL, cluster = NULL, scovs = NULL, slen = NULL, pp = NULL, pc = NULL)


##
#
if (nrow(part_blast) == 0) {
  print("Empty blast result")
} else {
  for (strain_in in unique(part_blast$strain)) {
    ###
    ##
    print(paste0("Call plasmids for: ", strain_in))
    
    ## First ordered mapping
    #
    part_map <- part_blast[strain == strain_in, red_len(sstart, send), by = .(strain, pls_id, cluster, slen)][, .(strain, pls_id, cluster, slen, scovs = V1/slen)][scovs >= 0.8][order(strain, -slen)]
    
    #
    part_map <- part_map[part_map[ , .I[which.max(slen * scovs)], by = .(strain, cluster)]$V1][order(strain, -slen)]
    
    #
    part_out_in <- merge(part_blast, part_map[, .(strain, pls_id)], by = c("strain", "pls_id"))
    
    
    
    ###
    ##
    if (nrow(part_map) == 0) {
      print("No plasmid detected now")
    } else {
      while (nrow(part_map) > 0) {
        ## load the iteration hit
        #
        filename_in <- part_map[1, strain]
        cluster_in <- part_map[1, cluster]
        pls_id_in <- part_map[1, pls_id]
        print(paste0("Check: ", filename_in, " ", cluster_in))
        
        
        ## check the mge
        #
        mge_ctg <- mge_ngs[filename == filename_in][(qstart <= 100) | (qend >= (qlen - 100)), qseqid]
        mge_ctg <- part_out_in[qseqid %in% mge_ctg, .(N = uniqueN(pls_id)), .(strain, qseqid)][N > 1, qseqid]
        part_out_in <- part_out_in[!qseqid %in% mge_ctg]
        
        
        ## include the hit
        #
        target_per <- part_out_in[strain == filename_in][pls_id == pls_id_in][, red_len(sstart, send), by = .(strain, pls_id, cluster, slen)][, .(strain, pls_id, cluster, scovs = V1/slen, slen)]
        
        #
        if (target_per$scovs >= 0.8) {
          #
          target_out <- rbind(target_out, target_per)
          
          #
          print(paste0("Include: ", filename_in, " ", pls_id_in))
        } else {
          print(paste0("Removed (mge): ", filename_in, " ", pls_id_in))
        }
        
        
        ## check sub 
        #
        sub_part_in <- part_out_in[(qseqid %in% part_out_in[strain == filename_in][pls_id == pls_id_in, qseqid])][strain == filename_in][pls_id != pls_id_in]
        sub_part_out <- sub_part_in[, red_len(sstart, send), by = .(strain, pls_id, cluster, slen)][, .(strain, pls_id, cluster, scovs = V1/slen, slen)][scovs >= 0.9]
        
        #
        if (nrow(sub_part_out) != 0) {
          #
          part_out_blast <- part_out_in[strain == filename_in][pls_id == pls_id_in]
          part_out_blast_um <- unique(part_out_blast[, .(strain, pls_id, qseqid, qcovs, qlen)])[, .(strain, pls_id, qseqid, um_p = (1-qcovs/100)*qlen)]
          
          #
          sub_part_out_blast <- merge(sub_part_out[, .(strain, pls_id, scovs)], sub_part_in, by = c("strain", "pls_id"))
          sub_part_out_blast_um <- unique(sub_part_out_blast[, .(strain, pls_id, qseqid, qcovs, qlen)])[, .(strain, pls_id, qseqid, um_s = (1-qcovs/100)*qlen)]
          
          #
          for (subpls in unique(sub_part_out_blast_um$pls_id)) {
            print(paste0("check sub pls:", subpls))
            sub_um <- merge(sub_part_out_blast_um[strain == filename_in][pls_id == subpls], part_out_blast_um[, .(qseqid, um_p)], by = c("qseqid"))
            if (sub_um[, sum(um_s)] > sub_um[, sum(um_p)]) {
              print(paste0("Removed sub pls: ", filename_in, " ", subpls))
            } else {
              print(paste0("Keep sub pls: ", filename_in, " ", subpls))
              sub_out <- rbind(sub_out, sub_part_out[strain == filename_in][pls_id == subpls][, .(strain, pls_id, cluster, scovs, slen, pp = pls_id_in, pc = cluster_in)])
            }
          }
        }
        
        #
        gc()
        
        
        ## Filter and re-assign
        #
        part_out_in <- part_out_in[!(qseqid %in% part_out_in[strain == filename_in][pls_id == pls_id_in, qseqid])]
        
        # 
        rm(part_map)
        part_map <- part_out_in[, red_len(sstart, send), by = .(strain, pls_id, cluster, slen)][, .(strain, pls_id, cluster, slen, scovs = V1/slen)][scovs >= 0.8][order(strain, -slen)]
        
        #
        part_map <- part_map[part_map[ , .I[which.max(slen * scovs)], by = .(strain, cluster)]$V1][order(strain, -slen)]
        
        #
        rm(part_out_in)
        part_out_in <- merge(part_blast, part_map[, .(strain, pls_id)], by = c("strain", "pls_id"))
        
        #
        gc()
        
      }
    }
  }
}


##### __Sec 3.3 output #####
## write
if (nrow(sub_out) != 0) {
  sub_out_dt <- sub_out[, .(strain, cluster = pc)][, .N, by = c("strain", "cluster")][N > 1]
  if (file.exists(paste0(out_dir, "/", strain_name, "_subpls.tsv"))) file.remove(paste0(out_dir, "/", strain_name, "_subpls.tsv"))
  fwrite(sub_out, paste0(out_dir, "/", strain_name, "_subpls.tsv"), sep = '\t', eol = '\n')
  gc()
} else {
  print(paste0("No sub plasmids"))
}

# remove sub
if (nrow(sub_out_dt) != 0) {
  target_out <- target_out[!sub_out_dt, on=.(strain, cluster)]
} else {
  print(paste0("No sub plasmid pair"))
}



# write HC
if (nrow(target_out) != 0) {
  #
  if (file.exists(paste0(out_dir, "/", strain_name, "_pls.tsv"))) file.remove(paste0(out_dir, "/", strain_name, "_pls.tsv"))
  fwrite(target_out, paste0(out_dir, "/", strain_name, "_pls.tsv"), sep = '\t', eol = '\n')
  #
  if (file.exists(paste0(out_dir, "/", strain_name, "_pls_raw.tsv"))) file.remove(paste0(out_dir, "/", strain_name, "_pls_raw.tsv"))
  fwrite(merge(target_out, part_blast, by = c("strain", "cluster", "pls_id", "slen")),
         paste0(out_dir, "/", strain_name, "_pls_raw.tsv"), sep = '\t', eol = '\n')
} else {
  print(paste0("No HC plasmids"))
}

#
gc()


##
#
print(paste0("END: ", strain_name, " in ", Sys.time()))



