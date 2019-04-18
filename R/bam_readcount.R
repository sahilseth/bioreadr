
# INSTALL -------

# module load conda_/3.6
# conda install -c bioconda bam-readcount 
# 



# Usage: bam-readcount [OPTIONS] <bam_file> [region]
# Generate metrics for bam_file at single nucleotide positions.
# Example: bam-readcount -f ref.fa some.bam
# Available options:
#   -h [ --help ]                         produce this message
# -v [ --version ]                      output the version
# -q [ --min-mapping-quality ] arg (=0) minimum mapping quality of reads used
# for counting.
# -b [ --min-base-quality ] arg (=0)    minimum base quality at a position to
# use the read for counting.
# -d [ --max-count ] arg (=10000000)    max depth to avoid excessive memory
# usage.
# -l [ --site-list ] arg                file containing a list of regions to
# report readcounts within.
# -f [ --reference-fasta ] arg          reference sequence in the fasta format.
# -D [ --print-individual-mapq ] arg    report the mapping qualities as a comma
# separated list.
# -p [ --per-library ]                  report results by library.
# -w [ --max-warnings ] arg             maximum number of warnings of each type
# to emit. -1 gives an unlimited number.
# -i [ --insertion-centric ]            generate indel centric readcounts.
# Reads containing insertions will not be
# included in per-base counts

# error/warnings -------
#' Question: How To Fix A Bam-Readcount Sm Error: "Couldn'T Grab Single-End Mapping Quality For Read"
#' This error message has to do with the SM tags on the read which is the single-ended mapping quality. 
#' Some aligners do not report this so bam-readcount reports this warning when it cannot grab the value. 
#' It does not invalidate your results and you can safely ignore this error, 
#' but you should not use the single ended mapping quality field of the output.



#' Run BAM readcount
#'
#' @param bam something
#' @param bed something
#' @param bamreadcount_exe something
#' @param samplename something
#' @param outfile something
#' @param fa_fl something
#'
#' @export
bam_readcount <- function(bam, 
                          samplename,
                          bed,
                          outfile,
                          fa_fl = "~/ref/human/b37/fastas/Homo_sapiens_assembly19.fasta",
                          bamreadcount_exe = "bam-readcount"){
  
  
  # bam-readcount [OPTIONS] <bam_file> [region]
  # -l [ --site-list ] arg                file containing a list of regions to
  # report readcounts within
  # The list of regions should be formatted as chromosome start and end. 
  # Each field should be tab separated and coordinates should be 1-based.
  cmd = glue("{bamreadcount_exe} -l {bed} -f {fa_fl} {bam} > {outfile}") %>% as.character()
  
  #system(cmd)
  flowmat = to_flowmat(list(bamreadcount = cmd), samplename = samplename)
  
  ret = list(flowmat = flowmat)
  return(ret)
}


#' parseline
#'
#' @param x something
#'
#' @export
bam_readcount.parseline <- function(x){
  cols = "base:count:avg_mapping_quality:avg_basequality:avg_se_mapping_quality:num_plus_strand:num_minus_strand:avg_pos_as_fraction:avg_num_mismatches_as_fraction:avg_sum_mismatch_qualities:num_q2_containing_reads:avg_distance_to_q2_start_in_q2_reads:avg_clipped_length:avg_distance_to_effective_3p_end"
  cols = strsplit(cols, ":")[[1]]
  
  # split line
  rw = strsplit(x, split = "\t")[[1]]
  
  # split each base counts
  tmp <- lapply(rw[-c(1:4)], function(bc){
    bc = strsplit(bc, split = ":")[[1]]
  }) %>% do.call(rbind, .)
  
  colnames(tmp) = cols
  
  dat = data.frame(chr = rw[1], pos = rw[2], 
                   ref = rw[3], dp = rw[4],
                   tmp, stringsAsFactors = F)
  
}



# chr	position	reference_base	depth	base:count:avg_mapping_quality:avg_basequality:avg_se_mapping_quality:num_plus_strand:num_minus_strand:avg_pos_as_fraction:avg_num_mismatches_as_fraction:avg_sum_mismatch_qualities:num_q2_containing_reads:avg_distance_to_q2_start_in_q2_reads:avg_clipped_length:avg_distance_to_effective_3p_end   ...
#x = '~/projects2/av_pdac_ct/analysis/exome/set1/snv/bamreadcount/GDraetts-CSMS-G001-PATC53R1-T_C8P7HACXX-1-AACGCTTA.bwa_recalibed_bamreadcnt.tsv'

#' parse read count
#'
#' @param x something
#' @param samplename something
#' @param bed something
#'
#' @export
bam_readcount.parse <- function(x, 
                                samplename,
                                bed){
  
  library(params)
  library(dplyr)
  
  tmp = scan(x, what = "character", sep = "\n") %>% unique()
  
  dat = lapply(tmp, bam_readcount.parseline) %>% bind_rows()
  dat2 = select(dat, chr, start = pos, ref_allele = ref, alt_allele = base, alt_count = count) %>% 
    mutate(start = as.integer(start), 
           end = start, 
           alt_count = as.integer(alt_count), 
           samplename = samplename) %>% 
    dplyr::group_by(chr, start) %>% 
    mutate(dp = sum(alt_count), 
           af = alt_count / dp) %>% 
    ungroup()
  
  # add ref counts
  dat_ref = dplyr::select(dat2, chr, start, ref_allele, alt_allele, alt_count) %>% 
    filter(ref_allele == alt_allele) %>% mutate(ref_count = alt_count) %>% 
    select(-alt_count, -ref_allele, -alt_allele)
  
  dat2 = left_join(dat2, dat_ref, by = c("chr", "start"))
  
  # combine with bed
  if(!is.data.frame(bed)){
    if(file.exists(bed)){
      bed = read_tsv(bed, col_types = cols(.default = col_character())) %>% 
        mutate(start = as.integer(start), 
               end = as.integer(end))
    }
  }
    
  
  dat3 = left_join(bed, dat2, by = c("chr", "start", "end", 
                                     "ref_allele", "alt_allele"))
  

}



#' bam_readcount_r
#' 
#' read add mutect files, filter and create bed. quick R func, for somatic variants
#'
#' @param trk can be a data.frame, or a tsv file with columns: sample1, mutect_fl, file1 (recalibed bam)
#' @param execute something
#' @param bamreadcount_exe something 
#' @param fa_fl something
#' 
#' @import flowr
#' @import parallel
#' 
#' @export
bam_readcount_r <- function(trk, 
                            col_fl = "MUT", 
                            col_bam = "BAM",
                            col_samp = "NAME",
                            execute = F,
                            bamreadcount_exe = "~/apps/conda/3.6/bin/bam-readcount",
                            force_redo = F,
                            fa_fl = "~/ref/human/b37/fastas/Homo_sapiens_assembly19.fasta"
                            
                            ){
  
  # pacman::p_load(flowr, parallel)
  
  if(!is.data.frame(trk))
    trk = read_tsv(trk)
  
  trk = data.frame(trk, check.names = F, stringsAsFactors = F)
  
  # # make sure trk has the reqd columns
  # cols_expected = c("mutect_fl", "sample1", "file1")
  # testthat::expect_named(trk, cols_expected, ignore.order = TRUE, ignore.case = TRUE)

  message("reading mutect ...")
  df_mutect = mutect.read(trk, 
                          col_fl = col_fl,
                          col_samp = col_samp)
  # create well annotated bed, judgement is always KEEP
  message("\ncreating uniq bed ...")
  df_mutect_bed = dplyr::select(df_mutect, chr, start, end, ref_allele, alt_allele, 
                                context, key:entrez_gene_id, 
                                # should add aaannotation
                                aaannotation) %>% unique()
  
  # write out bed, and run bam read count on each bam file:
  write_tsv(df_mutect_bed, "df_bed.tsv")
  
  # get bamreadcount fls:
  message("get bamreacount cmds ...")
  bamreadcnt_fl = paste0(trk[, "NAME"], ".bamreadcount.tsv")
  # run bam read count
  # source('~/Dropbox/public/flow-r/my.ultraseq/my.ultraseq/R/bam_readcount.R')
  bamreadcnt_fl = paste0(trk[, "NAME"], ".bamreadcount.tsv")
  bamreadcnt_ann_fl = paste0(trk[, "NAME"], ".bamreadcount_ann.tsv")
  # bam-readcount [OPTIONS] <bam_file> [region]
  # cmd = glue("{bamreadcount_exe} -l {bed} -f {fa_fl} {bam} > {outfile}") %>% 
  #   as.character()
  out = bam_readcount(bam = trk[, "BAM"], 
                      samplename = trk[, "NAME"], 
                      outfile = bamreadcnt_fl, 
                      fa_fl = fa_fl,
                      bamreadcount_exe = bamreadcount_exe,
                      bed = "df_bed.tsv")
  
  message("run bamreacount ...")
  out$flowmat$cmd = paste0(out$flowmat$cmd, " 2> ", trk[, "NAME"], ".bamreadcount.log");
  print(out$flowmat$cmd)
  
  # if execute and some files do not exists
  # OR if redo == T
  if(execute & any(!file.exists(bamreadcnt_fl)) | force_redo)
    tmp = parallel::mclapply(out$flowmat$cmd, system, mc.cores = length(out$flowmat$cmd))
  
  message("parse each of the files ...")
  tmp = lapply(1:nrow(trk), function(i){
    df_bamreadcount = bam_readcount.parse(x = bamreadcnt_fl[i], 
                        samplename = trk[i, "NAME"], 
                        bed = "df_bed.tsv")
    write_tsv(df_bamreadcount, bamreadcnt_ann_fl[i])
  }) %>% bind_rows()
  
  trk$MUT_RECALL = bamreadcnt_ann_fl
  write_tsv(trk, "trk.tsv")
  write_rds(tmp, "df_bamreadcount.rds")
  
  list(trk = trk, df_bamreadcount = tmp)
  
}















# END