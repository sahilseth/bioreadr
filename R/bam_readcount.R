


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
#' @param bam 
#' @param bed 
#' @param bamreadcount_exe 
#'
#' @return
#' @export
#'
#' @examples
bam_readcount <- function(bam, 
                          samplename,
                          bed,
                          outfile,
                          fa_fl = "~/ref/human/b37/fastas/Homo_sapiens_assembly19.fasta",
                          bamreadcount_exe = "bam-readcount"){
  
  
  # bam-readcount [OPTIONS] <bam_file> [region]
  cmd = glue("{bamreadcount_exe} -l {bed} -f {fa_fl} {bam} > {outfile}")
  
  #system(cmd)
  flowmat = to_flowmat(list(bamreadcount = cmd), samplename = samplename)
  
  ret = list(flowmat = flowmat)
  return(ret)
  
  
  
}


#' Title
#'
#' @param x 
#'
#' @return
#' @export
#'
#' @examples
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
  if(!is.data.frame(bed))
    if(file.exists(bed))
      bed = read_tsv(bed)
  
  dat3 = left_join(bed, dat2, by = c("chr", "start", "end", 
                                     "ref_allele", "alt_allele"))
  

}



#' bam_readcount_r
#' 
#' # read add mutect files, filter and create bed
#' 
#' quick R func, for somatic variants
#'
#' @return
#' @export
#'
#' @examples
bam_readcount_r <- function(trk, execute = F){
  
  p_load(flowr, parallel)
  
  message("reading mutect ...")
  df_mutect = mutect.read(trk)
  # create well annotated bed, judgement is always KEEP
  message("\ncreating uniq bed ...")
  df_mutect_bed = dplyr::select(df_mutect, chr, start, end, ref_allele, alt_allele, context, key:entrez_gene_id) %>% unique()
  
  # write out bed, and run bam read count on each bam file:
  write_tsv(df_mutect_bed, "df_bed.tsv")
  
  # get bamreadcount fls:
  message("get bamreacount cmds ...")
  bamreadcnt_fl = paste0(trk$samplelbl1, ".bamreadcount.tsv")
  # run bam read count
  source('~/Dropbox/public/flow-r/my.ultraseq/my.ultraseq/R/bam_readcount.R')
  bamreadcnt_fl = paste0(trk$samplelbl1, ".bamreadcount.tsv")
  bamreadcnt_ann_fl = paste0(trk$samplelbl1, ".bamreadcount_ann.tsv")
  out = bam_readcount(file.path("/rsrch1/iacs/ngs_runs", trk$runid, "bams", trk$bam1), 
                          samplename = trk$samplelbl1, 
                          outfile = bamreadcnt_fl, 
                          bed = "df_bed.tsv")
  
  message("run bamreacount ...")
  out$flowmat$cmd = paste0(out$flowmat$cmd, " 2> ", trk$samplelbl1, ".bamreadcount.log");
  print(out$flowmat$cmd)
  
  if(execute)
    tmp = mclapply(out$flowmat$cmd, system, mc.cores = 5)
  
  # parse each of the files
  tmp = lapply(1:nrow(trk), function(i){
    df_bamreadcount = bam_readcount.parse(x = bamreadcnt_fl[i], 
                        samplename = trk$samplelbl1[i], 
                        bed = "df_bed.tsv")
    write_tsv(df_bamreadcount, bamreadcnt_ann_fl[i])
  }) %>% bind_rows()
  
  trk$bamreadcnt_ann_fl = bamreadcnt_ann_fl
  
  list(trk = trk, df_bamreadcount = tmp)
  
}















# END