
# INSTALL -------

# module load conda_/3
# conda create --name bamreadcount
# conda activate bamreadcount
# conda install -c bioconda bam-readcount
# exe: ~/apps/conda/3/envs/bamreadcount/bin/bam-readcount



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
#' https://github.com/genome/bam-readcount/issues/25
#' Some aligners do not report this so bam-readcount reports this warning when it cannot grab the value.
#' It does not invalidate your results and you can safely ignore this error,
#' but you should not use the single ended mapping quality field of the output.


# main ------

# ** test ------

if (FALSE) {
  # test the function
  # cd /rsrch3/scratch/iacs/sseth/flows/SS/tnbc/ms51_wex/ss_cl_het_gb/runs/185_012/tmp
  # R CMD INSTALL --no-multiarch --with-keep.source ~/Dropbox/public/flowr/my.ultraseq/my.ultraseq
  # funr my.ultraseq::bam_readcount_r.ssm individual=185_012,185_012,185_012 name=IPCT-S4012-MOON0051-Cap2509-8-HTID301,IPCT-S4013-MOON0051-Cap2527-4-HTID229,IPCT-S4014-MOON0051-Cap2544-4-HTID003 normal=YES,NO,NO bam=/rsrch3/scratch/iacs/sseth/flows/SS/tnbc/ms51_wex/ss_cl_het_gb/bam/IPCT-S4012-MOON0051-Cap2509-8-HTID301_ngs-pipe-2-GATACTGCACGG.bwa_recalibed.bam,/rsrch3/scratch/iacs/sseth/flows/SS/tnbc/ms51_wex/ss_cl_het_gb/bam/IPCT-S4013-MOON0051-Cap2527-4-HTID229_ngs-pipe-1-AGGTCGAA.bwa_recalibed.bam,/rsrch3/scratch/iacs/sseth/flows/SS/tnbc/ms51_wex/ss_cl_het_gb/bam/IPCT-S4014-MOON0051-Cap2544-4-HTID003_ngs-pipe-1-TAGATGCCGTCC.bwa_recalibed.bam timepoint=GB,T0,T1 gender=f,f,f outprefix_paired=IPCT-S4012-MOON0051-Cap2509-8-HTID301___matched,IPCT-S4013-MOON0051-Cap2527-4-HTID229___matched,IPCT-S4014-MOON0051-Cap2544-4-HTID003___matched conf=~/Dropbox/public/flowr/my.ultraseq/projects/ss_tnbc/ms51/wex/flowr_conf.toml execute=T
  # # try 001
  # cd /rsrch3/scratch/iacs/sseth/flows/SS/tnbc/ms51_wex/ss_cl_het_noref/runs/185_001/tmp
  # R CMD INSTALL --no-multiarch --with-keep.source ~/Dropbox/public/flowr/my.ultraseq/my.ultraseq
  # funr my.ultraseq::bam_readcount_r.ssm bamreadcntdir=ssm/bamreadcnt individual=185_001,185_001 name=IPCT-S4013-MOON0051-Cap2527-4-HTID193,IPCT-S4013-MOON0051-Cap2527-4-HTID205 normal=NO,NO bam=IPCT-S4013-MOON0051-Cap2527-4-HTID193_ngs-pipe-1-CTGATATG.bwa_recalibed.bam,IPCT-S4013-MOON0051-Cap2527-4-HTID205_ngs-pipe-1-ATCCAGTA.bwa_recalibed.bam timepoint=T0,T1 gender=f,f outprefix=IPCT-S4013-MOON0051-Cap2527-4-HTID193,IPCT-S4013-MOON0051-Cap2527-4-HTID205 outprefix_paired=IPCT-S4013-MOON0051-Cap2527-4-HTID193___pon,IPCT-S4013-MOON0051-Cap2527-4-HTID205___pon conf=~/Dropbox/public/flowr/my.ultraseq/projects/ss_tnbc/ms51/wex/flowr_conf.toml execute=T
  # R CMD INSTALL --no-multiarch --with-keep.source ~/Dropbox/public/flowr/my.ultraseq/my.ultraseq
  # bash ../009.bamreadcnt/bamreadcnt_cmd_1.sh
}

#' bam_readcount_ssm
#'
#' @param trk
#' @param samplename
#' @param bamreadcntdir
#' @param variant_calling_mode
#' @param flow_conf
#'
#' @export
#'
#' @examples
bam_readcount_ssm <- function(trk, samplename, bamreadcntdir,
                              variant_calling_mode,
                              flow_conf = "~/Dropbox/public/flowr/my.ultraseq/projects/ss_tnbc/ms51/wex/flowr_conf.toml") {
  # using trk
  # write out the trk in bamreadcnt folder:
  # this already includes all inputs using prep_filenms

  # trk_fl = glue("{bamreadcountdir}/trk.tsv")
  # write_tsv(trk, trk_fl)
  # get cmd#' @export
  # supply a list of files
  individual <- paste0(trk$individual, collapse = ",")
  name <- paste0(trk$name, collapse = ",")
  normal <- paste0(trk$normal, collapse = ",")
  outprefix <- paste0(trk$outprefix, collapse = ",")
  bam <- paste0(trk$bam, collapse = ",")
  timepoint <- paste0(trk$timepoint, collapse = ",")
  gender <- paste0(trk$gender, collapse = ",")
  outprefix_paired <- paste0(trk$outprefix_paired, collapse = ",")
  cmd_bamreadcnt <- glue("funr my.ultraseq::bam_readcount_r.ssm bamreadcntdir={bamreadcntdir} individual={individual} variant_calling_mode={variant_calling_mode} name={name} normal={normal} bam={bam} timepoint={timepoint} gender={gender} outprefix={outprefix} outprefix_paired={outprefix_paired} conf={flow_conf} execute=TRUE force_redo=FALSE")
  cmds <- list(bamreadcnt = cmd_bamreadcnt)

  flowmat <- to_flowmat(cmds, samplename = samplename)

  list(flowmat = flowmat, outfiles = list(), trk = trk)
}

# funcs ------
bamrdcnt_prep_filenms <- function(trk, bamreadcntdir, variant_calling_mode) {
  params::check_args()

  trk %<>% metadata_for_dnaseq(trk)
  wranglr::expect_columns(trk, c("outprefix", "outprefix_paired", "gender"))

  # add dir to outprefix
  trk <- dplyr::mutate(trk,
    bamreadcnt_prefix_paired = file.path(bamreadcntdir, outprefix_paired)
  )

  # ceate a new trk with ALL reqd files
  trk <- trk %>%
    dplyr::mutate(
      # add bam readcount files
      # combined pt level bed file
      bed_fl = glue("{bamreadcntdir}/df_bed.tsv"),
      bamreadcnt_fl = glue("{bamreadcnt_prefix_paired}.bamreadcount.tsv"),
      bamreadcnt_log = glue("{bamreadcnt_prefix_paired}.bamreadcount.log"),
      bamreadcnt_ann_fl = glue("{bamreadcnt_prefix_paired}.bamreadcount_ann.tsv"),
      # combined parsed bam readcount file
      bamreadcnt_rds = glue("{bamreadcntdir}/df_bamreadcount.rds"),
      # filter bamreadcount:
      bamreadcnt_rds_f1 = glue("{bamreadcntdir}/df_bamreadcount_f1.rds"),
      bamreadcnt_ann_fl_f1 = glue("{bamreadcnt_prefix_paired}.bamreadcount_ann_f1.tsv")
    )
  # add SSM mut files (inputs)
  if (variant_calling_mode == "matched") {
    trk %<>% mutate(pcgr_combvcf_rds = glue("ssm/m1_m2_muse/pcgr/{outprefix_paired}_combvcf_snp.pcgr_acmg.grch37_combvcf.rds"))
  } else {
    trk %<>% mutate(pcgr_combvcf_rds = glue("ssm/m1_m2/pcgr/{outprefix_paired}_combvcf_snp.pcgr_acmg.grch37_combvcf.rds"))
  }

  trk
}

#' Run BAM readcount
#'
#' @param bam something
#' @param bed something
#' @param bamreadcount_exe something
#' @param samplename something
#' @param outfile something
#' @param fa_fl something
#'
#' @details should be 1-based:
#'
#' The list of regions should be formatted as chromosome start and end. Each field should be tab separated and coordinates should be 1-based.
#'
#' @export
bam_readcount <- function(bam,
                          samplename,
                          bed,
                          outfile,
                          fa_fl = "~/ref/human/b37/fastas/Homo_sapiens_assembly19.fasta",
                          bamreadcount_exe = "~/apps/conda/3/envs/bamreadcount/bin/bam-readcount") {


  # bam-readcount [OPTIONS] <bam_file> [region]
  # -l [ --site-list ] arg                file containing a list of regions to
  # report readcounts within
  # The list of regions should be formatted as chromosome start and end.
  # Each field should be tab separated and coordinates should be 1-based.
  cmd <- glue("{bamreadcount_exe} -l {bed} -f {fa_fl} {bam} > {outfile}") %>% as.character()

  # system(cmd)
  flowmat <- to_flowmat(list(bamreadcount = cmd), samplename = samplename)

  ret <- list(flowmat = flowmat)
  return(ret)
}


#' parseline
#'
#' @param x something
#'
#' @export
bam_readcount.parseline <- function(x) {
  cols <- "base:count:avg_mapping_quality:avg_basequality:avg_se_mapping_quality:num_plus_strand:num_minus_strand:avg_pos_as_fraction:avg_num_mismatches_as_fraction:avg_sum_mismatch_qualities:num_q2_containing_reads:avg_distance_to_q2_start_in_q2_reads:avg_clipped_length:avg_distance_to_effective_3p_end"
  cols <- strsplit(cols, ":")[[1]]

  # split line
  rw <- strsplit(x, split = "\t")[[1]]

  # split each base counts
  tmp <- lapply(rw[-c(1:4)], function(bc) {
    bc <- strsplit(bc, split = ":")[[1]]
  }) %>% do.call(rbind, .)

  colnames(tmp) <- cols

  dat <- data.frame(
    chr = rw[1], pos = rw[2],
    ref = rw[3], dp = rw[4],
    tmp, stringsAsFactors = F
  )
}



# chr	position	reference_base	depth	base:count:avg_mapping_quality:avg_basequality:avg_se_mapping_quality:num_plus_strand:num_minus_strand:avg_pos_as_fraction:avg_num_mismatches_as_fraction:avg_sum_mismatch_qualities:num_q2_containing_reads:avg_distance_to_q2_start_in_q2_reads:avg_clipped_length:avg_distance_to_effective_3p_end   ...
# x = '~/projects2/av_pdac_ct/analysis/exome/set1/snv/bamreadcount/GDraetts-CSMS-G001-PATC53R1-T_C8P7HACXX-1-AACGCTTA.bwa_recalibed_bamreadcnt.tsv'

#' parse read count
#'
#' @param x something
#' @param samplename something
#' @param bed something
#'
#' @export
bam_readcount.parse <- function(x,
                                samplename,
                                bed,
                                type = c("snv", "indel")) {
  library(params)
  library(dplyr)

  tmp <- scan(x, what = "character", sep = "\n") %>% unique()

  dat <- lapply(tmp, bam_readcount.parseline) %>% bind_rows()
  if (nrow(dat) == 0) {
    return(data.frame(key = NA, alt_count = NA, dp = NA))
  }

  dat2 <- dplyr::select(dat, chr, start = pos, ref_allele = ref, alt_allele = base, alt_count = count) %>%
    mutate(
      start = as.integer(start),
      end = start,
      alt_count = as.integer(alt_count),
      samplename = samplename
    ) %>%
    dplyr::group_by(chr, start) %>%
    mutate(
      dp = sum(alt_count),
      af = alt_count / dp
    ) %>%
    ungroup()

  # add ref counts
  dat_ref <- dplyr::select(dat2, chr, start, ref_allele, alt_allele, alt_count) %>%
    filter(ref_allele == alt_allele) %>%
    mutate(ref_count = alt_count) %>%
    select(-alt_count, -ref_allele, -alt_allele)

  dat2 <- left_join(dat2, dat_ref, by = c("chr", "start"))

  # combine with bed
  if (!is.data.frame(bed)) {
    if (file.exists(bed)) {
      bed <- read_tsv(bed, col_types = cols(.default = col_character()))
      colnames(bed)[1:5] <- c("chr", "start", "end", "ref_allele", "alt_allele")
      bed %<>% mutate(
        start = as.integer(start),
        end = as.integer(end)
      )
    }
  }

  dat3 <- tidylog::left_join(bed, dat2,
    by = c("chr", "start", "end", "ref_allele", "alt_allele")
  )
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
bam_readcount_r.mutect <- function(trk,
                                   col_fl = "MUT",
                                   col_bam = "BAM",
                                   col_samp = "NAME",
                                   execute = F,
                                   bamreadcount_exe = "~/apps/conda/3.6/bin/bam-readcount",
                                   force_redo = F,
                                   fa_fl = "~/ref/human/b37/fastas/Homo_sapiens_assembly19.fasta") {

  # pacman::p_load(flowr, parallel)

  if (!is.data.frame(trk)) {
    trk <- read_tsv(trk)
  }

  trk <- data.frame(trk, check.names = F, stringsAsFactors = F)

  # # make sure trk has the reqd columns
  # cols_expected = c("mutect_fl", "sample1", "file1")
  # testthat::expect_named(trk, cols_expected, ignore.order = TRUE, ignore.case = TRUE)

  message("reading mutect ...")
  df_mutect <- mutect.read(trk,
    col_fl = col_fl,
    col_samp = col_samp
  )
  # create well annotated bed, judgement is always KEEP
  message("\ncreating uniq bed ...")
  df_mutect_bed <- dplyr::select(
    df_mutect, chr, start, end, ref_allele, alt_allele,
    context, key:entrez_gene_id,
    # should add aaannotation
    aaannotation
  ) %>% unique()

  # write out bed, and run bam read count on each bam file:
  write_tsv(df_mutect_bed, "df_bed.tsv")

  # get bamreadcount fls:
  message("get bamreacount cmds ...")
  bamreadcnt_fl <- paste0(trk[, "NAME"], ".bamreadcount.tsv")
  # run bam read count
  # source('~/Dropbox/public/flow-r/my.ultraseq/my.ultraseq/R/bam_readcount.R')
  bamreadcnt_fl <- paste0(trk[, "NAME"], ".bamreadcount.tsv")
  bamreadcnt_ann_fl <- paste0(trk[, "NAME"], ".bamreadcount_ann.tsv")
  # bam-readcount [OPTIONS] <bam_file> [region]
  # cmd = glue("{bamreadcount_exe} -l {bed} -f {fa_fl} {bam} > {outfile}") %>%
  #   as.character()
  out <- bam_readcount(
    bam = trk[, "BAM"],
    samplename = trk[, "NAME"],
    outfile = bamreadcnt_fl,
    fa_fl = fa_fl,
    bamreadcount_exe = bamreadcount_exe,
    bed = "df_bed.tsv"
  )

  message("run bamreacount ...")
  out$flowmat$cmd <- paste0(out$flowmat$cmd, " 2> ", trk[, "NAME"], ".bamreadcount.log")
  print(out$flowmat$cmd)

  # if execute and some files do not exists
  # OR if redo == T
  if (execute & any(!file.exists(bamreadcnt_fl)) | force_redo) {
    tmp <- parallel::mclapply(out$flowmat$cmd, system, mc.cores = length(out$flowmat$cmd))
  }

  message("parse each of the files ...")
  tmp <- lapply(1:nrow(trk), function(i) {
    df_bamreadcount <- bam_readcount.parse(
      x = bamreadcnt_fl[i],
      samplename = trk[i, "NAME"],
      bed = "df_bed.tsv"
    )
    write_tsv(df_bamreadcount, bamreadcnt_ann_fl[i])
  }) %>% bind_rows()

  trk$MUT_RECALL <- bamreadcnt_ann_fl
  write_rds(tmp, "df_bamreadcount.rds")

  list(trk = trk, df_bamreadcount = tmp)
}

#' bam_readcount_r.ssm
#'
#' @param trk
#' @param conf
#' @param execute
#' @param bamreadcount_exe
#' @param force_redo
#' @param fa_fl
#'
#' @export
bam_readcount_r.ssm <- function(trk, conf,
                                bamreadcntdir,
                                individual,
                                name,
                                variant_calling_mode,
                                normal,
                                bam,
                                timepoint,
                                gender,
                                outprefix,
                                outprefix_paired,
                                execute = F,
                                col_fl = "pcgr_combvcf_rds",
                                bamreadcount_exe = "~/apps/conda/3/envs/bamreadcount/bin/bam-readcount",
                                force_redo = F,
                                fa_fl = "~/ref/human/b37/fastas/Homo_sapiens_assembly19.fasta") {
  source("~/Dropbox/public/flowr/my.ultraseq/my.ultraseq/R/bam_readcount.R")
  p_load(futile.logger)
  flog.info("create trk file")
  trk <- data.frame(
    individual = individual,
    name = name,
    normal = normal,
    bam = bam,
    timepoint = timepoint,
    gender = gender,
    outprefix = outprefix,
    outprefix_paired = outprefix_paired
  ) %>%
    bamrdcnt_prep_filenms(
      bamreadcntdir = bamreadcntdir,
      variant_calling_mode = variant_calling_mode
    )

  # pacman::p_load(flowr, parallel)
  trk %<>% metadata_for_dnaseq(trk) %>% data.frame(stringsAsFactors = FALSE)

  opts_flow$load_toml(conf)
  pacman::p_load(tidyverse, glue, janitor, magrittr, testit)

  # # make sure trk has the reqd columns
  # cols_expected = c("mutect_fl", "sample1", "file1")
  # testthat::expect_named(trk, cols_expected, ignore.order = TRUE, ignore.case = TRUE)

  flog.info("reading ssm ...")
  source("~/Dropbox/public/flowr/my.ultraseq/my.ultraseq/R/ssm_het.R")
  df_ssm <- het_ssm.read(trk,
    col_fl = col_fl,
    col_samp = "name"
  )
  # create well annotated bed, judgement is always KEEP
  flog.info("\ncreating uniq bed ...")
  # select only SSMs
  df_ssm_bed <- dplyr::mutate(df_ssm, end = pos) %>%
    dplyr::select(
      chr = chrom, start = pos, end,
      ref_allele = ref,
      alt_allele = alt,
      symbol, consequence,
      # values used later:
      key,
      # context, key:entrez_gene_id,
      # should add aaannotation
      protein_change
    ) %>%
    unique()
  # lets get bam read count for all, will parse SHIT later!
  # df_ssm_bed = tidylog::filter(df_ssm_bed, nchar(ref_allele) == 1, nchar(alt_allele) == 1)

  # write out bed, and run bam read count on each bam file:
  write_tsv(df_ssm_bed, trk$bed_fl[1])

  # get bamreadcount fls:
  flog.info("get bamreacount cmds ...")
  # run bam read count
  # source('~/Dropbox/public/flow-r/my.ultraseq/my.ultraseq/R/bam_readcount.R')
  # bam-readcount [OPTIONS] <bam_file> [region]
  # cmd = glue("{bamreadcount_exe} -l {bed} -f {fa_fl} {bam} > {outfile}") %>%
  #   as.character()
  out <- bam_readcount(
    bam = trk$bam,
    samplename = trk$name,
    outfile = trk$bamreadcnt_fl,
    fa_fl = fa_fl,
    bamreadcount_exe = bamreadcount_exe,
    bed = trk$bed_fl[1]
  )
  # single end mapping issue:
  # https://github.com/genome/bam-readcount/issues/25
  # This warning is thrown if the mapper doesn't populate the SM tag. Without the SM tag, the avg_se_mapping_quality field is unable to be calculated and shouldn't be considered. The rest of the program should be unaffected. I will add this point to the FAQ.
  message("run bamreacount ...")
  out$flowmat$cmd <- paste0(out$flowmat$cmd, " 2> ", unlist(trk$bamreadcnt_log))
  paste(out$flowmat$cmd, collapse = "\n") %>% cat()

  # if execute and some files do not exists
  # OR if redo == T
  print(trk$bamreadcnt_fl)
  if (execute & any(!file.exists(trk$bamreadcnt_fl)) | force_redo) {
    flog.info("running bamreadcnt")
    tmp <- parallel::mclapply(out$flowmat$cmd, system, mc.cores = length(out$flowmat$cmd))
  }

  flog.info("parse each of the files ...")
  # combine both files
  # undebug(bam_readcount.parse)
  tmp <- lapply(1:nrow(trk), function(i) {
    df_bamreadcount <- bam_readcount.parse(
      x = trk$bamreadcnt_fl[i],
      samplename = trk$name[i],
      bed = trk$bed_fl[1]
    )
    write_tsv(df_bamreadcount, trk$bamreadcnt_ann_fl[i])
  }) %>% bind_rows()

  write_rds(tmp, trk$bamreadcnt_rds[1])

  flog.info("filter bamreadcnt files")
  bam_readcount.filter(trk)

  list(trk = trk, df_bamreadcount = tmp)
}

bam_readcount.filter <- function(trk) {

  # trk = read_clone_arch_input_trk(trk, normalize_file_names = T)
  trk <- metadata_for_dnaseq(trk, normalize_file_names = T)

  # out = bam_readcount_r(trk, execute = execute)
  # df_mut_recall = mutect.read(trk,
  #                             col_samp = "NAME",
  #                             col_fl = "MUT_RECALL")
  flog.info("reading bamreadcount")
  df_bamreadcount <- read_rds(trk$bamreadcnt_rds[1])

  message("we have a total of ", n_distinct(df_bamreadcount$key), " mutations")
  # add a filtering step
  # get min DP, min/max AF, per mutation, across samples
  df_bamreadcount <- group_by(df_bamreadcount, key) %>%
    mutate(
      min_dp = min(dp),
      min_af = min(af),
      max_af = max(af),
      max_alt_count = max(alt_count)
    ) %>%
    ungroup()

  flog.info("filtering bamreadcount")
  # mut which represent at least 5% (in one of the samples)
  # at least 2 supporting reads!
  # minimum read depth (for good calls)
  df_bamreadcount <- filter(
    df_bamreadcount,
    min_dp > 20,
    max_af > 0.05,
    max_alt_count > 2
  ) # assuming, alt count is max, where AF is max...
  message("filtering # mutations: ", n_distinct(df_bamreadcount$key))
  write_rds(df_bamreadcount, trk$bamreadcnt_rds_f1[1])

  # write out files
  bamreadcount_f_fl <- lapply(unique(df_bamreadcount$samplename), function(samp) {
    outfl <- filter(trk, name == samp) %>%
      select(bamreadcnt_ann_fl_f1) %>%
      unlist()
    filter(df_bamreadcount, samplename == samp) %>% write_tsv(outfl)
  }) %>% unlist()
}

bam_readcount_pipe.df_mut <- function(df_mut) {


  # create bed


  out <- bam_readcount(
    bam = trk[, "BAM"],
    samplename = trk[, "NAME"],
    outfile = bamreadcnt_fl,
    fa_fl = fa_fl,
    bamreadcount_exe = bamreadcount_exe,
    bed = "df_bed.tsv"
  )

  message("run bamreacount ...")
  out$flowmat$cmd <- paste0(out$flowmat$cmd, " 2> ", trk[, "NAME"], ".bamreadcount.log")
  print(out$flowmat$cmd)

  # if execute and some files do not exists
  # OR if redo == T
  if (execute & any(!file.exists(bamreadcnt_fl)) | force_redo) {
    tmp <- parallel::mclapply(out$flowmat$cmd, system, mc.cores = length(out$flowmat$cmd))
  }

  message("parse each of the files ...")
  tmp <- lapply(1:nrow(trk), function(i) {
    df_bamreadcount <- bam_readcount.parse(
      x = bamreadcnt_fl[i],
      samplename = trk[i, "NAME"],
      bed = "df_bed.tsv"
    )
    write_tsv(df_bamreadcount, bamreadcnt_ann_fl[i])
  }) %>% bind_rows()

  trk$MUT_RECALL <- bamreadcnt_ann_fl
  write_tsv(trk, "trk.tsv")
  write_rds(tmp, "df_bamreadcount.rds")

  list(trk = trk, df_bamreadcount = tmp)
}










# END
