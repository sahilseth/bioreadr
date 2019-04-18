


# this ONE function, runs the entire pipeline
# in the end, we will just call this function
# bam-readcount only needs 2-3 cores.

# TEST --------
if(FALSE){
  # trk = df_trk %>% filter(path_patient_id == "185_057")
  
  # module load conda_/3.6
  # conda install -c bioconda bam-readcount
  # Install MiniConda distribution and setup a Python 2.7 environment.
  # conda install -c aroth85 pyclone
  library(pacman)
  p_load(tidyverse, flowr)
  
  pipe_src = "~/Dropbox/public/flow-r/my.ultraseq/my.ultraseq/R/bam_readcount_pyclone.R"
  flow_def = gsub(".R$", ".def", pipe_src)
  flow_conf = "~/Dropbox/public/flow-r/my.ultraseq/pipelines/dnaseq/facets_pyclone/bam_facets.conf"
  source(pipe_src);opts_flow$load(flow_conf)
  
  
  setwd("/rsrch2/iacs/iacs_dep/sseth/flows/SS/tnbc/ms51_wex_b1/pyclone/185_057")
  # setwd("/rsrch2/iacs/iacs_dep/sseth/flows/SS/tnbc/ms51_wex_b1/pyclone/185_198")
  
  # bam read count
  bam_readcount_pyclone_r_step1(trk = "trk.tsv")
  # filter bam read count, pyclone
  bam_readcount_pyclone_r_step2()
  # plot pyclone, citup
  bam_readcount_pyclone_r_step3()
  
  trk = readr::read_tsv("trk.tsv")
  # copy over files to this folder
  # this would be done, using a transfer script next time!
  file.copy(file.path("~/projects2/ss_tnbc/data/artemis/wex/2019_b1/mutect", trk$mutect_fl), ".")
  file.copy(file.path(trk$file1), ".")
  file.copy(paste0(trk$file1, ".bai"), ".")
  segs = trk$facets_fl_full %>% gsub("/home/sseth/projects2", "~/projects2", .)
  fits = gsub("_cncf.tsv", ".rds", segs)
  file.copy(segs, ".")
  file.copy(fits, ".")
  
  
}

read_clone_arch_input_trk <- function(trk, 
                                      normalize_file_names = F,
                                      check_files_exist = F){
  
  pacman::p_load(testit)
  
  if(!is.data.frame(trk))
    trk = readr::read_tsv(trk) %>% data.frame(stringsAsFactors = F)
  
  # bam file should already have been transferred
  if(normalize_file_names)
    trk %<>% 
      mutate(BAM = basename(BAM),
             MUT = basename(MUT),
             CNV = basename(CNV))
  
  # TESTS
  # make sure trk has the reqd columns
  cols_expected = c("MUT", "NAME", "BAM", "CNV")
  # testthat::expect_(cols_expected %in% colnames(trk))
  testit::assert("we have reqd columns in trk", {
    cols_expected %in% colnames(trk)
  })
  
  # make sure files exist
  bams = trk$BAM
  bais = gsub(".bam$", ".bam.bai", bams)
  muts = trk$MUT
  segs = trk$CNV
  # this will change based on input!
  fits = gsub("_cncf.tsv", ".rds", segs)
  
  # ONLY REMOVE TEMPORARILY!!
  if(check_files_exist)
    testit::assert("check if all files exists", {
      # file.exists(bams) &
      #   file.exists(bais) &
        file.exists(segs) &
        file.exists(fits) &
        file.exists(muts)
    })
  
  trk
}



#' this function actually runs the pipeline
#' 
#' we can possible, call this using the cmd line using ultraseq
#'
#' @param trk trk fl
#' @param bamreadcount_exe bamreadcount_exe
#' @param fa_fl fa_fl 
#' @param pyclone_params pyclone_params
#' @param segfl_type segfl_type 
#' @param mutfl_type mutfl_type 
#' @param conf conf
#' 
#'
#' @export
#'
pyclone_pipe_v1_step1_r <- function(trk,
                                    bamreadcount_exe = opts_flow$get("bamreadcount_exe"),
                                    fa_fl = opts_flow$get("ref_fasta"),
                                    
                                    pyclone_params = opts_flow$get("pyclone_params"),
                                    segfl_type = opts_flow$get("segfl_type"), 
                                    mutfl_type = opts_flow$get("mutfl_type"),

                                    conf = "~/Dropbox/public/flow-r/my.ultraseq/pipelines/dnaseq/facets_pyclone/facets_pyclone_citup_v1.conf"){
  
  opts_flow$load(conf)
  
  pacman::p_load(tidyverse, glue, janitor, magrittr, testit)

  trk = read_clone_arch_input_trk(trk, normalize_file_names = T)
  
  # cmd_bam_readcount <- glue("funr my.ultraseq::bam_readcount_r")
  
  # count them!
  # source('~/Dropbox/public/flow-r/my.ultraseq/my.ultraseq/R/pyclone.R')
  # source('~/Dropbox/public/flow-r/my.ultraseq/my.ultraseq/R/gm_mutect.R')
  # source('~/Dropbox/public/flow-r/my.ultraseq/my.ultraseq/R/bam_readcount.R')
  out = bam_readcount_r(trk,
                        col_fl = "MUT",
                        col_samp = "NAME",
                        bamreadcount_exe = bamreadcount_exe,
                        fa_fl = fa_fl,
                        execute = T)
  

}





#' pyclone_pipe_v2_step1_r
#' 
#' filter bamreadcount, and run pyclone
#'
#' @param pyclone_path tmp, path with all relavent files
#' @param trk trk 
#' @param conf conf
#' 
#' @import readr
#' @import dplyr
#' @import janitor
#' @import flowr
#'
#' @export
pyclone_pipe_v1_step2_r <- function(pyclone_path = ".", 
                                    trk = "trk.tsv",
                                    conf = "~/Dropbox/public/flow-r/my.ultraseq/pipelines/dnaseq/facets_pyclone/facets_pyclone_citup_v1.conf", 
                                    force_redo = FALSE,
                                    clean_bams = FALSE){
  
  pacman::p_load(tidyverse, glue, janitor, magrittr)
  flowr::opts_flow$load(conf)
  
  trk = read_clone_arch_input_trk(trk, normalize_file_names = T)
  
  # out = bam_readcount_r(trk, execute = execute)
  # df_mut_recall = mutect.read(trk, 
  #                             col_samp = "NAME",
  #                             col_fl = "MUT_RECALL")
  df_bamreadcount = read_rds("df_bamreadcount.rds")
  
  message("we have a total of ", n_distinct(df_bamreadcount$key), " mutations")
  # add a filtering step
  # get min DP, min/max AF, per mutation, across samples 
  df_bamreadcount = group_by(df_bamreadcount, key) %>% 
    mutate(min_dp = min(dp), 
           min_af = min(af),
           max_af = max(af),
           max_alt_count = max(alt_count)) %>% ungroup()
  
  # mut which represent at least 5% (in one of the samples)
  # at least 2 supporting reads!
  # minimum read depth (for good calls)
  df_bamreadcount = filter(df_bamreadcount, 
                           min_dp > 20, 
                           max_af > 0.05,
                           max_alt_count > 2) # assuming, alt count is max, where AF is max...
  message("filtering # mutations: ", n_distinct(df_bamreadcount$key))
  write_rds(df_bamreadcount, "df_bamreadcount_f1.rds")
  
  bamreadcount_f_fl = lapply(unique(df_bamreadcount$samplename), function(samp){
    fl = paste0(samp, ".bamreadcount_ann_f.tsv")
    filter(df_bamreadcount, samplename == samp) %>% write_tsv(fl)
    fl
  }) %>% unlist()
  trk$MUT_RECALL_F1 = bamreadcount_f_fl
  # update trk, with readcount file
  write_tsv(trk, "trk.tsv")
  
  # filter(df_bamreadcount, gene_knowngene == "EGFR")
  # filter(df_bamreadcount, gene_knowngene == "SETD2")
  
  message("bam readcount is all done, can clean up bams")
  if(clean_bams){
    system("rm *bam;rm *bai")
  }
  #trk = out$trk
  # pyclone inputs
  #source('~/Dropbox/public/flow-r/my.ultraseq/my.ultraseq/R/pyclone.R')
  # run pyclone ------
  # parse and annotate clusters
  # re-do basic plots
  # debug(pyclone_pipe_r)
  pyclone_pipe_r(trk, 
                 pyclone_params = opts_flow$get("pyclone_params"),
                 segfl_type = opts_flow$get("segfl_type"), 
                 mutfl_type = opts_flow$get("mutfl_type"), force_redo = force_redo)
  
  # its fine we can plot twice :)
  pyclone_plots(pyclone_path)
  
  
}


#' bam_readcount_pyclone_r_step3
#' 
#' pyclone plots citup/cloneevol
#'
#' @param pyclone_path pyclone_path
#' @param conf conf 
#' 
#' @import dplyr
#' 
#' @export
pyclone_pipe_v1_step3_r <- function(pyclone_path = ".", 
                                    conf = "~/Dropbox/public/flow-r/my.ultraseq/pipelines/dnaseq/facets_pyclone/facets_pyclone_citup_v1.conf"){
  setwd(pyclone_path)
  
  opts_flow$load(conf)
  
  # read in the input file
  
  # pyclone plots
  source('~/Dropbox/public/flow-r/my.ultraseq/my.ultraseq/R/pyclone.R')
  source('~/Dropbox/public/flow-r/my.ultraseq/my.ultraseq/R/citup.R')
  
  pyclone_plots(pyclone_path)
  
  # pyclone_parse_cloneevol
  
  citup_exe = "module load conda_/2.7;source activate citup;run_citup_qip.py" 
  citup_exe = "module load citup/0.1.0;/risapps/rhel6/citup/anaconda2-4.2/bin/run_citup_qip.py" 
  citup_params = "--submit local --loglevel DEBUG --maxjobs 8"
  # pyclone citup
  pyclone_citup(pyclone_path,
                citup_exe = citup_exe,
                citup_params = citup_params)
  
  # timescape
  
  
}



