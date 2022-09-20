
pyclone_prep_filenms <- function(trk, pyclonedir = NULL, flowr_path = NULL) {
  check_args(ignore = c("flowr_path"))

  # check generic columns
  trk <- metadata_for_dnaseq(trk)
  expect_columns(trk, c("outprefix", "outprefix_paired", "gender"))

  # add dir to outprefix
  trk <- dplyr::mutate(trk,
    pyclone_prefix_paired = file.path(pyclonedir, outprefix_paired),
    pyclone_muts = file.path(pyclonedir, "pyclone_muts.tsv"),
    pyclone_h5_fl = file.path(pyclonedir, "pyclone_out.h5"),
    pyclone_out_fl = file.path(pyclonedir, "pyclone_out.tsv")
  )
  # trk = dplyr::mutate(seg_fl = "xx", mut_fl = "xx")

  # contintious path of writing out the trk file!
  # if flowr_path is supplied write it out into the runs folder
  if (!is.null(flowr_path)) {
    pyclonedir_full <- file.path(flowr_path, trk$individual[1], "tmp", pyclonedir)
    wranglr::mkdir(pyclonedir_full)
    write_tsv(trk, file.path(pyclonedir_full, "trk.tsv"))
  }

  trk
}
  

#' this function is invoked from ss_het
#' creates all the commands to run pyclone and subsequent steps
pyclone_pipe_v1_purecn_flo <- function(trk,
                                       samplename,
                                       odir,
                                       flow_conf,
                                       pyclonedir,
                                       bam_readcount_r_func = bam_readcount_r.ssm,
                                       fetch_files.exe = opts_flow$get("fetch_files.exe")) {
  library(pacman)
  p_load(tidyverse, flowr)

  # prep file names, and write out the trk file
  trk <- pyclone_prep_filenms(trk, pyclonedir = pyclonedir, flowr_path = flowr_path)

  # transfer fls
  # trk = readr::read_tsv(trk_fl)
  trk %<>% metadata_for_dnaseq(trk)

  # create all pyclone trk files
  # ** step 1 create pyclone inputs and run pyclone ---------
  cmd_step1 <- glue("funr my.ultraseq::pyclonevi_pipe_r trk={pyclonedir}/trk.tsv flowr_conf={flowr_conf}") %>%
    as.character()
  cmd_step1

  #
  #   # this will ideally, read file from the flowr rundir
  # # create a new trk, with bam readcount files
  # # ** step 2 pyclone ---------
  # # execute: run pyclone (or only parse results)
  # cmd_step2 <- glue("funr my.ultraseq::pyclonevi_pipe_v1_step2_r trk=trk.tsv conf={flow_conf} clean_bams=TRUE force_redo=FALSE") %>%
  #   as.character()
  # # ** step 3 citup ---------
  # cmd_step3 <- glue("funr my.ultraseq::pyclonevi_pipe_v1_step3_r conf={flow_conf}") %>%
  #   as.character()

  cmds <- list(
    pyclone_mrg = cmd_step1
  )
  flowmat_pyclone <- to_flowmat(cmds, samplename)

  # ** transfer out ---------
  # fls_out <- c(
  #   "tables", "plots", "yaml",
  #   "*tsv", # includes trk, bamreadcount, mat, vaf etc.
  #   "citup_vaf", "citup_ccf" # citup stuff
  # )
  # wranglr:::mkdir(odir)
  # cmd_transfer_out <- glue("rsync -avP {fls_out} {odir}/") %>%
  #   as.character()
  # flowmat_transfer_out <- to_flowmat(
  #   list(transfer_out = as.character(cmd_transfer_out)),
  #   samplename = samplename
  # )
  #
  flowmat <- bind_rows(
    # flowmat_transfer_in,
    flowmat_pyclone
    # flowmat_transfer_out
  )

  # these are a final set of outfiles
  # lst_outfiles <- c(
  #   pyclone_table = "tables/loci_ann.tsv",
  #   se_pyclone = "tables/se_pyclone.rds"
  # )

  return(list(
    flowmat = flowmat,
    trk = trk
    # outfiles = lst_outfiles
  ))
}


# ------ run analysis -----
# Here we run allowing for up to 40 clusters (clones), 
# using the Beta-Binomial distribution and performing 10 random restarts. 
# This should take a under five minutes.
pyclonevi.run_analysis_pipeline <- function(df,
                                          oprefix = "",
                                          pyclonevi_exe = "~/apps/conda/3/envs/pyclone-vi/bin/pyclone-vi",
                                          pyclonevi_params = "-c 40 -d beta-binomial -r 10") {

  # combine both files:
  # issue with some mutation are repeated:.default
  # https://github.com/Roth-Lab/pyclone-vi/issues/12
  #   variant_case                          mutation_id        n
  #   <chr>                                 <chr>          <int>
  # 1 IPCT-S4013-MOON0051-Cap2527-4-HTID229 chr:6:86253434     2
  # 2 IPCT-S4014-MOON0051-Cap2544-4-HTID003 chr:6:86253434     2
  # 
  i=1
  df_muts = lapply(1:2, function(i){
    read_tsv(df$pyclone_inp_tsv[i]) %>% 
      mutate(tumour_content = round(df$purity[i], 2))
   }) %>% bind_rows() %>% 
   mutate(major_cn = round(major_cn, 0),
         alt_counts = var_counts)
  write_tsv(df_muts, "pyclone/pyclone_test.tsv")
  # run pyclone
  # mutfls <- paste(df$pyclone_inp_tsv, collapse = " ")
  # tum_contents <- paste(round(df$purity, 2), collapse = " ")
  cmd_pyc <- glue("{pyclonevi_exe} fit -i pyclone/pyclone_test.tsv -o pyclone/test.h5 {pyclonevi_params}")
  system(cmd_pyc)
  cmd_pyc2 = glue("{pyclonevi_exe} write-results-file -i pyclone/test.h5 -o pyclone/test.tsv")
  system(cmd_pyc2)

  # parse loci
  # pyclone_parse_cluster(df, odir)
  # run clonal evol

  # get trees

  # get sig muts

  # run pathway analysis

  return(list(cmd = cmd_pyc))
}
