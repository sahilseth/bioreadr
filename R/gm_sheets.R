

# source('~/projects2/ss_tnbc/analysis/R/split_art_path_id.R')

# x = "IPCT-FC217-MOON0051-Cap1678-6-ID48_180910-STJ00106-0553-AHWNGNBBXX-4-TCGGCA"
gm_split_propel_single_fls <- function(x){
  
  p_load(glue)
  x = x %>% as.character() %>% basename()
  
  # there should be NO other sample_types
  # x = "IPCT-FC217-MOON0051-Cap1678-6-ID48_180910-STJ00106-0553-AHWNGNBBXX-4-TCGGCA"
  #gsub("(.?)_(.*)-([0-9]*)-([ATGC]*|NoIndex)(.?)", "\\1,\\2,\\3,\\4", x)
  if(grepl("\\.", x[1]))
    fmt = "{{proj_samp}}_{{runid}}-{{lane}}-{{index}}\\.{{ext}}"
  else
    fmt = "{{proj_samp}}_{{runid}}-{{lane}}-{{index}}"      
  
  df = ultraseq:::split_names_fastq2(x, 
                                     format = fmt, 
                                     # regex pattern for each piece
                                     lst_patterns = list(
                                       proj_samp = "(.?)",
                                       runid = "(.*)",
                                       lane = "([0-9]*)",
                                       index = "([ATGC]*|NoIndex)",
                                       ext = "(.?)"),
                                     strict_format_checking = FALSE) %>% tbl_df()
  df = mutate(df, oprefix = glue("{proj_samp}_{runid}-{lane}-{index}"))
  
  df
}

gm_split_cgl_fqs <- function(x){
  
  x = x %>% as.character() %>% basename()
  
  # there should be NO other sample_types
  # x = "IPCT-S4009-MOON0051-Cap2373-5-NT96_191001-A00422-0109-AHMCFWDSXX-3-GCCAAGAC.R2.fastq.gz"
  # gsub("(.?)_(.*)-([0-9]*)-([ATGC]*|NoIndex)(.?)", "\\1,\\2,\\3,\\4", x)
  if(grepl("\\.", x[1]))
    fmt = "{{proj_samp}}_{{runid}}-{{lane}}-{{index}}\\.{{read}}\\.{{ext}}"
  else
    fmt = "{{proj_samp}}_{{runid}}-{{lane}}-{{index}}"      
  
  df = ultraseq:::split_names_fastq2(x, 
                                     format = fmt, 
                                     # regex pattern for each piece
                                     lst_patterns = list(
                                       proj_samp = "(.?)",
                                       runid = "(.*)",
                                       lane = "([0-9]*)",
                                       index = "([ATGC]*|NoIndex)",
                                       read = "R([1-2]{1})",
                                       ext = "(.?)"),
                                     strict_format_checking = FALSE) %>% tbl_df()
  df = mutate(df, oprefix = glue("{proj_samp}_{runid}-{lane}-{index}"))
  
  df
}

# IPCT-FC213-MOON0051-Cap1646-6-ID42
gm_split_ipct_name <- function(x){
  
  x = x %>% as.character() %>% basename()
  
  fmt = "{{center}}-{{batch}}-{{project}}-{{samplename}}"
  
  df = ultraseq:::split_names_fastq2(x, 
                                     format = fmt, 
                                     # regex pattern for each piece
                                     lst_patterns = list(
                                       center = "(.?)",
                                       batch = "(.?)",
                                       project = "(.?)",
                                       samplename = "(.*)"),
                                     strict_format_checking = FALSE) %>% tbl_df()

  df
  
  
}






if(FALSE){
  x = "IPCT-S2001-MOON0051-Cap1747-4-HTID004_190103-A00422-0027-AHGMW2DMXX-2-CGCTCCTTCCTT--Cap1710-8-HTID325_190103-A00422-0027-AHGMW2DMXX-1-TCTGGCTTGAAA.011120191547195658.annotated.tsv"
  suffix2=".annotated.tsv"
}

gm_split_propel_paired_fls <- function(x, suffix2 = ""){
  x = basename(x)
  df = ultraseq:::split_names_fastq2(x, 
                                     format = "{{t_proj_samp}}_{{t_runid}}--{{n_proj_samp}}_{{n_runid}}{{suffix}}", 
                                     # regex pattern for each piece
                                     lst_patterns = list(
                                       t_proj_samp = "(.?)",
                                       t_runid = "(.*)",
                                       n_proj_samp = "(.*)",
                                       n_runid = "(.*)", 
                                       suffix = paste0("(", suffix2, ")")),
                                     strict_format_checking = FALSE) %>% tbl_df()
  df = mutate(df, oprefix = glue("{t_proj_samp}_{t_runid}--{n_proj_samp}_{n_runid}"))
  df
}


# 'IPCT-S2001-MOON0051-Cap1747-4-HTID004_190103-A00422-0027-AHGMW2DMXX-2-CGCTCCTTCCTT--Cap1710-8-HTID325_190103-A00422-0027-AHGMW2DMXX-1-TCTGGCTTGAAA.011120191547195658.annotated.tsv'
gm_trk_mutect.ann <- function(mutect_path, recursive = F){
  
  suffix = ".[0-9]?.annotated.tsv"
  suffix2 = "\\.[0-9]*\\.annotated.tsv"
  mutfls = list.files(mutect_path, "", full.names = T, recursive = recursive) %>% 
    grep(suffix, .,  value=TRUE)
  mutfls 
  
  df_trk = gm_split_propel_paired_fls(mutfls, suffix2 = suffix2)
  data.frame(df_trk, fl_full = mutfls, stringsAsFactors = F)
  
}

# 'IPCT-S2001-MOON0051-Cap1747-4-HTID004_190103-A00422-0027-AHGMW2DMXX-2-CGCTCCTTCCTT--Cap1710-8-HTID325_190103-A00422-0027-AHGMW2DMXX-1-TCTGGCTTGAAA.pindel.011120191547213838.All.annotated.tsv'
gm_trk_pindel.ann <- function(pindel_path){
  
  suffix = "[0-9]?.All.annotated.tsv"
  suffix2 = "\\.pindel\\.[0-9]*\\.All\\.annotated.tsv"
  mutfls = list.files(pindel_path, "", full.names = T) %>% 
    grep(suffix, .,  value = TRUE)
  
  df_trk = gm_split_propel_paired_fls(mutfls, suffix2 = suffix2)
  
  data.frame(df_trk, fl_full = mutfls, stringsAsFactors = F)
  
}

# .WES.segList.20190110205550.xls
gm_trk_exomecn.ann <- function(exomcn_path){
  
  suffix = ".WES.segList"
  suffix_regexp = "\\.WES\\.segList\\.[0-9]*\\.xls"
  
  fls = list.files(exomcn_path, "", full.names = T) %>% 
    grep(suffix, .,  value = TRUE)
  head(fls)
  
  df_trk = gm_split_propel_paired_fls(fls, suffix2 = suffix_regexp)
  head(df_trk)
  
  data.frame(df_trk, fl_full = fls, stringsAsFactors = F)
  
}

# IPCT-S2001-MOON0051-Cap1710-8-HTID289_190103-A00422-0027-AHGMW2DMXX-1-ACGACGTCACGG.bwa_recalibed.bam.platypus.011020191547178559.annotated.tsv
# IPCT-S2001-MOON0051-Cap1710-8-HTID289_190103-A00422-0027-AHGMW2DMXX-1-ACGACGTCACGG.bwa_recalibed.bam.platypus.vcf
# .WES.segList.20190110205550.xls
# platypus_path = "/rsrch2/iacs/ngs_runs/190103_A00422_0027_AHGMW2DMXX/platypus"
gm_trk_platypus.ann <- function(platypus_path){
  p_load(tidyverse)
  # for searching
  suffix = ".platypus." 
  # to be removed to get oprefix
  # suffix_regexp = "\\.bwa_recalibed.bam.platypus.\\.[0-9]*\\.annotated.tsv"
  
  fls = list.files(platypus_path, "tsv", full.names = T) %>% 
    grep(suffix, .,  value = TRUE)
  head(basename(fls))
  
  df_trk = gm_split_propel_single_fls(fls)
  # add VCF file as well
  df_trk %<>% mutate(fl_vcf = paste0(oprefix, ".bwa_recalibed.bam.platypus.vcf"))
  head(df_trk) %>% data.frame()
  
  data.frame(df_trk, fl_full = fls, stringsAsFactors = F)
  
}


# IPCT-S2001-MOON0051-Cap1747-4-HTID004_190103-A00422-0027-AHGMW2DMXX-2-CGCTCCTTCCTT--Cap1710-8-HTID325_190103-A00422-0027-AHGMW2DMXX-1-TCTGGCTTGAAA-20190228-15-08-49-igYHtbvu_fit.rds'
gm_trk_facets.fit <- function(facets_path){
  
  suffix = "_fit.rds"
  # -20190228-15-08-49-igYHtbvu_fit.rds
  # this will work for the next 100 years!
  suffix_regexp = "\\-20[0-9,A-Z,a-z,\\-]*\\_fit.rds"
  
  fls = list.files(facets_path, "", full.names = T) %>% 
    grep(suffix, .,  value = TRUE)
  head(fls)
  
  df_trk = gm_split_propel_paired_fls(fls, suffix2 = suffix_regexp)
  head(df_trk)
  
  data.frame(df_trk, fl_full = fls, stringsAsFactors = F)
  
}

# IPCT-S2001-MOON0051-Cap1747-4-HTID004_190103-A00422-0027-AHGMW2DMXX-2-CGCTCCTTCCTT--Cap1710-8-HTID325_190103-A00422-0027-AHGMW2DMXX-1-TCTGGCTTGAAA-20190228-15-08-49-igYHtbvu_fit.rds'
# IPCT-S2001-MOON0051-Cap1747-4-HTID004_190103-A00422-0027-AHGMW2DMXX-2-CGCTCCTTCCTT--Cap1710-8-HTID325_190103-A00422-0027-AHGMW2DMXX-1-TCTGGCTTGAAA-20190228-15-08-49-igYHtbvu_fit_cncf.tsv
gm_trk_facets.seg <- function(facets_path){
  
  suffix = "_fit_cncf.tsv"
  # -20190228-15-08-49-igYHtbvu_fit.rds
  # this will work for the next 100 years!
  suffix_regexp = "\\-20[0-9,A-Z,a-z,\\-]*\\_fit_cncf.tsv"
  
  fls = list.files(facets_path, "", full.names = T) %>% 
    grep(suffix, .,  value = TRUE)
  head(fls)
  
  df_trk = gm_split_propel_paired_fls(fls, suffix2 = suffix_regexp)
  head(df_trk)
  
  data.frame(df_trk, fl_full = fls, stringsAsFactors = F)
  
}





gm_read_tnbc_samplesheet <- function(fl = "~/projects/samplesheets/tnbc/seq_tracker.xlsx"){
  df_trk = read_sheet(fl, start_row = 2)
  # this may not be required later
  df_sample_type = split_art_pathid(df_trk$art_path_id, col_nm = "art_path_id") %>% unique()
  df_trk_ann = left_join(df_trk, df_sample_type, by = c("art_path_id" = "art_path_id"))
  df_trk_ann
}






.get_bam <- function(propel_single_id, runid, verbose = F){
  runpath1 = "~/.rsrch1/genomic_med/omics/Womens/MS51"
  runpath2 = "~/.rsrch2/iacs/ngs_runs"
  runpath3 = "~/.rsrch1/iacs/ngs_runs"
  
  bam1 = glue("{runpath1}/bam/{propel_single_id}.bwa_recalibed.bam")
  bam2 = glue("{runpath2}/{runid}/bams/{propel_single_id}.bwa_recalibed.bam")
  bam3 = glue("{runpath3}/{runid}/bams/{propel_single_id}.bwa_recalibed.bam")
  
  # for rnaseq
  bam4 = glue("{runpath1}/bam/{propel_single_id}.star_marked.bam")
  
  if(verbose)
    message(bam1, bam2, bam3)
  
  if(file.exists(bam1)){
    bam = bam1
  }else if(file.exists(bam2)){
    bam = bam2
  }else if(file.exists(bam3)){
    bam = bam3
  }else if(file.exists(bam4)){
    bam = bam4
  }else{
    bam = ""
  }
  
  as.character(bam)
}


# lets start fetching files
fetch_fls <- function(df_samp, 
                      col_oprefix = "oprefix", 
                      col_runid = "runid"){
  
  df_samp = data.frame(df_samp, stringsAsFactors = F)
  # for each file
  tmp = lapply(1:nrow(df_samp), function(i){
    
    propel_single_id = df_samp[, col_oprefix][i]
    runid = df_samp[, col_runid][i]
    runid = gsub("-", "_", runid)
    
    bam = .get_bam(propel_single_id, runid)
    #         mutect = .get_mutect(propel_single_id, runid)
    #         exomecn = .get_exomecn(propel_single_id, runid)
    
    data.frame(bam = bam)
    
  }) %>% bind_rows()
  bind_cols(df_samp, tmp)
}


# using this: 
# x="IPCT-S2001-MOON0051-Cap1747-4-HTID004_190103-A00422-0027-AHGMW2DMXX-2-CGCTCCTTCCTT--Cap1710-8-HTID325_190103-A00422-0027-AHGMW2DMXX-1-TCTGGCTTGAAA-20190228-15-08-49-igYHtbvu_fit_cncf.tsv"
# GET:
# IPCT-S2001-MOON0051-Cap1747-4-HTID004_190103-A00422-0027-AHGMW2DMXX-2-CGCTCCTTCCTT--Cap1710-8-HTID325_190103-A00422-0027-AHGMW2DMXX-1-TCTGGCTTGAAA-20190228-15-08-49-igYHtbvu_fit_cncf.tsv
# IPCT-S2001-MOON0051-Cap1747-4-HTID004_190103-A00422-0027-AHGMW2DMXX-2-CGCTCCTTCCTT--Cap1710-8-HTID325_190103-A00422-0027-AHGMW2DMXX-1-TCTGGCTTGAAA-20190228-15-08-49-igYHtbvu_fit_genomewide.pdf
# IPCT-S2001-MOON0051-Cap1747-4-HTID004_190103-A00422-0027-AHGMW2DMXX-2-CGCTCCTTCCTT--Cap1710-8-HTID325_190103-A00422-0027-AHGMW2DMXX-1-TCTGGCTTGAAA-20190228-15-08-49-igYHtbvu_fit.rds
# IPCT-S2001-MOON0051-Cap1747-4-HTID004_190103-A00422-0027-AHGMW2DMXX-2-CGCTCCTTCCTT--Cap1710-8-HTID325_190103-A00422-0027-AHGMW2DMXX-1-TCTGGCTTGAAA-20190228-15-08-49-igYHtbvu_oo.rds
# IPCT-S2001-MOON0051-Cap1747-4-HTID004_190103-A00422-0027-AHGMW2DMXX-2-CGCTCCTTCCTT--Cap1710-8-HTID325_190103-A00422-0027-AHGMW2DMXX-1-TCTGGCTTGAAA-20190228-15-08-49-igYHtbvu_xx.rds
get_facet_files <- function(x){
  
  
  fl_plot_genomewide = gsub("_fit_cncf.tsv$", "_fit_genomewide.pdf", x);#fl_plot_genomewide
  fl_fit = gsub("_fit_cncf.tsv$", "_fit.rds", x);fl_fit
  fl_oo = gsub("_fit_cncf.tsv", "_oo.rds", x);fl_plot_genomewide
  fl_xx = gsub("_fit_cncf.tsv", "_xx.rds", x);fl_plot_genomewide
  
  fl_prefix = gsub("_fit_cncf.tsv$", "", x)
  
  list(
    # fl_prefix = fl_prefix, 
    fl_plot_genomewide = fl_plot_genomewide,
    fl_fit = fl_fit, 
    fl_oo = fl_oo, 
    fl_xx = fl_xx,
    fl_cncf = x)
}







# END
