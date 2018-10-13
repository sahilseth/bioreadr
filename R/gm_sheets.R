# x = '185_002_SG'
# x = '185_002_X0_P1'
#ultraseq:::
# source('~/Dropbox/public/flow-r/ultraseq/ultraseq/R/sheets_fastqs.R')
# split_names_fastq2(x, 
#                               "{{trialid}}_{{path_patient_id}}_{{sample_type}}_{{pdx_passage}}", 
#                               lst_patterns = list(trialid = ".{3}", 
#                                                   path_patient_id = ".{3}",
#                                                   sample_type = ".{3}",
#                                                   pdx_passage = ".?"))
# 
# library(whisker)
# library(flowr)
gm_split_art_pathid <- function(x, col_nm = "art_pathid", sep = "_"){
  
  # there should be NO other sample_types
  df = data.frame(x = x, stringsAsFactors = F) %>% 
    tidyr::separate(col = x, into = c("trialid", "path_patient_id", "sample_type", "pdx_passage"), 
             sep = sep, fill = "right", remove = F) %>% 
    dplyr::mutate(sample_type = factor(sample_type, levels = c("SG", "BG", "T0", "T1", "TS", "X0", "X1", "XS")))
  colnames(df)[1] = col_nm
  df
}


gm_split_propel_single_fls <- function(x){
  
  x = x %>% as.character() %>% basename()
  
  # there should be NO other sample_types
  # x = "IPCT-FC217-MOON0051-Cap1678-6-ID48_180910-STJ00106-0553-AHWNGNBBXX-4-TCGGCA"
  #gsub("(.?)_(.*)-([0-9]*)-([ATGC]*|NoIndex)(.?)", "\\1,\\2,\\3,\\4", x)
  df = ultraseq:::split_names_fastq2(x, 
                                     format = "{{proj_samp}}_{{runid}}-{{lane}}-{{index}}\\.{{ext}}", 
                                     # regex pattern for each piece
                                     lst_patterns = list(
                                       proj_samp = "(.?)",
                                       runid = "(.*)",
                                       lane = "([0-9]*)",
                                       index = "([ATGC]*|NoIndex)",
                                       ext = "(.?)"),
                                     strict_format_checking = FALSE) %>% tbl_df()
  df
  
  
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
  df
}


gm_trk_mutect.ann <- function(mutpath){
  
  suffix = ".[0-9]?.annotated.tsv"
  suffix2 = "\\.[0-9]*\\.annotated.tsv"
  mutfls = list.files(mutpath, "", full.names = T) %>% 
    grep(suffix, .,  value=TRUE)
  mutfls 
  
  df_trk = gm_split_propel_paired_fls(mutfls, suffix2 = suffix2)
  data.frame(df_trk, fl_full = mutfls, stringsAsFactors = F)
  
}


gm_read_tnbc_samplesheet <- function(fl = "~/projects/samplesheets/tnbc/seq_tracker.xlsx"){
  df_trk = read_sheet(fl, start_row = 2)
  # this may not be required later
  df_sample_type = split_art_pathid(df_trk$art_path_id, col_nm = "art_path_id") %>% unique()
  df_trk_ann = left_join(df_trk, df_sample_type, by = c("art_path_id" = "art_path_id"))
  df_trk_ann
}









# END
