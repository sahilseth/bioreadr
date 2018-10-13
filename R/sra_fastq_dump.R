

if(FALSE){
  
  library(pacman)
  p_load(readr, dplyr, flowr, glue, janitor, tidyr, readxl)
  
  # read in SRR accessions
  wd = "~/rsrch1_data/public_data/tnbc/SRP114962";setwd(wd)
  srr_ids = scan("SRR_Acc_List.txt", what = "character");srr_ids
  df_srr = read_tsv("SraRunTable.txt") %>% clean_names() %>% 
    separate(sample_name, c("sample_name", "cell_id"), "cell", remove = F, fill = "right") %>% 
    mutate(tissue = ifelse(tissue == "tumor", sample_name, tissue),
           tissue = gsub("Blood", "", tissue),
           tissue = gsub("OP", "", tissue),
           tissue = ifelse(nchar(tissue) > 6, substr(tissue, 1, 6), tissue), 
           treatment = gsub("TX", "", treatment),
           treatment = gsub("2cycleschemo", "mid", treatment),
           treatment = gsub("operative", "post", treatment))
  
  df_srr_summ = count(df_srr, assay_type, tissue, library_layout, treatment, library_source)
  df_srr_summ_wd = spread(df_srr_summ, treatment, n)
  
  write_sheet(df_srr_summ, "df_srr_summ.xlsx")  
  write_sheet(df_srr_summ_wd, "df_srr_summ_wd.xlsx")  
  
  
  
  fqpath = "~/rsrch1_data/public_data/tnbc/SRP114962/fastq"
  
  out = sra_fastq_dump(srr_ids, samplename = "dn_fq", fqpath = fqpath)
  out$flowmat %>% head()
  flowdef = to_flowdef(out$flowmat, 
             queue = "transfer", sub_type = "scatter", 
             platform = "lsf", walltime = "24:00", memory_reserved = "8192")
  
  write_tsv(out$flowmat, "flowmat.tsv")
  write_tsv(head(out$flowmat), "flowmat_h.tsv")
  write_tsv(flowdef, "flowdef.tsv")
  
  #' to_flow
  #' flowr to_flow x=flowmat.tsv def=flowdef.tsv execute=TRUE submit=TRUE flow_run_path=.
  #' flowr status x=/rsrch1/iacs/iacs_dep/sseth/data/public_data/tnbc/SRP114962/flowname-dn_fq-20180502-14-02-19-X5YQWR89
  
}

sra_fastq_dump <- function(srr_ids, 
                           samplename = "",
                           fqpath, 
                           fastq_dump_exe = "~/apps/conda/3.5/bin/fastq-dump"){
  
  cmd = glue("{fastq_dump_exe} --outdir {fqpath} --gzip --skip-technical  --readids --read-filter pass --dumpbase --split-3 --clip {srr_ids}")
  cmd
  
  flowmat = to_flowmat(list(dn = cmd), samplename = samplename)
  
  return(list(flowmat = flowmat))
  
}







# END