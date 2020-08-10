.get_bam <- function(propel_single_id, runid, verbose = F){
  runpath1 = "~/.rsrch1/genomic_med/omics/Womens/MS51"
  runpath2 = "~/.rsrch2/iacs/ngs_runs"
  runpath3 = "~/.rsrch1/iacs/ngs_runs"
  
  bam1 = glue("{runpath1}/bam/{propel_single_id}.bwa_recalibed.bam")
  bam2 = glue("{runpath2}/{runid}/bams/{propel_single_id}.bwa_recalibed.bam")
  bam3 = glue("{runpath3}/{runid}/bams/{propel_single_id}.bwa_recalibed.bam")
  
  if(verbose)
    message(bam1, bam2, bam3)
  
  if(file.exists(bam1)){
    bam = bam1
  }else if(file.exists(bam2)){
    bam = bam2
  }else if(file.exists(bam3)){
    bam = bam3
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

