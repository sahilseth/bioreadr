.get_bam <- function(propel_single_id, runid = "", verbose = F, ext = ".bwa_recalibed.bam"){
  runpath1.1 = "/rsrch1/genomic_med/omics/Womens/MS51/bam"
  runpath1.2 = "/rsrch1/genomic_med/omics/Womens/MS51/bams"
  runpath2 = "/rsrch3/scratch/genomic_med/omics/Womens/MS51/bams"
  runpath3 = "/rsrch3/scratch/iacs/sseth/flows/SS/tnbc/ms51_wex/ss_cl_het_gb/bam"
  runpath4 = "/rsrch2/iacs/ngs_runs"
  runpath5 = "/rsrch1/iacs/ngs_runs"

  bam1.1 = glue("{runpath1.1}/{propel_single_id}{ext}")
  bam1.2 = glue("{runpath1.2}/{propel_single_id}{ext}")
  bam2 = glue("{runpath2}/{propel_single_id}{ext}")
  bam3 = glue("{runpath3}/{propel_single_id}{ext}")
  bam4 = glue("{runpath4}/{runid}/bams/{propel_single_id}{ext}")
  bam5 = glue("{runpath3}/{runid}/bams/{propel_single_id}{ext}")
  
  if(verbose)
    message(bam1.1, bam1.2, bam2, bam3, bam4, bam5)
  
  if(file.exists(bam1.1)){
    bam = bam1.1
  }else if(file.exists(bam1.2)){
    bam = bam1.2
  }else if(file.exists(bam2)){
    bam = bam2
  }else if(file.exists(bam3)){
    bam = bam3
  }else if(file.exists(bam4)){
    bam = bam4
  }else if(file.exists(bam5)){
    bam = bam5
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

