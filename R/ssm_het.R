het_ssm.read <- function(trk,
                           col_fl = "pcgr_combvcf_rds",
                           col_samp = "name", 
                           cores = 1){
  
  trk = data.frame(trk, stringsAsFactors = F)
  i=1
  df_mutect = mclapply(1:nrow(trk), function(i){
    message(".", appendLF = F)
    df_ssm = read_rds(trk[i, col_fl]) %>% 
      mutate(samplename = trk[i, col_samp])
  }, mc.cores = cores) %>% bind_rows()
}
