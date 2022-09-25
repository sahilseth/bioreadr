

#' to_pvacseq.phlat
#' 
#' convert optitype for pvacseq
#'
#' @param path 
#' @param outfile1 outfile
#' @param outfile2 outfile
#'
#' @export
to_pvacseq.phlat <- function(phlat_fl, outfile1, outfile2){
  
  # setwd("~/.rsrch3/scratch/iacs/sseth/flows/SS/tnbc/ms51_wex/neoantigen/runs/185_003/tmp")
  # list.files(path)
  # fl = list.files(path, pattern = "result.tsv", recursive = T) %>% tail(1)
  df_phlat = data.table::fread(phlat_fl, data.table = F) %>% 
    mutate(gene = gsub("HLA_", "", Locus)) %>% 
    clean_names() %>% 
    mutate(allele1 = hla_trim_d4(allele1),
           allele2 = hla_trim_d4(allele2))
  
  hla_pvac_mhc1 = df_phlat %>% 
    dplyr::filter(gene %in% LETTERS) %>%
    mutate(allele1 = paste0("HLA-", allele1), 
           allele2 = paste0("HLA-", allele2)) %>% 
    select(allele1, allele2) %>% unlist() 
  
  # DRB1*11:01 OR DRB1*11:01-DRB1*11:01
  hla_pvac_mhc2 = df_phlat %>% 
    dplyr::filter(!gene %in% LETTERS) %>%
    mutate(allele = paste0(allele1, "-", allele2)) %>% 
    select(allele1, allele2, allele) %>% unlist() 
  
  # hla_pvac = c(hla_pvac_mhc1, hla_pvac_mhc2)
  
  paste0(hla_pvac_mhc1, collapse = ",") %>% write(outfile1)
  paste0(hla_pvac_mhc2, collapse = ",") %>% write(outfile2)
  
  c(outfile1, outfile2)

}