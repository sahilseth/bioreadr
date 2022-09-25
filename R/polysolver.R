



#' to_pvacseq.polysolver
#' 
#' convert optitype for pvacseq
#'
#' @param path 
#' @param outfile outfile
#'
#' @export
to_pvacseq.polysolver <- function(x, outfile){
  
  # x = "~/.rsrch3/scratch/iacs/sseth/flows/SS/tnbc/ms51_wex/neoantigen/runs/185_003/tmp/polysolver/winners.hla_netmhc.txt"
  df_poly = data.table::fread(x, data.table = F, header = F)
  hla_alleles = df_poly$V1 %>% unlist() %>% paste0(collapse = ",")
  
  write(hla_alleles, outfile)
  outfile
  
}