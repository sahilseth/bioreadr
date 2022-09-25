


#' to_pvacseq.optitype
#' 
#' convert optitype for pvacseq
#'
#' @param path 
#' @param outfile outfile
#'
#' @export
to_pvacseq.optitype <- function(path, outfile){
  
  # path = "~/.rsrch3/scratch/iacs/sseth/flows/SS/tnbc/ms51_wex/neoantigen/runs/185_003/tmp/optitype"
  # fl = list.files(path, pattern = "result.tsv", recursive = T) %>% tail(1)
  fl = fs::dir_ls(path, glob = "*/*_result.tsv", recurse = T) %>% tail(1)
  df_opti = data.table::fread(fl, data.table = F)
  hla_alleles =   dplyr::select(df_opti, A1:C2) %>% unlist() %>% 
    paste0("HLA-", ., collapse = ",") %>% 
    mutate(allele1 = hla_trim_d4(allele1),
           allele2 = hla_trim_d4(allele2))
  
  write(hla_alleles, outfile)
  outfile

}