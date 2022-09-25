

#' to_hla_pvacseq.hlavbseq
#'
#' @param x d4 file
#' @param outfile1 
#'
#' @details more details here:
#' http://nagasakilab.csml.org/hla/
#' 
#' @export
to_pvacseq.hlavbseq <- function(x, 
                                outfile1, outfile2,
                                #allele_db_fl = "~/.rsrch3/home/iacs/sseth/ref/human/b37/annotations/hlavbseq/all_alleles.txt"
                                allele_db_fl = opts_flow$get("hla_allele_db")
){
  pacman::p_load(janitor, readr, dplyr)
  # x="WEX-1004-N_nochr_hla_umap_hlavb_d4.txt"
  # x = "~/.rsrch3/scratch/iacs/sseth/flows/SS/tnbc/ms51_wex/neoantigen/runs/185_003/tmp/IPCT-S2006-MOON0051-Cap2023-8-ID01_190503-A00728-0034-BHJ3W3DMXX-1-ATCACG.bwa_recalibed_hla_umap_hlavb_d4.txt"
  # hla_nethmc_db = read_tsv("nethmc_alleles.txt", col_names = "gene")$gene
  
  # pvac/netmhc recognized alleles
  hla_db = read_tsv(allele_db_fl, col_names = "gene")$gene
  
  df_hla = read_tsv(x, col_types = cols(
    Gene = col_character(),
    Allele1 = col_character(),
    Allele2 = col_character()
  )) %>% clean_names()
  head(df_hla)
  
  # "HLA-A*02:01,HLA-B*35:01",
  hla_pvac_mhc1 = df_hla %>% 
    dplyr::filter(gene %in% LETTERS) %>%
    mutate(allele1 = paste0("HLA-", allele1), 
           allele2 = paste0("HLA-", allele2)) %>% 
    dplyr::select(allele1, allele2) %>% unlist() 
  
  # DRB1*11:01 OR DRB1*11:01-DRB1*11:01
  hla_pvac_mhc2 = df_hla %>% 
    dplyr::filter(!gene %in% LETTERS) %>%
    mutate(allele = paste0(allele1, "-", allele2)) %>% 
    select(allele1, allele2, allele) %>% unlist() 
  
  # filter mhc2:
  hla_pvac_mhc2 = hla_pvac_mhc2[hla_pvac_mhc2 %in% hla_db]
  # hla_pvac = c(hla_pvac_mhc1, hla_pvac_mhc2)
  # # keep those recognized by pvacseq
  # hla_pvac = hla_pvac[hla_pvac %in% hla_db]
  # hla_pvac = hla_pvac %>% paste(collapse = ",")
  
  paste0(hla_pvac_mhc1, collapse = ",") %>% write(outfile1)
  paste0(hla_pvac_mhc2, collapse = ",") %>% write(outfile2)
  
}

to_hla_pvacseq.hlavbseq = to_pvacseq.hlavbseq




# create ref -------
# wget http://nagasakilab.csml.org/hla/hla_all_v2.fasta
# wget http://nagasakilab.csml.org/hla/Allelelist_v2.txt



# END
