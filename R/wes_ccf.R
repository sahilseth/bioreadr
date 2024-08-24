
# get df_mut_cnv:
# ideally, we would intersect 
# df_seg, with integer CNV data
# with df_mut
# both would be granges object


# add cnv data for these mutations
# df_snv_f2 = df_snv_f1 %>% 
#   tidylog::left_join(df_cnv_ann_lng2, 
#                      by = c("hugo_symbol" = "gene_name", 
#                             "tumor_sample_barcode" = "Tumor_Sample_Barcode"))

if(FALSE){
  # this sample has 5 mutations (missense mutations)
  df_mut_cnv = df_snv_f2 %>% tidylog::filter(grepl("12011101021114", tumor_sample_barcode)) %>% 
    dplyr::select(gene_name = hugo_symbol, key, 
                  chr = chromosome, start = start_position, 
                  ref = reference_allele, alt = tumor_seq_allele2, vaf = taf, 
                  cnv_state = cnv_state_ploidy_corr)
  df_mut_cnv
  
  out = out_sclust_merged$df_summ %>% tidylog::filter(grepl("12011101021114", sample_name))
  purity = out$purity
  ploidy = out$ploidy
  
  vaf = df_mut_cnv$vaf[1]
  cnv_state = df_mut_cnv$cnv_state[1]
  # 
  ccf_cnv = 1
  
}


# Function to calculate CCF
# @param VAF - Variant allele frequency observed in reads; 
# @param ploidy - ploidy in position of reported variant (optional, default = 2 ). In other words, is this variant together with CNV;
# @param ccf_cnv - Cancer Cell Fraction of this ploidy. For germline CNVs its 1, and for somatic CNVs it can take values in interval (0,1] (optional, default = 1);
# @param purity - purity of cancer tissue and it is value in interval (0,1] but is expected to be high, much closer to 1 then 0.  (optional, default = 1)
# https://www.sciencedirect.com/science/article/pii/S2405471221002842#sec2
# diploid: 
# c = (2*vaf)/purity
# for example: 0.25 vaf; with purity 0.5 --> CCF = 1
# C = ( purity*N_total + (1-purity)*2 ) * vaf /(purity*M)
calc_ccf <- function (vaf, ploidy = 2, ccf_cnv = 1, purity = 1) {
  if (sum(is.na(ploidy))){
    ploidy[is.na(ploidy)] <- 2
  }
  if (sum(is.na(ccf_cnv))){
    ccf_cnv[is.na(ccf_cnv)] <- 1
  }  
  if (sum(is.na(purity))){
    purity[is.na(purity)] <- 1
  }  
  ccf <- ((2 + (ploidy-2)*ccf_cnv)*vaf)/purity
  return(ccf)
}

ccf <- function(df_mut_cnv, purity, ploidy){
  
  # we expect the following columns in the df_mut
  wranglr::expect_columns(df_mut_cnv, c("gene_name", "chr", "start", "ref", "alt", "vaf"))
  
  sample.mutations$CCF <- ccfPloidy(sample.mutations$vaf, sample.mutations$ploidy, sample.mutations$ccf_cnv)
  
  
}