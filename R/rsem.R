
# Following tutorial from Bo Li
# https://github.com/bli25ucb/RSEM_tutorial




# x = "/Users/sahilseth/rsrch1_data/public_data/tnbc/art/rnaseq/rsem/IPCT-FC217-MOON0051-Cap1678-6-ID48_180910-STJ00106-0553-AHWNGNBBXX-4-TCGGCA.genes.results"
read_rsem <- function(x, reader = read_tsv, ...){
  df = reader(x, ...) %>% tbl_df() %>% clean_names()
  
  # summarize gene level counts
  # df = separate(df, col = target_id, 
  #               into = c("transcript_id", "gene_id"), sep = "\\|", remove = FALSE, extra = "drop")
  # head(df) %>% as.data.frame()
  
  # sum est_counts and transcripts
  df = group_by(df, gene_id) %>%
    mutate(gene_expected_count = sum(expected_count), 
           gene_tpm = sum(tpm), 
           gene_fpkm = sum(fpkm))
  df
}





