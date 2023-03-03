# ** ens ann -----

# ens_query = "Ensembl 91 EnsDb for Homo sapiens", ah_id = "AH60773"
get_ens_df <- function(ens_query = NULL, ah_id = NULL){
  p_load(AnnotationHub)
  ah = AnnotationHub()
  # snapshotDate(ah) = "2021-02-24"
  human_ens <- query(ah, c("Homo sapiens", "EnsDb"))
  human_ens
  
  # q_ens = query(ah, ens_query)
  # print(q_ens)
  
  # ah_id = "AH60773"
  # ah_id = q_ens$ah_id

  if(!is.null(ah_id))
    ens <- human_ens[[ah_id]]
  if(is.null(ah_id)){
    q_ens = query(ah, ens_query)
    ens <- human_ens[[q_ens$ah_id]]
  }
  
  
  # extract gene level data
  # Extract gene-level information
  df_ens = genes(ens, return.type = "data.frame") %>% 
    data.frame(row.names = .$gene_id)
   df_ens 
}