gm_coverage.read <- function(x){

  # read the file
  df_trk = read_tsv(x) %>% clean_names()

}
