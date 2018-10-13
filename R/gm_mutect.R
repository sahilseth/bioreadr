




#' read mutect
#'
#' @param muts 
#'
#' @return
#' @export
#'
#' @examples
mutect.read <- function(trk,
                        col_fl = "mutect",
                        col_samp = "samplelbl1"){
  
  trk = data.frame(trk, stringsAsFactors = F)
  df_mutect = lapply(1:nrow(trk), function(i){
    message(".", appendLF = F)
    df_mutect = read_tsv(trk[i, col_fl], progress = F, col_types = cols(.default = col_character())) %>% 
      clean_names() %>% 
      mutate(samplename = trk[i, col_samp])
    df_mutect = filter(df_mutect, covered == "COVERED")
  }) %>% bind_rows()
  
  # parse a few columns
  df_mutect = mutate(df_mutect, 
                     start = as.integer(start),
                     end = as.integer(end),
                     
                     t_ref_count = as.integer(t_ref_count),
                     t_alt_count = as.integer(t_alt_count),
                     
                     n_ref_count = as.integer(n_ref_count),
                     n_alt_count = as.integer(n_alt_count),
                     
                     tumor_f = as.numeric(tumor_f)) %>% 
    mutate(t_depth = (t_alt_count + t_ref_count), 
           t_af = t_alt_count/t_depth,
           
           n_depth = (n_alt_count + n_ref_count), 
           n_af = n_alt_count/n_depth)
  
  # create intersected bed
}


gm_get_aachange <- function(x){
  
  tmp = sapply(x, function(xi){
    tmp = strsplit(xi, ":")[[1]]
    paste0(tmp[1], "_", tail(tmp, 1))
  })
  as.character(tmp)
}

