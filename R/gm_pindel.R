

#' .gm_pindel.read
#'
#' @param fl read a single pindel file
#'
#' @return
#' @export
#'
#' @examples
.gm_pindel.read <- function(fl, fix_deletions = T){
  
  if(!file.exists(fl)){
    message("file does not exist: ", fl)
    return(data.frame())
  }
    
    
  df_pindel = readr::read_tsv(fl, progress = F, col_types = cols(.default = col_character())) %>% 
    janitor::clean_names()
  colnames(df_pindel)
  
  df_pindel = mutate(df_pindel, 
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

  
}


#' read pindel
#'
#' @param trk something
#' @param col_fl something
#' @param col_samp something
#'
#' @export
gm_pindel.read <- function(trk,
                        col_fl = "pindel_fl",
                        col_samp = "NAME", 
                        cores = 1){
  
  trk = data.frame(trk, stringsAsFactors = F)
  df_pindel = mclapply(1:nrow(trk), function(i){
    message(".", appendLF = F)
    df_pindel = .gm_pindel.read(trk[i, col_fl]) %>% 
      mutate(samplename = trk[i, col_samp])
  }, mc.cores = cores) %>% bind_rows()
  

  # create intersected bed
  df_pindel
}


gm_get_aachange <- function(x){
  
  tmp = sapply(x, function(xi){
    tmp = strsplit(xi, ":")[[1]]
    paste0(tmp[1], "_", tail(tmp, 1))
  })
  as.character(tmp)
}

