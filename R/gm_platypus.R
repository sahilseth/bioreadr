
.gm_plat.read <- function(fl){
  
  if(!file.exists(fl)){
    message("file does not exist: ", fl)
    return(data.frame())
  }
  
  df_plat = readr::read_tsv(fl, progress = F, 
           col_types = readr::cols(.default = col_character())) %>% 
    clean_names() 
  message("> filter PASS")
  df_plat %<>% tidylog::filter(filter == "PASS")
  
  # parse a few columns
  df_plat = mutate(df_plat, 
                     start = as.integer(start),
                     end = as.integer(end),
                     
                     t_ref_count = as.integer(t_ref_count),
                     t_alt_count = as.integer(t_alt_count)) %>% 
    mutate(t_depth = (t_alt_count + t_ref_count), 
           t_af = t_alt_count/t_depth)
  
  df_plat
}



#' read mutect
#'
#' @param trk something
#' @param col_fl something
#' @param col_samp something
#'
#' @export
gm_plat.read <- function(trk,
                        col_fl = "mutect",
                        col_samp = "samplelbl1", 
                        cores = 1){
  
  trk = data.frame(trk, stringsAsFactors = F)
  df_mutect = mclapply(1:nrow(trk), function(i){
    message(".", appendLF = F)
    df_mutect = .gm_plat.read(trk[i, col_fl]) %>% 
      mutate(samplename = trk[i, col_samp])
  }, mc.cores = cores) %>% bind_rows()
}


gm_get_aachange <- function(x){
  
  tmp = sapply(x, function(xi){
    tmp = strsplit(xi, ":")[[1]]
    paste0(tmp[1], "_", tail(tmp, 1))
  })
  as.character(tmp)
}

gm_get_clinsig <- function(x){
  
  tmp = sapply(x, function(xi){
    tmp = strsplit(xi, ";")[[1]][[1]]
    gsub("CLINSIG=", "", tmp)
  })
  as.character(tmp)
}

