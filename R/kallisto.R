


#' parse_kallisto_metrics
#'
#' @param pattern something 
#' @param wd path for flowr dir
parse_kallisto_metrics <- function(wd, pattern = "kallisto"){
  
  if(missing(pattern))
    kal.dir <- grep("kallisto", list.dirs(wd, recursive = FALSE, full.names = FALSE), value = TRUE)
  
  aln.dir = file.path(wd, aln.dir)
  
  if(!file.exists(aln.dir)){
    message("aln.dir does not exist: ", basename(aln.dir))
    return()
  }
  if(length(list.files(aln.dir, pattern = "out"))==0){
    message("aln dir does not have *.out files yet, still running ?", basename(aln.dir))
    return()
  }
  
  message("Using aln.dir: ", basename(aln.dir))
  if(missing(files))
    files = list.files(aln.dir, pattern = ".out", full.names = TRUE)
  else
    message("Files supplied, will use those instead")
}


#' Title
#'
#' @param x one metric file
.parse_kallisto_metrics <- function(x){
  
  if(!file.exists(x))
    return(NA)
  
  y = scan(x, what = "character", sep = "\n", quiet = TRUE)
  
  # extract relevant columns
  y1 = grep("[quant]", y, value = T, fixed = T)
  y2 = grep("[index]", y, value = T, fixed = T)
  y3 = grep("[   em]", y, value = T, fixed = T)
  y = c(y1, y2, y3)
  
  y = stringr:::str_trim(y) ## trim whitespace
  
  # extract specific values
  vars = c(kmer_length = "k-mer length: ", 
           est_avg_frag_length = "estimated average fragment length: ",
           num_kmers = "number of k-mers: ")
  
  tmp = lapply(seq_along(vars), function(i){
    var = vars[i]
    df = grep(var, y, value = T) %>% strsplit(":") %>% 
      do.call(rbind, .)
    df[, 2] = stringr::str_trim(df[, 2])
    df[, 1] = names(var)
    data.frame(df, stringsAsFactors = F)
  }) %>% bind_rows()
  names(tmp) = c("var", "value")
  
  tmp2 = gsub(".* processed (.*) reads, (.*) reads pseudoaligned", "\\1;\\2", y[5]) %>% 
    strsplit(";")
  tmp2 = data.frame(var = c("total_reads", "pseudoaligned"), value = tmp2[[1]], stringsAsFactors = F)
  
  df_metrics = bind_rows(tmp2, tmp)
  
  df_metrics = mutate(df_metrics, 
                      value = gsub(",", "", value), 
                      value = as.numeric(value))
  df_metrics
}












# END

