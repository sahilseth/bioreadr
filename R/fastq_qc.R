

#' Custom fastqc (R)
#'
#' @param fq 
#' @param outfile optional outfile
#' 
#' @import Biostrings ShortRead
#'
fastq_qc <- function(x, outfile){
  message("working on ", x)
  
  # read fq
  fq = ShortRead::readFastq(x)
  
  # get consensus
  out = Biostrings::consensusMatrix(sread(fq), as.prob = TRUE)
  consensus_seq = Biostrings::consensusString(out)
  
  # get avg quality
  qual = as(quality(fq), "matrix")
  qual = round(colMeans(qual, na.rm=TRUE), 0)
  
  df_qual = tibble(position = 1:length(qual), qual = qual)
  
  # get other metrics (WIP)
  
  ret1 = list(consensus_seq = consensus_seq, df_qual)
  as.data.frame(ret1) %>% write_tsv(ret1, outfile)
  
  invisible(ret1)
}

if(FALSE){
  
  x = ""
}
