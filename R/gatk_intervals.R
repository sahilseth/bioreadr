# given a fasta file, get the intervals  
#
#


get_intervals <- function(ref_fasta = opts_flow$get("ref_fasta_path"), 
                          type = c("gatk")){
  chrs_info <- get_fasta_chrs(ref_fasta)
  intervals_opts = paste0(" -L ", chrs_info) # interval files
  list(intervals_opts = intervals_opts)
}