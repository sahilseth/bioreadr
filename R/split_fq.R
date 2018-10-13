

#' split fastq into pre-determined number of splits
#' 
#' @description 
#' 
#' This calls a internal module
#' 
#' @param x file to split
#' @param chunks approx size of resulting file, in MB
#'
#' @export
#'
#'
#' @details
#' estimate the number of splits required
#' pre-compute the number of splits it would generate
#' make commands to make the split
#' approx, make 200mb files, for splitting...
split_fq <- function(x, samplename, jobname,
  split_fq_exe = system.file("scripts/split_fq.sh", package = "ultraseq"),
  chunks = 200){

  outbase = gsub("fastq.gz$|fastq$", "", x)
  splt_out = sprintf("%s%03d", outbase, 0:(chunks-1))
  final_out = paste0(splt_out, ".fastq")

  cmds = list(split_fq = sprintf("%s %s %s",
    split_fq_exe, chunks, x))

  if(!missing(jobname))
    names(cmds) = jobname
  
  flowmat = to_flowmat(cmds, samplename = samplename)
  return(list(outfiles = final_out, flowmat = flowmat))
}




## fqz=/rsrch2/iacs/ngs_runs/Y76I6Y76/Project_LTP9X-AMLLR/Sample_Futreal-LTP9X-AMLLR-240656/Futreal-LTP9X-AMLLR-240656_CTTGTA_L001_R1_001.fastq.gz
## perl ~/Dropbox/public/github_ngsflows/inst/scripts/fastq-splitter.pl --n-parts 278
##
##
##
##
##


if(FALSE){
  .detect_chunk_size <- function(){
    sz = file.size(x[1])/10^6
    
    ## if more than 1gb, only then splitting is worth it !
    if(split_size*2 > sz){
      message("No need to split")
      return()
    }
    chunks = round(sz/split_size)
    
  }
  
  x = "/rsrch1/iacs/tmp/illumina_platinum/50x/NA12877/ERR194146.fastq.gz"
  split_fq(x)

}


