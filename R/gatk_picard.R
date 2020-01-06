
#' Use picard's MergeSamFiles tool to merge bam/sam files
#'
#' @description
#' The resulting bam file is sorted and with a indexed created.
#' Validation stringency of inputs is kept as lenient.
#' Multi-threading is turned on by default, though in our experience this does
#' not seem to use a lot of threads.
#' 
#' https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.0.0/picard_sam_MergeSamFiles.php
#'
#' @param bams a vectors of files to merge
#' @param outfile name of the merged bam file
#' @param samplename something
#' @param java_exe path to java
#' @param gatk4_exe something
#' @param gatk4.opts something
#' @param mergesam_opts something
#' @param java_mem java_memory
#'
#' @export
#' 
#' @examples \dontrun{
#' 
#' picard_mergesamfiles(bams = c("a.bam", "b.bam"), "merged.bam", samplename = "ab")
#' 
#' }
picard_mergesamfiles.gatk4 <- function(bams,
                                       outfile,
                                       samplename = opts_flow$get("samplename"),
                                       java_exe = opts_flow$get("java_exe"),
                                       java_mem = opts_flow$get("java_mem_str"),
                                       gatk4_exe = opts_flow$get("gatk4_exe"),
                                       # gatk4.opts = opts_flow$get("gatk4.opts"),
                                       mergesam_opts = "--ASSUME_SORTED --USE_THREADING --CREATE_INDEX"
){
  
  
  
  check_args()
  
  bam_list = paste("--INPUT ", bams, sep = "", collapse = " ")
  cmds = list(merge = glue("{gatk4_exe} --java-options {java_mem} MergeSamFiles {bam_list} --OUTPUT {outfile} {mergesam_opts}"))
  
  # INPUT is a NAMED list
  flowmat = to_flowmat(cmds, samplename)
  return(list(outfiles = outfile, flowmat = flowmat))
  
}




