
#' picard_rg
#' 
#' @export
#' 
#' @importFrom tools file_path_sans_ext
#' 
picard_rg <- function(x, 
                      samplename = opts_flow$get("samplename"),
                      ## convert these into get option also, only for this flow
                      rg_lane = opts_flow$get("rg_lane"),
                      rg_platform = opts_flow$get("rg_platform"),
                      rg_center = opts_flow$get("rg_center"),
                      
                      java_exe = opts_flow$get("java_exe"),
                      java_mem = opts_flow$get("java_mem"),
                      java_tmp = opts_flow$get("java_tmp"),
                      picard_jar = opts_flow$get("picard_jar")){
  
  
  check_args()

  ## make this editable later ....
  rgid = rglb = rgsm = samplename
  rgpu = rg_lane
  
  ## add RG to the orignal bam name
  bamrg_files = sprintf("%s_rg.bam", tools::file_path_sans_ext(x))
  cmds = list(fixrg = sprintf("%s %s -Djava.io.tmpdir=%s -jar %s AddOrReplaceReadGroups INPUT=%s OUTPUT=%s SORT_ORDER=coordinate RGID='%s' RGLB='%s' RGPL='%s' RGPU='%s' RGSM='%s' RGCN='%s' VALIDATION_STRINGENCY=LENIENT",
                              java_exe, java_mem, java_tmp, picard_jar, 
                              x, bamrg_files, rgid, rglb, 
                              rg_platform, rgpu, rgsm, rg_center))
  
  flowmat = to_flowmat(cmds, samplename)
  ret = list(outfiles = bamrg_files, flowmat = flowmat)
  return(ret)
  
}

#' Use picard's MergeSamFiles tool to merge bam/sam files
#'
#' @description
#' The resulting bam file is sorted and with a indexed created.
#' Validation stringency of inputs is kept as lenient.
#' Multi-threading is turned on by default, though in our experience this does
#' not seem to use a lot of threads.
#'
#' @param bams a vectors of files to merge
#' @param outfile name of the merged bam file
#' @param samplename 
#' @param java_exe path to java
#' @param java_mem java memory
#' @param java_tmp
#' @param picard_jar
#'
#' @export
#' 
#' @examples \dontrun{
#' 
#' picard_mergesamfiles(bams = c("a.bam", "b.bam"), "merged.bam", samplename = "ab")
#' 
#' }
picard_mergesamfiles <- function(bams,
                                 outfile,
                         samplename = opts_flow$get("samplename"),
                         java_exe = opts_flow$get("java_exe"),
                         java_mem = opts_flow$get("java_mem"),
                         java_tmp = opts_flow$get("java_tmp"),
                         picard_jar = opts_flow$get("picard_jar"), 
                         mergesam_opts = "ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true USE_THREADING=true"
                         ){
  
  
  
  check_args()
  
  bam_list = paste("INPUT=", bams, sep = "", collapse = " ")
  cmds = list(merge = sprintf("%s %s -Djava.io.tmpdir=%s -jar %s MergeSamFiles %s OUTPUT=%s %s",
                              java_exe, java_mem, java_tmp, picard_jar, bam_list, outfile, mergesam_opts))
  
  # INPUT is a NAMED list
  flowmat = to_flowmat(cmds, samplename)
  return(list(outfiles = outfile, flowmat = flowmat))
  
}


#' Title
#'
#' @param bam
#' @param mergedbam
#' @param samplename
#' @param java_exe
#' @param java_mem
#' @param java_tmp
#' @param picard_dir
#'
#' @export
#' 
picard_samtofastq <- function(bam,
                             samplename = opts_flow$get("samplename"),
                             paired = TRUE,
                             split = FALSE,
                             num_reads = 8000000,
                             split_fq_exe = system.file("scripts/split_fq", package = "ngsflows"),
                             java_exe = opts_flow$get("java_exe"),
                             java_mem = opts_flow$get("java_mem"),
                             java_tmp = opts_flow$get("java_tmp"),
                             picard_dir = opts_flow$get("picard_dir"),
                             bam_fastq_opts = "INCLUDE_NON_PF_READS=true INCLUDE_NON_PRIMARY_ALIGNMENTS=false RE_REVERSE=true VALIDATION_STRINGENCY=LENIENT"){
  
  check_args()
  
  
  if(paired){
    fq1 = gsub(".bam$", "_1.fastq", bam)
    fq2 = gsub(".bam$", "_2.fastq", bam)
    fq3 = gsub(".bam$", "_unpaired.fastq", bam)
    fqs = list(fq1 = fq1, fq2 = fq2, fq3 = fq3)
    fqout = sprintf("FASTQ=%s SECOND_END_FASTQ=%s UNPAIRED_FASTQ=%s",
                    fq1, fq2, fq3)
    # we would stream picard into creating splitted fq files. A second flow would read these and start the next step.
    
  }else{
    fq1 = gsub(".bam$", ".fastq", bam)
    fqout = paste0("FASTQ=", fq1)
    fqs = list(fq1 = fq1)
  }
  
  if(split){
    fq = gsub(".bam$", "", bam)
    fqout = sprintf("FASTQ=/dev/stdout INTERLEAVE=true | bash %s -n %s -f /dev/stdin -o %s",
                    split_fq_exe, as.integer(num_reads), fq)
  }
  
  
  bam_fastq = sprintf("%s %s -Djava.io.tmpdir=%s -jar %s/picard.jar SamToFastq INPUT=%s %s %s",
                      java_exe, java_mem, java_tmp, picard_dir, bam, bam_fastq_opts, fqout)
  
  
  cmd = list(bam_fastq = bam_fastq)
  #system(unlist(cmd))
  
  # INPUT is a NAMED list ------
  flowmat = to_flowmat(cmd, samplename)
  return(list(flowmat = flowmat, outfiles = fqs))
}

#' Use picard's reorder bam file
#'
#' @description
#' The resulting file is sorted and index is created for it.
#' Validation stringency of inputs is kept as lenient.
#' Multi-threading is turned on by default, though in our experience this does
#' not seem to use a lot of threads.
#'
#' @param x a vectors of files to merge
#' @param mergedbam
#' @param samplename
#' @param java_exe
#' @param java_mem
#' @param java_tmp
#' @param picard_dir
#'
#' @export
#' 
picard_reorder <- function(x, outfile,
                           samplename = opts_flow$get("samplename"),
                           java_exe = opts_flow$get("java_exe"),
                           java_mem = opts_flow$get("java_mem"),
                           java_tmp = opts_flow$get("java_tmp"),
                           picard_dir = opts_flow$get("picard_dir"),
                           picard_reorder_opts = opts_flow$get("picard_reorder_opts"),
                           ref_fasta = opts_flow$get('ref_fasta'),
                           samtools_exe = opts_flow$get("samtools_exe")){
  
  if(missing(outfile))
    outfile = gsub(".bam$", "_reorder.bam", x)
  
  check_args()
  
  reorder = sprintf("%s %s -Djava.io.tmpdir=%s -jar %s/picard.jar ReorderSam INPUT=%s OUTPUT=%s REFERENCE=%s;%s index %s",
                    java_exe, java_mem, java_tmp, picard_dir, x, outfile, ref_fasta, samtools_exe, outfile)
  
  ## --- INPUT is a NAMED list
  flowmat = to_flowmat(list(reorder = reorder), samplename)
  return(list(flowmat = flowmat, outfiles = outfile))
}


