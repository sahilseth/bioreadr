
#java -jar GenomeAnalysisTK.jar \
#-T DepthOfCoverage \
#-R reference.fasta \
#-o file_name_base \
#-I input_bams.list


#' a wrapper around gatk haplotyper
#' @param x input bam
#' @details For more details refer to: https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_CatVariants.php
#' 
#' @export
#' 
depthofcoverage <- function(x,
                       samplename = opts_flow$get("samplename"),
                       outfile,
                       split_by_chr = TRUE,

                       ref_fasta = opts_flow$get("ref_fasta"),
                       java_exe = opts_flow$get("java_exe"),
                       java_mem = opts_flow$get("java_mem"),
                       java_tmp = opts_flow$get("java_tmp"),
                       gatk_jar = opts_flow$get("gatk_jar"),
                       haplotyper_opts = opts_flow$get("haplotyper_opts"),
                       cpu_haplotyper = opts_flow$get("cpu_haplotyper")) {

  ## no args should be null
  check_args()

  pipename = match.call()[[1]]
  message("Generating a ", pipename, " flowmat for sample: ", samplename)

  bam_prefix <- gsub(".bam", "", basename(x))
  if(split_by_chr){
    chrs_info <- get_fasta_chrs(ref_fasta)
    chrs_prefix <- paste(bam_prefix, chrs_info, sep = "_")
    intervals_opts = paste0(" -L ", chrs_info)             ## interval files

  }else{
    chrs_prefix = paste0(bam_prefix, ".")
    intervals_opts = ""

  }

  outvcf <- paste0(chrs_prefix, ".haplotyper.vcf")

  cmd_haplotyper <- sprintf("%s %s -Djava.io.tmpdir=%s -jar %s -T HaplotypeCaller -R %s -I %s -o %s -nct %s %s %s",
                        java_exe, java_mem, java_tmp, gatk_jar,
                        ref_fasta, x, outvcf, cpu_haplotyper, haplotyper_opts, intervals_opts)

  if(missing(outfile))
    outfile = paste0(bam_prefix, "_merged_haplotyper.vcf")

  cmd_merge = sprintf("%s %s -Djava.io.tmpdir=%s -cp  %s org.broadinstitute.gatk.tools.CatVariants -R %s -out %s -assumeSorted %s",
                      java_exe, java_mem, java_tmp, gatk_jar, ref_fasta, outfile,
                      paste(" -V", outvcf, collapse = ""))

  cmds <- list(haplotyper = cmd_haplotyper, merge_haplotyper = cmd_merge)

  flowmat = to_flowmat(cmds, samplename = samplename)
  return(list(flowmat=flowmat, outfiles=list(all = outvcf, merged = outfile)))
}
