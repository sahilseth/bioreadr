#' Title
#'
#' @param bam 
#' @param bed 
#' @param df_bed 
#'
#'
mpileup <- function(bam, 
                    samplename = "samp",
                    ref_fa = "/rsrch2/iacs/iacs_dep/sseth/ref/human/annotations/gencode/v19/GRCh37.p13.genome.fa",
                    regions_bed = "/rsrch2/iacs/iacs_dep/sseth/projects2/ss_pdac_subtyping/20161012_rna_calls/regions.bed",
                    outfile_raw,
                    outfile_call,
                    samtools_exe = "samtools", 
                    bcftools_exe = "bcftools",
                    execute = FALSE){
  # check for index, if not index it
  # needs absolute path
  #bam="Sample_PDX-PDAC-153-FX-R11_vs_human_after_remove_mouse_reads_Aligned.out.WithReadGroup.sorted.bam"
  
  if(missing(outfile_raw)){
    outfile_raw = tempfile()
  }
  
  if(!missing(outfile_call)){
    if(!dir.exists(dirname(outfile_call)))
      dir.create(dirname(outfile_call))
  }
  
  #-t output DP: depth
  #-u: uncompressed
  #-v: vcf
  #-r regions
  # https://github.com/samtools/samtools/issues/79
  # -r does not handle insertion and deletions correctly
  # also
  #DP4 Number of 1) forward ref alleles; 2) reverse ref; 3) forward non-ref; 4) reverse non-ref alleles, used in variant calling.
  #Sum can be smaller than DP because low-quality bases are not counted.
  # -t DP these are already included-t DP,DP4
  cmd_mpileup = sprintf("%s mpileup -u -l %s -f %s %s > %s", 
                        samtools_exe, regions_bed, ref_fa, bam, outfile_raw)
  
  #filter to get only non-ref
  cmd_call = sprintf("%s call -vc %s > %s", 
                     bcftools_exe, outfile_raw, outfile_call)
  
  cmds = list(cmd_mpileup = cmd_mpileup, cmd_call = cmd_call)
  flowmat = to_flowmat(cmds, samplename)
  
  if(execute)
    sapply(flowmat$cmd, system)

  ret = list(outfile_call = outfile_call, flowmat = flowmat)
  ret
}
