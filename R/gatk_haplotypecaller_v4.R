# source("~/Dropbox/public/flowr/my.ultraseq/my.ultraseq/R/bam_set.R")


#' a wrapper around gatk haplotyper
#' @param x input bam
#' @param haplotyper_opts all additional arguments supported by GATK
#' @details For more details refer to: https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_CatVariants.php
#' 
#' @export
#' 
haplotype_caller <- function(trk, 
                            split_by_chr = TRUE,
                            
                            ref_fasta = opts_flow$get('ref_fasta'),
                            
                            java_exe = opts_flow$get("java_exe"),
                            #java_opts = opts_flow$get("java_opts"),
                            java_mem = opts_flow$get('java_mem_str'),
                            
                            gatk4_exe = opts_flow$get('gatk4_exe'),
                            
                            gatk_hc_cpu = opts_flow$get('gatk_hc_cpu'),
                            # gatk_hc_mem = opts_flow$get("gatk_hc_mem"),
                            gatk_hc_dir = pipe_str$gatk_hc$dir,
                            gatk_hc_opts = opts_flow$get("gatk_hc_opts"),
                            
                            germline_vcf = opts_flow$get('germline_variants_vcf'),
                            
                            # we can use these generic intervals for all kinds of files and 
                            # filter later
                            interval_split = opts_flow$get("capture_bi_wex_booster_intervals_split")

                            ) {

  # no args should be null
  check_args()

  pipename = match.call()[[1]]
  flog.debug(paste0("Generating a haplotypeCaller flowmat for sample: "))

  # gatk_hc_dir = "gatk_hc"
  # gatk_hc_opts = "-ERC GVCF -G StandardAnnotation -G AS_StandardAnnotation -G StandardHCAnnotation"

  # https://gatk.broadinstitute.org/hc/en-us/articles/360037593911-CombineGVCFs
  trk %<>% 
    mutate(gatk_hc_outprefix = file.path(gatk_hc_dir, outprefix), 
           gatk_hc_vcf = paste0(gatk_hc_outprefix, ".hc.g.vcf.gz"),
           gatk_hc_geno_vcf = paste0(gatk_hc_outprefix, ".hc.geno.vcf.gz")) %>% 
    metadata_for_dnaseq()

  # for now this will ONLY work for single tumor sample
  # this maybe needs to be a loop!
  i=1
  tmp = lapply(1:nrow(trk), function(i){
    samplename = trk$samplename[i]
    bam = trk$bam[i]
    flog.debug(glue("working on {i} {samplename} {bam}"))
    # tumor_name = trk_tum$name[i]
    # normal_name = trk_norm$name[1]
    outprefix = trk$gatk_hc_outprefix[i]
    gatk_gvcf = trk$gatk_hc_vcf[i]
    gatk_geno_vcf = trk$gatk_hc_geno_vcf[i]
    
    bamset = bam_set(bam = bam, 
                    outprefix = outprefix,
                    ref_fasta = ref_fasta, 
                    interval_split = interval_split,
                    split_by = "interval_split")

    gatk_gvcfs = paste0(bamset$outprefix_interval, ".vcf.gz")
    # -bamout {mutect_bams}
    gatk_intervals = bamset$intervals
    cmd_gatk_hc = glue("mkdir -p {gatk_hc_dir};module load singularity/3.5.2;{gatk4_exe} --java-options '-Xmx32g' HaplotypeCaller ",
      "-I {bam} -O {gatk_gvcfs} --dbsnp {germline_vcf} ",
      "-R {ref_fasta} --native-pair-hmm-threads 4 {gatk_hc_opts} {gatk_intervals}") 

    gatk_gvcfs_i = paste0(gatk_gvcfs, collapse = " -I ")
    cmd_mergevcf = glue("module load singularity/3.5.2;{gatk4_exe} --java-options {java_mem} GatherVcfs ",
      "-I {gatk_gvcfs_i} -O {gatk_gvcf}; ", 
      "{gatk4_exe} IndexFeatureFile -I {gatk_gvcf}")
    
    # can be used for peddy, and some other QC stuff
    cmd_geno = glue("module load singularity/3.5.2;{gatk4_exe} --java-options {java_mem} GenotypeGVCFs ",
                    "-V {gatk_gvcf} -R {ref_fasta} -O {gatk_geno_vcf}")

    cmds <- list(gatk_hc.splt = cmd_gatk_hc, gatk_hc.mrg = c(cmd_mergevcf, cmd_geno))
    # mutect_gather_bams = cmd_mutect_gather_bams

    flowmat = to_flowmat(cmds, samplename = samplename) %>%
      dplyr::mutate(cmd = as.character(cmd))

    ret = list(flowmat = flowmat)
    ret
  })
  lst = purrr::transpose(tmp)
  lst$flowmat %<>% bind_rows()
  lst$trk = trk

  return(lst)
}



