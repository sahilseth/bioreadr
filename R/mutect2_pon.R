



# step1: call all variants ------
mutect2_pon_call_variants <- function(bam, 
                                      samplename, 
                                      ref_fasta = opts_flow$get("ref_fasta"), 
                                      gatk4_exe = opts_flow$get("gatk4_exe"), 
                                      mutect.cpu = opts_flow$get("mutect2_cpu"), 
                                      mutect.mem = opts_flow$get("mutect2_mem"),
                                      # mutect.java_opts = opts_flow$get("mutect2.java_opts"), 
                                      mutect2.pon.opts = opts_flow$get("mutect2.pon.opts"),
                                      java_mem_str = opts_flow$get("java_mem_str"),
                                      germline_vcf = opts_flow$get('germline_variants_vcf'),
                                      split_by = "interval_split",
                                      interval_split = NULL
                                      
                                      ){
  
  check_args()
  bamset  = bam_set(bam = bam, 
                    outprefix = samplename,
                    interval_split = interval_split,
                    ref_fasta = ref_fasta, 
                    split_by = split_by)
  
  #mutect_opts = "--artifact_detection_mode -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS -L $regions_bed_fl"
  mutect_vcfs = paste0(bamset$outprefix_interval, ".mutect.vcf.gz")
  intervals = bamset$intervals
  
  cmd_sampname = glue("seqname=$(gatk GetSampleName -I {bam} -O /dev/stdout 2> /dev/null)") %>% as.character()
  # https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_walkers_mutect_CreateSomaticPanelOfNormals.php
  cmd_mutect <- glue("{cmd_sampname}; ",
                     "{gatk4_exe} --java-options '-Xmx{mutect.mem}g -XX:+UseParallelGC -XX:ParallelGCThreads={mutect.cpu}' ",
                     "Mutect2 -R {ref_fasta} -I {bam} ",
                     "--germline-resource {germline_vcf} ",
                     "-tumor $seqname ",
                     "-O {mutect_vcfs} {mutect2.pon.opts} {intervals}") %>% as.character()
  # if(length(cmd_mutect))
  #   stop("some variables are null, cmd_mutect not formed")
  
  # ** merge VCFs ---------
  mutect_vcf = paste0(samplename, ".vcf.gz")
  mutect_vcfs_i = paste0(mutect_vcfs, collapse = " -I ")
  cmd_mergevcf = glue("{gatk4_exe} --java-options {java_mem_str} MergeVcfs -I {mutect_vcfs_i} -O {mutect_vcf};", 
                      "bash ~/Dropbox/public/flow-r/my.ultraseq/pipelines/bin/m2_merge_vcf_stats.sh {mutect_vcf}.stats") %>% as.character()
  
  flowmat = to_flowmat(list(
    mutect = cmd_mutect,
    merge_vcf = cmd_mergevcf
  ), samplename = samplename)
  
  list(flowmat = flowmat, outfiles = list(mutect_vcf = mutect_vcf))
}


mutect2_pon_create <- function(vcfs, oufile){
  
  # create a string of vcfs
  vcfs_str = paste0("-vcfs ", vcfs)
  cmd = glue("gatk CreateSomaticPanelOfNormals \
  -vcfs 3_HG00190.vcf.gz \
  -vcfs 4_NA19771.vcf.gz \
  -vcfs 5_HG02759.vcf.gz \
  -O 6_threesamplepon.vcf.gz")
  
  
}


# example -------
if(FALSE){
  p_load(tidyr)
  out = lapply(1:nrow(df_trk), function(i){
    mutect2_pon_call_variants(bam = df_trk$bam[i], 
                              samplename = df_trk$samplename[i], 
                              ref_fasta = ref_fasta, 
                              gatk4_exe = gatk4_exe, 
                              mutect.mem = mutect.mem, 
                              mutect.cpu = mutect.cpu, 
                              mutect_opts = "", 
                              java_mem = java_mem)
  })
  names(out) = df_trk$samplename
  
  # out = transform(out)
  flowmat = lapply(out, "[[", "flowmat") %>% bind_rows()
  
  write_tsv(flowmat, "flowmat.tsv")
}
