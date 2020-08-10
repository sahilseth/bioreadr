



# step1: call all variants ------
mutect2_pon_call_variants <- function(bam, 
                                      samplename, 
                                      ref_fasta = opts_flow$get("ref_fasta"), 
                                      gatk4_exe = opts_flow$get("gatk4_exe"), 
                                      mutect2_cpu = opts_flow$get("mutect2_cpu"), 
                                      mutect2_mem = opts_flow$get("mutect2_mem"),
                                      mutect2_merge_stats_exe = opts_flow$get("mutect2_merge_stats_exe"),
                                      java_mem = opts_flow$get("java_mem_str"),
                                      split_by = "interval_split",
                                      interval_split = opts_flow$get("capture_bi_wex_booster_intervals_split")
                                      
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
                     "{gatk4_exe} --java-options '-Xmx{mutect2_mem}g -XX:+UseParallelGC -XX:ParallelGCThreads={mutect2_cpu}' ",
                     "Mutect2 -R {ref_fasta} -I {bam} ",
                     "-O {mutect_vcfs} --max-mnp-distance 0 {intervals}") %>% as.character()
  # if(length(cmd_mutect))
  #   stop("some variables are null, cmd_mutect not formed")
  
  # ** merge VCFs ---------
  mutect_vcf <- paste0(samplename, ".vcf.gz")
  mutect_vcfs_i <- paste0(mutect_vcfs, collapse = " -I ")
  # GatherVcfs
  # https://software.broadinstitute.org/gatk/documentation/tooldocs/current/picard_vcf_GatherVcfs.php
  # gather should be faster
  cmd_mergevcf <- glue(
    "{gatk4_exe} --java-options {java_mem} GatherVcfs -I {mutect_vcfs_i} -O {mutect_vcf}; ",
    "bash {mutect2_merge_stats_exe} '*vcf.gz.stats' {mutect_vcf}.stats;",
    "{gatk4_exe} IndexFeatureFile -I {mutect_vcf}")
  # we need two arguments, pattern and output file filename for the stats
  
  cmd_rm_tmp_vcf = glue("rm -f {mutect_vcfs}")
  cmd_rm_tmp_vcf = ""

  flowmat = to_flowmat(list(
    mutect2.splt = cmd_mutect,
    mutect2.mrg = c(cmd_mergevcf, cmd_rm_tmp_vcf)
  ), samplename = samplename)
  
  list(flowmat = flowmat, outfiles = list(mutect_vcf = mutect_vcf))
}

# https://gatk.broadinstitute.org/hc/en-us/articles/360042479112-CreateSomaticPanelOfNormals-BETA-
mutect2_pon_create <- function(vcfs,
                               ponm_db = "pon_db",
                               ponm_vcf = "ponm.vcf.gz",
                               ponm_vcf2 = "ponm.vcf.gz",
                               gatk4_exe = opts_flow$get('gatk4_exe'),
                               java_mem = opts_flow$get("java_mem_str"),
                               ref_fasta = opts_flow$get('ref_fasta'),
                               germline_vcf = opts_flow$get('germline_variants_vcf')

                               ){
  
  # create a string of vcfs
  vcfs_v = paste0(vcfs, collapse = " -V ")

  # create ponm_db
  tmp <- get_fasta_chrs(opts_flow$envir$ref_fasta)
  intervals = paste0(" -L ", tmp$chrs[1:25]) %>% paste0(collapse = "")
  # paste0(" -L ", paste0(tmp$chrs[1:25], collapse = ","))
  cmd1 = glue("{gatk4_exe} --java-options {java_mem} GenomicsDBImport -R {ref_fasta} ",
        "--genomicsdb-workspace-path {ponm_db} --overwrite-existing-genomicsdb-workspace ",
        # limit batch size to prevent memory leak
        "-V {vcfs_v} {intervals} --reader-threads 24 --batch-size 200")

  cmd2 = glue("{gatk4_exe} CreateSomaticPanelOfNormals -R {ref_fasta} ",
              "--germline-resource {germline_vcf} ",
              "-V gendb://{ponm_db} ",
              "-O {ponm_vcf} ")

  # older style of doing this ( does not work with latest GATK )
  cmd3 = glue("{gatk4_exe} CreateSomaticPanelOfNormals -R {ref_fasta} ",
              "--germline-resource {germline_vcf} ",
              "-V {vcfs_v} ",
              "-O {ponm_vcf2} ")

  return(c(cmd1, cmd2))
}

# update DB instead
# gatk GenomicsDBImport \
#     --genomicsdb-update-workspace-path existing_database
#     --output-interval-list-to-file /path/to/output/file


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
