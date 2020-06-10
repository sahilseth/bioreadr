# OVERVIEW ------
# GATK version: gatk_/4.1.1.0 samtools_/1.9R/3.5.2 
# tutorial:
# https://software.broadinstitute.org/gatk/documentation/article?id=11682
# https://software.broadinstitute.org/gatk/documentation/article?id=11683

# TOOLS:
# https://software.broadinstitute.org/gatk/documentation/tooldocs/4.1.1.0/org_broadinstitute_hellbender_tools_copynumber_ModelSegments.php

# ** IGV view -----
# https://gatkforums.broadinstitute.org/gatk/discussion/11224/visualize-cnv-results-from-wgs-analysis

# ** default params borrowed from detitin -----
# https://github.com/broadinstitute/deTiN
# --ascna_probe_number_filter (default = 200) We require 200 probes based on empirical results using GATK4CNV that segments smaller than this tend to be enriched for artifacts. For WGS this parameter can be set to 0.
# --ascna_SNP_number_filter (default = 20) We require 20 SNPs based on empirical results using GATK4CNV that segments smaller than this tend to be enriched for artifacts.

# prepare intervals ------
# module load gatk/4.1.0.0
# bed="/rsrch3/home/iacs/sseth/ref/az_ref_beds/hg19/bed/Exome-NGv3.bed"
# bed2="/rsrch3/home/iacs/sseth/ref/az_ref_beds/hg19/bed/Exome-NGv3_nochr.bed"
# #sed 's/^chr\|%$//g' $bed > $bed2
# gatk_intervals="/rsrch3/home/iacs/sseth/ref/az_ref_beds/hg19/bed/Exome-NGv3_gatk.intervals"
# ref_fasta="/rsrch3/home/iacs/sseth/ref/human/b37/fastas/Homo_sapiens_assembly19.fasta"
# ref_dict="/rsrch3/home/iacs/sseth/ref/human/b37/fastas/Homo_sapiens_assembly19.dict"
# 
# gatk PreprocessIntervals \
# -L $bed2 \
# -R ${ref_fasta}\
# --bin-length 0 \
# --interval-merging-rule OVERLAPPING_ONLY \
# -O $gatk_intervals



# ** het sites gnomad ------
# more details: gatkcnv_example.Rmd
# https://gatkforums.broadinstitute.org/dsde/discussion/11683#ref9
# https://gatkforums.broadinstitute.org/gatk/discussion/11683/how-to-part-ii-sensitively-detect-copy-ratio-alterations-and-allelic-segments


# vcf_af="af-only-gnomad.raw.sites.vcf.gz"
# vcf_af_norm="af-only-gnomad.raw.sites.vcf_norm.gz"
# vcf_af_norm_rehead="af-only-gnomad.raw.sites.vcf_norm_rehead.gz"
# vcf_af_biallelic="af-only-gnomad.raw.sites_biallelic.vcf.gz"
# ref_fasta="Homo_sapiens_assembly19.fasta"
# minimum_allele_frequency=0.001
## Need to have an index because SelectVariants requires it
# gatk IndexFeatureFile -F $vcf_af
# 
# gatk SelectVariants \
# -V  $vcf_af \
# -O $vcf_af_biallelic \
# -R $ref_fasta \
# -select-type SNP -restrict-alleles-to BIALLELIC \
# -select "AF > ${minimum_allele_frequency}"

# source('~/Dropbox/public/flow-r/my.ultraseq/my.ultraseq/R/run_cmd.R')

.collect_read_counts <- function(bam, 
                                 counts,
                                 gatkcnv_intervals = opts_flow$get("gatkcnv_intervals"),
                                 gatk4_exe = opts_flow$get("gatk4_exe"),
                                 run_cmds = F,
                                 redo = F){
  check_args()
  cmd = glue("mkdir -p `dirname {counts}`;{gatk4_exe} --java-options '-Xmx32g' CollectReadCounts -I {bam} -L {gatkcnv_intervals} -O {counts} --format TSV --interval-merging-rule OVERLAPPING_ONLY")
  
  if(run_cmds){
    # else, actually run the cmd
    run_cmd(cmd, target = counts, cmdname = "collectreadcounts", redo = redo)
  }
  
  list(cmds = cmd, outfiles = counts)
}


.denoise_counts <- function(counts,
                            std_cr,
                            dn_cr,
                            gatkcnv_pon_hdfs = opts_flow$get("gatkcnv_pon_hdfs"),
                            gatk4_exe = opts_flow$get("gatk4_exe"),
                            run_cmds = F, redo = F){
  check_args()
  cmd = glue("{gatk4_exe} --java-options '-Xmx12g' DenoiseReadCounts ", 
             "-I {counts} --count-panel-of-normals {gatkcnv_pon_hdfs} ",
             "--standardized-copy-ratios {std_cr} ",
             "--denoised-copy-ratios {dn_cr}")
  
  # if(test) return(cmd)
  # else, actually run the cmd
  if(run_cmds){
    run_cmd(cmd, target = dn_cr, cmdname = "denoise_counts", redo = redo)
  }
  
  list(cmds = cmd, outfiles = dn_cr)
  
}

.plot_denoise_cr <- function(std_cr,
                             dn_cr,
                             outprefix,
                             ref_dict = opts_flow$get("ref_dict"),
                             gatk4_exe = opts_flow$get("gatk4_exe"),
                             java_tmp = opts_flow$get("java_tmp"),
                             plot_minimum_contig_length = "46709983",
                             run_cmds = F, redo = F){
  
  check_args();cmdname = "plot_denoise_cr"
  cmd = glue("{gatk4_exe} --java-options '-Xmx12g -Djava.io.tmpdir={java_tmp}' PlotDenoisedCopyRatios ", 
             "--standardized-copy-ratios {std_cr} ",
             "--denoised-copy-ratios {dn_cr} ",
             "--sequence-dictionary {ref_dict} ",
             "--minimum-contig-length 46709983 ",
             "--output . ",
             "--output-prefix {outprefix}")
  
  # if(test) return(cmd)
  # else, actually run the cmd
  outplot = glue("{outprefix}.denoisedLimit4.png")
  if(run_cmds){
    run_cmd(cmd, target = outplot, cmdname = cmdname, 
            stderr = glue("{outprefix}_{cmdname}.log"), redo = redo)
  }
  
  list(cmds = cmd, outfiles = outplot)
  
}



.collect_allelic_counts <- function(bam,
                                    acounts,
                                    gatk4_exe = opts_flow$get("gatk4_exe"),
                                    ref_fasta = opts_flow$get("ref_fasta"),
                                    germline_variants_biallelic_vcf = opts_flow$get("germline_variants_biallelic_vcf"),
                                    
                                    run_cmds = T, redo = F){
  check_args()
  cmd = glue("mkdir -p `dirname {acounts}`;{gatk4_exe} --java-options '-Xmx32g' CollectAllelicCounts ", 
             "-L {germline_variants_biallelic_vcf} ", 
             "-I {bam} ", 
             "-R {ref_fasta} ", 
             "-O {acounts}")
  
  # if(test) return(cmd)
  if(run_cmds){
    # else, actually run the cmd
    run_cmd(cmd, target = acounts, 
            cmdname = "collect_allelic_counts", 
            stderr = glue("{outprefix}_collect_allelic_counts.log"), redo = redo)
  }
  list(cmds = cmd, outfiles = acounts)
}

.model_segments <- function(bam,
                            dn_cr,
                            outprefix, # usually {outprefix}_clean_wimn OR {outprefix}_clean_womn (with or without matched normal)
                            acounts,
                            normal_acounts = NULL,
                            
                            # output of ModelSegments
                            hets,
                            cr_seg,
                            final_seg,
                            
                            # output of CallCopyRatioSegments
                            called_seg,
                            
                            gatk4_exe = opts_flow$get("gatk4_exe"),
                            ref_dict = opts_flow$get("ref_dict"),
                            # ref_fasta = opts_flow$get("ref_fasta"),
                            plot_minimum_contig_length = "46709983",
                            
                            gatkcnv_modelseg_params = "--minimum-total-allele-count 30 --maximum-number-of-smoothing-iterations 100 ",
                            redo = F, 
                            test = F, run_cmds = F
){
  
  check_args(ignore = "normal_acounts")
  cmd1 = glue("{gatk4_exe} --java-options '-Xmx16g' ModelSegments ", 
              "--allelic-counts {acounts} ",
              "--denoised-copy-ratios {dn_cr} ",
              # 30 is default: alleleic copy ratios
              # default 25 for WGS
              "{gatkcnv_modelseg_params} ",
              "--output . --output-prefix {outprefix}")
  if(!is.null(normal_acounts)){
    # message("adding matched normal counts to ModelSegments")
    cmd1 = glue("{cmd1} --normal-allelic-counts {normal_acounts}")
  }
  
  # not sure about the outputs of this one
  cmd2 = glue("{gatk4_exe} --java-options '-Xmx16g' CallCopyRatioSegments ", 
              "--input {cr_seg} --output {called_seg}")
  
  cmd3 = glue("{gatk4_exe} --java-options '-Xmx16g' PlotModeledSegments ", 
              "--denoised-copy-ratios {dn_cr} ", 
              "--allelic-counts {hets} ", 
              "--segments {final_seg} ", 
              "--sequence-dictionary {ref_dict} ", 
              "--minimum-contig-length {46709983} ", 
              "--output . ", 
              "--output-prefix {outprefix}")
  
  
  # else, actually run the cmd
  if(run_cmds){
    run_cmd(cmd1, target = final_seg, cmdname = "modelsegments", 
            stderr = glue("{outprefix}_modelsegments.log"), redo = redo)
    run_cmd(cmd2, target = called_seg, cmdname = "CallCopyRatioSegments", 
            stderr = glue("{outprefix}_callcopyratiosegments.log"), redo = redo)
    # 185_145_GB-D_clean_plotmodeledsegments.log
    run_cmd(cmd3, target = glue("{outprefix}.modeled.png"), cmdname = "PlotModeledSegments", 
            stderr = glue("{outprefix}_plotmodeledsegments.log"), redo = redo)
  }
  
  list(cmds = c(cmd1, cmd2, cmd3), outfiles = c(called_seg, glue("{outprefix}.modeled.png")))
  
}


gatkcnv_somatic_pon <- function(df_trk,
                                num_cores = nrow(df_trk),
                                run_cmds = F){
  
  # ceate a new trk with ALL reqd files
  df_trk %<>% 
    mutate(
      bam = basename(bam),
      oprefix = glue("gatkcnv/{name}"),
      gatkcnv_counts = glue("{oprefix}.counts.tsv"), 
      gatkcnv_std_cr = glue("{oprefix}.standardizedCR.tsv"),
      gatkcnv_dn_cr = glue("{oprefix}.denoisedCR.tsv"),
      gatkcnv_acounts = glue("{oprefix}.allelicCounts.tsv"),
      
      # will always run WOMN
      ## ModelSegments outputs
      gatk_womn_oprefix = glue("gatkcnv/{name}_clean_womn"),
      gatkcnv_womn_hets = glue("{gatk_womn_oprefix}.hets.tsv"),
      gatkcnv_womn_cr_seg = glue("{gatk_womn_oprefix}.cr.seg"),
      gatkcnv_womn_final_seg = glue("{gatk_womn_oprefix}.modelFinal.seg"),
      ## CallCopyRatioSegments outputs
      gatkcnv_womn_called_seg = glue("{gatk_womn_oprefix}.called.seg")
      
      
    )
  df_trk
  
  wranglr::mkdir("gatkcnv")
  
  
  # ** .collectreadcounts -------
  # need for ALL samples  i=1
  tmp = mclapply(1:nrow(df_trk), function(i){
    .collect_read_counts(df_trk$BAM[i], 
                         df_trk$gatkcnv_counts[i],
                         run_cmds=run_cmds, redo = F)
  }, mc.cores = num_cores)
  # tmp
  warnings()
  
  # ** denoise -------
  tmp = mclapply(1:nrow(df_trk), function(i){
    .denoise_counts(bam = df_trk$BAM[i], 
                    counts = df_trk$gatkcnv_counts[i],
                    std_cr = df_trk$gatkcnv_std_cr[i],
                    dn_cr = df_trk$gatkcnv_dn_cr[i],
                    run_cmds=run_cmds)
  }, mc.cores = num_cores)
  # tmp
  warnings()
  
  # ** plot denoise -------
  tmp = mclapply(1:nrow(df_trk), function(i){
    .plot_denoise_cr(std_cr = df_trk$gatkcnv_std_cr[i],
                     dn_cr = df_trk$gatkcnv_dn_cr[i],
                     outprefix = df_trk$oprefix[i], 
                     run_cmds=run_cmds, redo = F)
  }, mc.cores = num_cores)
  # tmp
  warnings()
  
  # ** .collectallelicounts -------
  i=1
  tmp = mclapply(1:nrow(df_trk), function(i){
    .collect_allelic_counts(bam = df_trk$BAM[i], 
                            acounts = df_trk$gatkcnv_acounts[i],
                            run_cmds=run_cmds)
  }, mc.cores = 1)
  # tmp
  warnings()
  
  # ** modelSegments -------
  # womn
  i=3
  tmp = mclapply(1:nrow(df_trk), function(i){
    .model_segments(bam = df_trk$BAM[i], 
                    dn_cr = df_trk$gatkcnv_dn_cr[i],
                    outprefix = df_trk$gatk_womn_oprefix[i],
                    acounts = df_trk$gatkcnv_acounts[i],
                    hets = df_trk$gatkcnv_womn_hets[i],
                    cr_seg = df_trk$gatkcnv_womn_cr_seg[i],
                    final_seg = df_trk$gatkcnv_womn_final_seg[i],
                    gatkcnv_modelseg_params = "--minimum-total-allele-count-case 30 --maximum-number-of-smoothing-iterations 25 --number-of-smoothing-iterations-per-fit 1 --kernel-variance-allele-fraction 0.8 --kernel-variance-copy-ratio 0.8",
                    called_seg = df_trk$gatkcnv_womn_called_seg[i],
                    run_cmds=run_cmds, redo = F)
  }, mc.cores = num_cores)
  # tmp
  
  # Number of segmentation-smoothing iterations per MCMC model refit. 
  # (Increasing this will decrease runtime, but the final number of segments may be higher. 
  # Setting this to 0 will completely disable model refitting between iterations.)
  
  # --minimum-total-allele-count-case:Integer
  # --minimum-total-allele-count-case 30: min depth of seq in case
  # (but then how would you call COMPLETE LOSS?)
  # Minimum total count for filtering allelic counts in the case sample, if available.  The
  # default value of zero is appropriate for matched-normal mode; increase to an appropriate
  # value for case-only mode.  Default value: 0. 
  
  write_tsv(df_trk, path = "gatkcnv/df_trk.tsv")
}


gatkcnv_somatic_matched <- function(trk,
                                    samplename,
                                    # outpath = "gatkcnv/",
                                    num_cores = nrow(df_trk), 
                                    run_cmds = F){
  # check generic columns
  trk = metadata_for_dnaseq_tools(trk)
  expect_columns(trk, c("outprefix", "outprefix_paired"))
  
  trk = mutate(outprefix = file.path(pipe_str$gtkcnv$dir, outprefix), 
         outprefix_paired = file.path(pipe_str$gtkcnv$dir, outprefix_paired))
  
  # ceate a new trk with ALL reqd files
  trk %<>% 
    mutate(
      # skip taking basename of bam!
      # bam = basename(bam),
      # oprefix = glue("{outpath}{name}"),
      gatkcnv_counts = glue("{outprefix}.counts.tsv"), 
      gatkcnv_std_cr = glue("{outprefix}.standardizedCR.tsv"),
      gatkcnv_dn_cr = glue("{outprefix}.denoisedCR.tsv"),
      gatkcnv_acounts = glue("{outprefix}.allelicCounts.tsv"),
      
      ## ModelSegments outputs
      # oprefix_final = glue("{outpath}{name}_clean"),
      gatkcnv_hets = glue("{outprefix_paired}.hets.tsv"),
      gatkcnv_cr_seg = glue("{outprefix_paired}.cr.seg"),
      gatkcnv_final_seg = glue("{outprefix_paired}.modelFinal.seg"),
      ## CallCopyRatioSegments outputs
      gatkcnv_called_seg = glue("{outprefix_paired}.called.seg")
    )
  trk
  
  trk_tum = filter(trk, normal == "NO")
  trk_norm = filter(trk, normal == "YES")
  
  # wranglr::mkdir("gatkcnv")
  
  
  # ** .collectreadcounts -------
  # need for ALL samples  i=1
  out_readcnts = mclapply(1:nrow(trk), function(i){
    .collect_read_counts(trk$bam[i], 
                         trk$gatkcnv_counts[i],
                         run_cmds=run_cmds, redo = F)
  }, mc.cores = num_cores)
  # tmp
  warnings()
  
  # ** denoise -------
  out_dr = mclapply(1:nrow(trk), function(i){
    .denoise_counts(counts = trk$gatkcnv_counts[i],
                    std_cr = trk$gatkcnv_std_cr[i],
                    dn_cr = trk$gatkcnv_dn_cr[i],
                    run_cmds=run_cmds)
  }, mc.cores = num_cores)
  # tmp
  warnings()
  
  # ** plot denoise -------
  out_plot_dr = mclapply(1:nrow(trk), function(i){
    .plot_denoise_cr(std_cr = trk$gatkcnv_std_cr[i],
                     dn_cr = trk$gatkcnv_dn_cr[i],
                     outprefix = trk$outprefix[i], 
                     run_cmds=run_cmds, redo = F)
  }, mc.cores = num_cores)
  # tmp
  warnings()
  
  # ** .collectallelicounts -------
  i=1
  out_acounts = mclapply(1:nrow(trk), function(i){
    .collect_allelic_counts(bam = trk$bam[i], 
                            acounts = trk$gatkcnv_acounts[i],
                            run_cmds=run_cmds)
  }, mc.cores = 1)
  # tmp
  warnings()
  
  # ** modelSegments -------
  # womn
  i=3
  out_seg = mclapply(1:nrow(trk_tum), function(i){
    .model_segments(bam = trk_tum$bam[i], 
                    dn_cr = trk_tum$gatkcnv_dn_cr[i],
                    outprefix = trk_tum$outprefix_paired[i],
                    acounts = trk_tum$gatkcnv_acounts[i],
                    normal_acounts = trk_norm$gatkcnv_acounts[1],
                    
                    hets = trk_tum$gatkcnv_hets[i],
                    cr_seg = trk_tum$gatkcnv_cr_seg[i],
                    final_seg = trk_tum$gatkcnv_final_seg[i],
                    gatkcnv_modelseg_params = "--minimum-total-allele-count-case 30 --maximum-number-of-smoothing-iterations 25 --number-of-smoothing-iterations-per-fit 1 --kernel-variance-allele-fraction 0.8 --kernel-variance-copy-ratio 0.8",
                    called_seg = trk_tum$gatkcnv_called_seg[i],
                    run_cmds=run_cmds, redo = F)
  }, mc.cores = num_cores)
  # tmp
  
  # Number of segmentation-smoothing iterations per MCMC model refit. 
  # (Increasing this will decrease runtime, but the final number of segments may be higher. 
  # Setting this to 0 will completely disable model refitting between iterations.)
  
  # --minimum-total-allele-count-case:Integer
  # --minimum-total-allele-count-case 30: min depth of seq in case
  # (but then how would you call COMPLETE LOSS?)
  # Minimum total count for filtering allelic counts in the case sample, if available.  The
  # default value of zero is appropriate for matched-normal mode; increase to an appropriate
  # value for case-only mode.  Default value: 0. 
  
  # write_tsv(df_trk, path = "gatkcnv/df_trk.tsv")
  cmds <- list(gtkcnv.splt = unlist(out_readcnts), 
               gtkcnv.splt = unlist(out_acounts),
               gtkcnv.mrg = unlist(out_dr),
               gtkcnv.mrg = unlist(out_plot_dr),
               gtkcnv.mrg = unlist(out_seg))
  # mutect_gather_bams = cmd_mutect_gather_bams
  
  flowmat = to_flowmat(cmds, samplename = samplename) %>% 
    mutate(cmd = as.character(cmd))
  
  list(flowmat = flowmat, outfiles = list(), trk = trk)
}



# somatic: trial/example/debug --------
if(FALSE){
  # module load R_/3.5.3
  
  # R CMD INSTALL ~/Dropbox/public/flow-r/my.ultraseq/my.ultraseq
  # module load R/3.5.2 gatk_/4.1.1.0 samtools_/1.9
  
  library(pacman)
  p_load(dplyr, readr, janitor, glue, magrittr)
  p_load(flowr, my.ultraseq, parallel)
  
  outpath="/rsrch3/scratch/iacs/sseth/flows/SS/tnbc/ms51_wex_b1/ss_cl_het/185_145"
  setwd(outpath)
  
  df_trk_full = read_tsv("~/projects2/ss_tnbc/data/artemis/wex/2019_b1/superfreq_v2/df_sfq_metadata_full.tsv")
  df_trk_ms51 = filter(df_trk_full, subproject == "ms51_dna_b1")
  df_trk_ms51$INDIVIDUAL %>% unique() %>% sort()
  df_trk = filter(df_trk_ms51, INDIVIDUAL == "185_145")
  write_tsv(df_trk, glue("{outpath}/df_trk.tsv"))
  
  # ** load vars --------
  opts_flow$load("~/Dropbox/public/flow-r/my.ultraseq/pipelines/dnaseq/ss_cl_het/somatic_snv_indels.conf")
  opts_flow$set(
    gatk4_exe = "$HOME/apps/gatk/gatk-4.1.1.0/gatk",
    gatk_intervals = "$HOME/ref/az_ref_beds/hg19/bed/Exome-NGv3_gatk.intervals",
    # pon_counts = "$HOME/ref/human/b37/common_normal/ms51_v1/gatkcnv_pon/v1_sfq/cnvponC.pon.hdf5",
    pon_counts = "$HOME/ref/human/b37/common_normal/ms51_v1/gatkcnv_pon/cnvpon_ms51_tier1_tier2.pon.hdf5",
    vcf_af_biallelic="$HOME/ref/human/b37/annotations/gnomad/af-only-gnomad.raw.sites_biallelic.vcf.gz"  
  )
  opts_flow$get("ref_fasta")
  
  # load conf (instead of vars)
  
  
  # ** fetch fls --------
  fls = df_trk$BAM
  fetch_files.exe = "~/bin/fetch_files_v2.sh"
  cmd_transfer = glue("{fetch_files.exe} {fls}*") %>% as.character()
  message("> transferring files")
  # tmp = mclapply(cmd_transfer, system, mc.cores = 4)
  
  # ** full pipe ------
  source('~/Dropbox/public/flowr/my.ultraseq/my.ultraseq/R/run_cmd.R')
  source('~/Dropbox/public/flowr/my.ultraseq/pipelines/dnaseq/ss_cl_het/gatkcnv.R')
  gatkcnv_somatic(df_trk, num_cores = 1)
  
  
  
  
}



# https://gatkforums.broadinstitute.org/gatk/discussion/11684
gatkcnv_germline <- function(){
  
  # AnnotateIntervals (once) -------
  # gatk AnnotateIntervals \
  # -L chr20XY.interval_list \
  # -R ref/Homo_sapiens_assembly38.fasta \
  # -imr OVERLAPPING_ONLY \
  # -O chr20XY.annotated.tsv
  # Although optional for the tool, we recommend annotating mappability by providing a 
  # --mappability-track regions file in either .bed or .bed.gz format. Be sure to merge any overlapping intervals beforehand. The tutorial omits use of this resource.
  # https://bismap.hoffmanlab.org/
  
  # filter intervals -------
  # gatk FilterIntervals \
  # -L chr20XY.interval_list \
  # --annotated-intervals chr20XY.annotated.tsv \
  # -I cvg/HG00096.tsv -I cvg/HG00268.tsv -I cvg/HG00419.tsv -I cvg/HG00759.tsv \
  # -I cvg/HG01051.tsv -I cvg/HG01112.tsv -I cvg/HG01500.tsv -I cvg/HG01565.tsv \
  # -I cvg/HG01583.tsv -I cvg/HG01595.tsv -I cvg/HG01879.tsv -I cvg/HG02568.tsv \
  # -I cvg/HG02922.tsv -I cvg/HG03006.tsv -I cvg/HG03052.tsv -I cvg/HG03642.tsv \
  # -I cvg/HG03742.tsv -I cvg/NA18525.tsv -I cvg/NA18939.tsv -I cvg/NA19017.tsv \
  # -I cvg/NA19625.tsv -I cvg/NA19648.tsv -I cvg/NA20502.tsv -I cvg/NA20845.tsv \
  # -imr OVERLAPPING_ONLY \
  # -O chr20XY.cohort.gc.filtered.interval_list
  # imp default thresholds  
  # --low-count-filter-count-threshold default is 5
  # --low-count-filter-percentage-of-samples default is 90.0
  # --extreme-count-filter-minimum-percentile default is 1.0
  # --extreme-count-filter-maximum-percentile default is 99.0
  # --extreme-count-filter-percentage-of-samples default is 90.0
  
  # To disable counts based filtering, omit the read counts or, 
  # e.g. when using the v4.1.0.0 cnv_germline_cohort_workflow.wdl pipeline script, set the two percentage-of-samples parameters as follows.
  # --low-count-filter-percentage-of-samples 100 \
  # --extreme-count-filter-percentage-of-samples 100 \
  
  
  # collectreadcounts -------
  
  
}




"Support for Copy Number Variations (CNVs) with GATK4
https://software.broadinstitute.org/gatk/documentation/article?id=11682
https://gatkforums.broadinstitute.org/dsde/discussion/11683/

https://gatkforums.broadinstitute.org/gatk/discussion/7387/description-and-examples-of-the-steps-in-the-acnv-case-workflow
"

# prepare intervals:

# 1. Collect raw counts data with PreprocessIntervals and CollectFragmentCounts
# Before collecting coverage counts that forms the basis of copy number variant detection, we define the resolution of the analysis with a genomic intervals list. The extent of genomic coverage and the size of genomic intervals in the intervals list factor towards resolution.
# 
# Preparing a genomic intervals list is necessary whether an analysis is on targeted exome data or whole genome data. In the case of exome data, we pad the target regions of the capture kit. In the case of whole genome data, we divide the reference genome into equally sized intervals or bins. In either case, we use PreprocessIntervals to prepare the intervals list.
# 
# For the tutorial exome data, we provide the capture kit target regions in 1-based intervals and set --bin-length to zero.

# gatk PreprocessIntervals \
# -L targets_C.interval_list \
# -R /gatk/ref/Homo_sapiens_assembly38.fasta \
# --bin-length 0 \
# --interval-merging-rule OVERLAPPING_ONLY \
# -O sandbox/targets_C.preprocessed.interval_list




# ![](https://us.v-cdn.net/5019796/uploads/editor/3z/gim58s5j2wmk.png)

# cnv-C: panel of normals for cnv really helps
# this is a critical part here
# used for 

# Step 1. Het Pulldown
# ** These instructions describe one method for Het Pulldown for matched samples. For more options, including tumor-only, please see: http://gatkforums.broadinstitute.org/gatk/discussion/7719/overview-of-getbayesianhetcoverage-for-heterozygous-snp-calling **
#   
#   Inputs
# control_bam -- BAM file for control sample (normal).
# case_bam -- BAM file for case sample (tumor).
# reference_sequence -- FASTA file for b37 reference.
# snp_file -- Picard interval list of common SNP sites at which to test for heterozygosity in the control sample .
# Outputs
# normal_het_pulldown -- TSV file with M entries containing ref/alt counts, ref/alt bases, etc., where M is the number of hets called in the control sample.
# tumor_het_pulldown -- TSV file with M entries containing ref/alt counts, ref/alt bases, etc. for sites in the case sample that were called as het in the control sample, where M is the number of hets called in the control sample.
# Format for both output files:
#   
#   CONTIG  POSITION        REF_COUNT       ALT_COUNT       REF_NUCLEOTIDE  ALT_NUCLEOTIDE  READ_DEPTH
# 1       809876  5       16      A       G       21
# 1       881627  23      12      G       A       35
# 1       882033  9       10      G       A       19
# 1       900505  26      24      G       C       50
# ....snip....
# Invocation
# java -jar <path_to_gatk_protected_jar> GetBayesianHetCoverage --reference <reference_sequence>
#   --snpIntervals <snp_file> --tumor <case_bam> --tumorHets <tumor_het_pulldown> --normal <control_bam>
#   --normalHets <normal_het_pulldown> --hetCallingStringency 30

# Step 2. Allelic CNV
# Inputs
# tumor_het_pulldown -- Generated in step 1.
# coverage_profile -- Tangent-normalized coverage TSV file obtained in the GATK CNV case workflow.
# called_segments -- Called-segments TSV file obtained in the GATK CNV case workflow.
# output_prefix -- Path and file prefix for creating the output files. For example, /home/lichtens/my_acnv_output/sample1
# Outputs
# acnv_segments -- TSV file with name ending with -sim-final.seg containing posterior summary statistics for log_2 copy ratio and minor-allele fraction in each segment. Using the above output_prefix, /home/lichtens/my_acnv_output/sample1-sim-final.seg
# acnv_cr_parameters -- TSV file with name ending with -sim-final.cr.param containing posterior summary statistics for global parameters of the copy-ratio model. Using the above output_prefix, /home/lichtens/my_acnv_output/sample1-sim-final.cr.param
# acnv_af_parameters -- TSV file with name ending with -sim-final.af.param containing posterior summary statistics for global parameters of the allele-fraction model. Using the above output_prefix, /home/lichtens/my_acnv_output/sample1-sim-final.af.param
# Other files containing intermediate results of the calculation are also generated.
# 
# Invocation
# java -Xmx8g -jar <path_to_gatk_protected_jar> AllelicCNV  --tumorHets <tumor_het_pulldown>
#   --tangentNormalized <coverage_profile> --segments <called_segments> --outputPrefix <output_prefix>
#   Step 3. Call CNLoH and Balanced Segments
# ** WARNING: This tool is experimental and exists primarily for internal Broad use. **
#   
#   Inputs
# tumor_het_pulldown -- Generated in step 1.
# acnv_segments -- Generated in step 2 (*-sim-final.seg).
# coverage_profile -- Tangent-normalized coverage TSV file obtained in the GATK CNV case workflow
# output_dir -- Directory for creating the output files. For example, /home/lichtens/my_acnv_cnlohcalls_output/
#   Outputs
# GATK-CNV-formatted seg file -- TSV file ending with -sim-final.cnv.seg. This file is formatted identically as the output of GATK CNV. Note that this implies that the allelic fraction values are not captured in this file.
# AllelicCapSeg-formatted seg file -- TSV file ending with -sim-final.acs.seg. This file is formatted identically as the output of Broad CGA AllelicCapSeg. Note that this file can be used as input to Broad-internal versions of ABSOLUTE.
# TITAN-compatible het file --TSV file ending with -sim-final.titan.het.tsv. This file can be used as the input to TITAN for the het read counts.
# TITAN-compatible copy-ratio file -- TSV file ending with -sim-final.titan.tn.tsv. This file can be used as the input to TITAN for the per-target copy-ratio estimates.
# Invocation
# java -Xmx8g -jar <path_to_gatk_protected_jar> CallCNLoHAndSplits  --tumorHets <tumor_het_pulldown>
#   --segments <acnv_segments> --tangentNormalized <coverage_profile> --outputDir <output_dir>
#   --rhoThreshold 0.2 --numIterations 10  --sparkMaster local[*]  

# def _titan_cn_file(cnr_file, work_dir, data):
#   """Convert CNVkit or GATK4 normalized input into TitanCNA ready format.
#     """
# out_file = os.path.join(work_dir, "%s.cn" % (utils.splitext_plus(os.path.basename(cnr_file))[0]))
# support_cols = {"cnvkit": ["chromosome", "start", "end", "log2"],
#   "gatk-cnv": ["CONTIG", "START", "END", "LOG2_COPY_RATIO"]}
# cols = support_cols[cnvkit.bin_approach(data)]
# if not utils.file_uptodate(out_file, cnr_file):
#   with file_transaction(data, out_file) as tx_out_file:
#   iterator = pd.read_table(cnr_file, sep="\t", iterator=True, header=0, comment="@")
# with open(tx_out_file, "w") as handle:
#   for chunk in iterator:
#   chunk = chunk[cols]
# chunk.columns = ["chrom", "start", "end", "logR"]
# if cnvkit.bin_approach(data) == "cnvkit":
#   chunk['start'] += 1
# chunk.to_csv(handle, mode="a", sep="\t", index=False)
# return out_file


# def _titan_het_file(vrn_files, work_dir, paired):
#   assert vrn_files, "Did not find compatible variant calling files for TitanCNA inputs"
# from bcbio.heterogeneity import bubbletree
# class OutWriter:
#   def __init__(self, out_handle):
#   self.writer = csv.writer(out_handle, dialect="excel-tab")
# def write_header(self):
#   self.writer.writerow(["Chr", "Position", "Ref", "RefCount", "Nref", "NrefCount", "NormQuality"])
# def write_row(self, rec, stats):
#   if rec.qual and float(rec.qual) > 0:
#   self.writer.writerow([rec.chrom, rec.pos, rec.ref, stats["tumor"]["depth"] - stats["tumor"]["alt"],
#                         rec.alts[0], stats["tumor"]["alt"], rec.qual])
# return bubbletree.prep_vrn_file(vrn_files[0]["vrn_file"], vrn_files[0]["variantcaller"],
#                                 work_dir, paired, OutWriter)











# END
