
#' Flow following Broad's best practices for variant calling, starting from sorted bam
#'
#' @title Pre-process bam files following Broad's best practices for variant calling, starting from aligned BAM file
#' @description This function provides a wrapper around the best practices described on \href{https://www.broadinstitute.org/gatk/guide/bp_step.php?p=1}{GATK's website}.
#' If the link is broken google 'GATK best practices'
#'
#' This aims to perform the following steps ( for DNA ):
#'
#' \itemize{
#' \item mark duplicates
#' \item realign indels
#' \item recalibrate bases
#' \item current version: \emph{3.4-46}
#' }
#'
#' For RNA GATK recommends a additional step of split n trim, which is not currently supported (contributions welcome !).
#'
#' \strong{NOTE}:
#'
#' Some GATK tools use \href{https://www.broadinstitute.org/gatk/guide/article?id=1975}{cpu threads while others use data threads},
#' flowr tries to use efficiently make the best use of both/either depending on tool's compatibility.
#'
#' @param bam bam file path
#' @param samplename name of the sample
#' @param outfile output file name
#' 
#' @param java_exe path to java
#' @param java_tmp path to java_tmp, can leave blank
#' 
#' @param split_by_chr split processing by chromosomr where ever possible
#' 
#' @param gatk_jar_path path to gatk jar file
#' @param picard_jar_path path to picard jar file
#' @param samtools_exe path to samtools
#' 
#' @param ref_fasta_path reference fasta file
#' 
#' @param picard_markdup_opts a character vector of options for picard mark duplication step
#' @param gatk_target_opts a character vector of options for gatk target step
#' @param gatk_realign_opts a character vector of options for gatk realign step
#' @param gatk_baserecalib_opts a character vector of options for gatk baserecalib step
#' @param gatk_printreads_opts a character vector of options for gatk printreads step
#' 
#' 
#' @param mem_markdup memory used by java, example -Xmx1g
#' @param mem_target memory used by java, example -Xmx1g
#' @param mem_realign memory used by java, example -Xmx1g
#' @param mem_baserecalib memory used by java, example -Xmx1g
#' @param mem_printreads memory used by java, example -Xmx1g
#' 
#' 
#' @param cpu_markdup not used.
#' @param cpu_target number of threads used for GATK target creation step
#' @param cpu_realign number of cpu used
#' @param cpu_baserecalib number of cpu used
#' @param cpu_printreads number of cpu used
#' 
#' @param execute_cmds run commands, after creation. Useful for testing/debugging and running on local platforms.
#'
#' @export
#'
#' @examples \dontrun{
#' ## load options, including paths to tools and other parameters
#' opts_flow$load(flowr::fetch_conf("ultraseq.conf"), check = FALSE)
#' out = preprocess("my_wex.bam", samplename = "samp", split_by_chr = TRUE)
#'
#' }
preprocess.gatk_v2 <- function(bam,
                               outfile,
                               samplename = opts_flow$get("samplename"),
                               split_by_chr = opts_flow$get("split_by_chr"),
                               
                               java_exe = opts_flow$get("java_exe"),
                               java_tmp = opts_flow$get("java_tmp"),
                               
                               gatk_jar_path = opts_flow$get('gatk_jar_path'),
                               picard_jar_path = opts_flow$get('picard_jar_path'),
                               samtools_exe = opts_flow$get('samtools.exe'),
                               
                               cpu_markdup = 1,
                               mem_markdup = opts_flow$get("mem_markdup"),
                               
                               cpu_target = opts_flow$get("cpu_target"),  ## not used
                               mem_target = opts_flow$get("mem_target"),
                               
                               cpu_realign = opts_flow$get("cpu_realign"),
                               mem_realign= opts_flow$get("mem_realign"),
                               
                               ## scatter 8 per node nct=8
                               cpu_baserecalib = opts_flow$get("cpu_baserecalib"),
                               mem_baserecalib = opts_flow$get("mem_baserecalib"),
                               
                               ## scatter 8 per node nct=8
                               cpu_printreads = opts_flow$get("cpu_printreads"),
                               mem_printreads = opts_flow$get("mem_printreads"),
                               
                               ref_fasta_path = opts_flow$get('ref_fasta_path'),
                               
                               picard_markdup_opts = opts_flow$get('picard_markdup_opts'),
                               gatk_target_opts = opts_flow$get('gatk_target_opts'),
                               gatk_realign_opts = opts_flow$get('gatk_realign_opts'),
                               gatk_baserecalib_opts = opts_flow$get('gatk_baserecalib_opts'),
                               gatk_printreads_opts = opts_flow$get('gatk_printreads_opts'), 
                               
                               execute_cmds = FALSE
){
  
  check_args(ignore = "outfile")
  # source('~/Dropbox/public/flowr/ultraseq/ultraseq/R/bam_set.R')
  bamset = bam_set(bam = bam, 
                   outprefix = outfile, 
                   ref_fasta_path = ref_fasta_path, 
                   split_by_chr = split_by_chr)
  
  # get the name of the function
  pipename = match.call()[[1]]
  message("Generating a ", pipename, " flowmat for sample: ", samplename)
  
  # ------------ dedup; SINGLE FILE
  dedupbam <- paste0(bamset$outprefix, ".marked.bam")
  metricsfile <- paste0(bamset$outprefix, ".marked.metrics")
  cmd_markdup <- sprintf("%s %s -Djava.io.tmpdir=%s -jar %s MarkDuplicates INPUT=%s OUTPUT=%s METRICS_FILE=%s %s",
                         java_exe, mem_markdup, java_tmp, picard_jar_path, bamset$bam,
                         dedupbam, metricsfile, 
                         picard_markdup_opts)
  cmd_markdup
  
  ## ------------ realign; SINGLE FILE
  intervalsfiles <- paste0(bamset$outprefix, ".realign.intervals")
  ## ------------ do this for all chrs
  cmd_target <- sprintf("%s %s -Djava.io.tmpdir=%s -jar %s -T RealignerTargetCreator -R %s -I %s -o %s -nt %s %s",
                        java_exe, mem_target, java_tmp, gatk_jar_path, ref_fasta_path, dedupbam,
                        intervalsfiles, cpu_target, gatk_target_opts)
  
  realignedbams <- paste0(bamset$outprefix_chr ,".realigned.bam")
  cmd_realign <- sprintf("%s %s -Djava.io.tmpdir=%s -jar %s -T IndelRealigner -R %s -I %s -targetIntervals %s -o %s %s %s",
                         java_exe, mem_realign, java_tmp, gatk_jar_path, ref_fasta_path, dedupbam,
                         intervalsfiles, realignedbams, gatk_realign_opts, bamset$gatk_intervals)
  
  ## ------------ base recalibration
  ## explicity define intervals here as well, though not required.
  recalibbams <- paste0(bamset$outprefix_chr, ".recalibed.bam")
  recalibtabfile <- paste0(bamset$outprefix_chr, ".recalib.tab")
  cmd_baserecalib <- sprintf("%s %s -Djava.io.tmpdir=%s -jar %s -T BaseRecalibrator -R %s -I %s -o %s -nct %s %s %s",
                             java_exe, mem_baserecalib, java_tmp, gatk_jar_path, ref_fasta_path,
                             realignedbams, recalibtabfile, cpu_baserecalib,
                             gatk_baserecalib_opts, bamset$gatk_intervals)
  
  cmd_printreads1 <- sprintf("%s %s -Djava.io.tmpdir=%s -jar %s -T PrintReads -R %s -I %s -BQSR %s -o %s -nct %s %s %s",
                             java_exe, mem_printreads, java_tmp, gatk_jar_path, ref_fasta_path, realignedbams,
                             recalibtabfile, recalibbams, cpu_printreads,
                             gatk_printreads_opts, bamset$gatk_intervals)
  cmd_printreads2 <- sprintf("%s index %s", samtools_exe, recalibbams)
  cmd_printreads = sprintf("%s;%s", cmd_printreads1, cmd_printreads2)
  
  cmds <- list(markdup = cmd_markdup,
               target = cmd_target, realign = cmd_realign,
               baserecalib = cmd_baserecalib, printreads = cmd_printreads)
  sapply(cmds, length)
  
  
  if(execute_cmds)
    sapply(cmds, system)
  
  flowmat = to_flowmat(cmds, samplename = samplename)
  return(list(flowmat=flowmat, outfiles = recalibbams))
  
}

#' Flow following Broad's best practices for variant calling, starting from sorted bam
#'
#' @title Pre-process bam files following Broad's best practices for variant calling, starting from aligned BAM file
#' @description This function provides a wrapper around the best practices described on 
#' \href{https://software.broadinstitute.org/gatk/best-practices/workflow?id=11165}{GATK's website}.
#' If the link is broken, google 'GATK best practices'
#'
#' This aims to perform the following steps ( for DNA ):
#'
#' \itemize{
#' \item mark duplicates
#' \item realign indels
#' \item recalibrate bases
#' \item current version: \emph{3.4-46}
#' }
#'
#' For RNA GATK recommends a additional step of split n trim, which is not currently supported (contributions welcome !).
#'
#' \strong{NOTE}:
#'
#' Some GATK tools use \href{https://www.broadinstitute.org/gatk/guide/article?id=1975}{CPU threads while others use data threads},
#' flowr tries to use efficiently make the best use of both/either depending on tool's compatibility.
#' 
#' 1.   call MarkDuplicates ----> SortAndFixTags
#' 2. Perform Base Quality Score Recalibration (BQSR) on the sorted BAM in parallel (BaseRecalibrator)
#' 
#'
#' @param bam bam file path
#' @param samplename name of the sample
#' @param outfile output file name
#' 
#' @param java_exe path to java
#' @param java_tmp path to java_tmp, can leave blank
#' 
#' @param split_by_chr split processing by chromosomr where ever possible
#' 
#' @param gatk_jar_path path to gatk jar file
#' @param picard_jar_path path to picard jar file
#' @param samtools_exe path to samtools
#' 
#' @param ref_fasta_path reference fasta file
#' 
#' @param picard_markdup_opts a character vector of options for picard mark duplication step
#' @param gatk_target_opts a character vector of options for gatk target step
#' @param gatk_realign_opts a character vector of options for gatk realign step
#' @param gatk_baserecalib_opts a character vector of options for gatk baserecalib step
#' @param gatk_printreads_opts a character vector of options for gatk printreads step
#' 
#' 
#' @param mem_markdup memory used by java, example -Xmx1g
#' @param mem_target memory used by java, example -Xmx1g
#' @param mem_realign memory used by java, example -Xmx1g
#' @param mem_baserecalib memory used by java, example -Xmx1g
#' @param mem_printreads memory used by java, example -Xmx1g
#' 
#' 
#' @param cpu_markdup not used.
#' @param cpu_target number of threads used for GATK target creation step
#' @param cpu_realign number of cpu used
#' @param cpu_baserecalib number of cpu used
#' @param cpu_printreads number of cpu used
#' 
#' @param execute_cmds run commands, after creation. Useful for testing/debugging and running on local platforms.
#'
#' @export
#'
#' @examples \dontrun{
#' ## load options, including paths to tools and other parameters
#' opts_flow$load(flowr::fetch_conf("ultraseq.conf"), check = FALSE)
#' out = preprocess("my_wex.bam", samplename = "samp", split_by_chr = TRUE)
#'
#' }
preprocess.gatk_v4 <- function(bam,
                               outfile,
                               samplename = opts_flow$get("samplename"),
                               split_by_chr = opts_flow$get("split_by_chr"),
                               
                               samtools_exe = opts_flow$get('samtools.exe'),
                               
                               java_exe = opts_flow$get("java_exe"),
                               java_mem = opts_flow$get("java_mem_str"),
                               
                               ref_fasta = opts_flow$get('ref_fasta'),
                               
                               gatk_jar = opts_flow$get('gatk4.jar'),
                               gatk_opts = opts_flow$get("gatk4.opts"),
                               
                               
                               # dedup
                               #cpu_markdup = 1,
                               markdup_mem = opts_flow$get("gatk4.markdup.mem"),
                               markdup_opts = opts_flow$get('gatk4.markdup.opts'),
                               
                               ## scatter 8 per node nct=8
                               #cpu_baserecalib = opts_flow$get("gatk4.baserecalib.cpu"),
                               baserecalib_mem = opts_flow$get("gatk4.baserecalib.mem"),
                               baserecalib_opts = opts_flow$get('gatk4.baserecalib.opts'),
                               
                               execute_cmds = FALSE
){
  
  check_args(ignore = "outfile")

  # source('~/Dropbox/public/flowr/ultraseq/ultraseq/R/bam_set.R')
  bamset = bam_set(bam = bam, outfile = outfile, ref_fasta_path = ref_fasta, split_by_chr = split_by_chr)
  
  # get the name of the function
  pipename = match.call()[[1]]
  message("Generating a ", pipename, " flowmat for sample: ", samplename)
  
  # ------------ dedup; SINGLE FILE
  dedupbam <- paste0(bamset$outprefix, ".duplicates_marked.bam")
  metricsfile <- paste0(bamset$outprefix, ".duplicate_metrics")
  cmd_markdup <- glue("{java_exe} {gatk_opts} {markdup_mem} -jar {gatk_jar} MarkDuplicates ",
                      "--INPUT {bam} --OUTPUT {dedupbam} --METRICS_FILE {metricsfile} {markdup_opts}")
  cmd_markdup
  
  # ------------ base recalibration
  # https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_walkers_bqsr_BaseRecalibrator.php
  # First pass of the Base Quality Score Recalibration (BQSR) -- Generates recalibration table based on various
  #user-specified covariates (such as read group, reported quality score, machine cycle, and nucleotide context).
  #Version:4.0.10.1
  # explicity define intervals here as well, though not required.
  bqsr_reports <- paste0(bamset$outprefix_interval, ".recal_data.table")
  gatk_intervals = bamset$gatk_intervals
  cmd_baserecalib <- glue("{java_exe} {gatk_opts} {java_mem} -jar {gatk_jar} BaseRecalibrator -R {ref_fasta} ", 
                          "-I {dedupbam} -O {bqsr_reports} ",
                          "{baserecalib_opts} {gatk_intervals}")
  cmd_baserecalib[1]
  # ${gatk_path} --java-options "${java_opt}" BaseRecalibrator \
  # -R ${ref_fasta}  -I ${input_bam} \
  # --use-original-qualities \
  # -O ${recalibration_report_filename} \
  # --known-sites ${dbSNP_vcf} --known-sites ${sep=" --known-sites " known_indels_sites_VCFs} \
  # -L ${sep=" -L " sequence_group_interval}
  
  bqsr_reports_i <- paste0(bqsr_reports, collapse = " -I ");bqsr_reports_i
  bqsr_report_o = paste0(bamset$outprefix, ".recal_data.table")
  cmd_bqsr_reports <- glue("{java_exe} {gatk_opts} {java_mem} -jar {gatk_jar} GatherBQSRReports ", 
                           "-I {bqsr_reports_i} -O {bqsr_report_o}")
  cmd_bqsr_reports
  # ${gatk_path} --java-options "${java_opt}" \
  # GatherBQSRReports \
  # -I ${sep=' -I ' input_bqsr_reports} \
  # -O ${output_report_filename}
  # This tool performs the second pass in a two-stage process called Base Quality Score Recalibration (BQSR). Specifically, it recalibrates the base qualities of the input reads based on the recalibration table produced by the BaseRecalibrator tool, and outputs a recalibrated BAM or CRAM file
  recalibbams <- paste0(bamset$outprefix_interval, ".duplicates_marked.recalibrated.bam")
  cmd_apply_bqsr <- glue("{java_exe} {gatk_opts} {java_mem} -jar {gatk_jar} ApplyBQSR ", 
                         "-R {ref_fasta} -I {dedupbam} -O {recalibbams} ",
                         "{gatk_intervals} -bqsr {bqsr_reports} ",
                         "--static-quantized-quals 10 --static-quantized-quals 20 ",
                         "--static-quantized-quals 30 --add-output-sam-program-record ",
                         "--create-output-bam-md5 --use-original-qualities")
  cmd_apply_bqsr[1]
  # ${gatk_path} --java-options "${java_opt}" \
  # ApplyBQSR \
  # -R ${ref_fasta}  -I ${input_bam} -O ${output_bam_basename}.bam \
  # -L ${sep=" -L " sequence_group_interval} \
  # -bqsr ${recalibration_report} \
  # --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 \
  # --add-output-sam-program-record \
  # --create-output-bam-md5 --use-original-qualities
  # }
  
  # ** gather bams ----
  recalibbams_i = paste0(recalibbams, collapse = " -I ");recalibbams_i
  recalibbam_o = paste0(bamset$outprefix, ".duplicates_marked.recalibrated.bam");recalibbam_o
  cmd_gather_bams <- glue("{java_exe} {gatk_opts} {java_mem} -jar {gatk_jar} GatherBamFiles ", 
                          "-I {recalibbams_i} -O {recalibbam_o} ",
                          "--CREATE_INDEX true --CREATE_MD5_FILE true")
  cmd_gather_bams
  # ${gatk_path} --java-options "-Dsamjdk.compression_level=${compression_level} ${java_opt}" \
  # GatherBamFiles \
  # --INPUT ${sep=' --INPUT ' input_bams} \
  # --OUTPUT ${output_bam_basename}.bam \
  # --CREATE_INDEX true \
  # --CREATE_MD5_FILE true
  # }
  
  # cleanup job
  #fls_rm = c(recalibbams, bqsr_reports)
  
  # merge parallel jobs
  cmd_baserecalib[1];cmd_apply_bqsr[1]
  cmds <- list(markdup = cmd_markdup,
               baserecalib = paste(cmd_baserecalib, cmd_apply_bqsr, sep=";"),
               gather_bams = paste(cmd_bqsr_reports, cmd_gather_bams, sep = ";"))
  sapply(cmds, length)
  
  if(execute_cmds)
    sapply(cmds, system)
  
  flowmat = to_flowmat(cmds, samplename = samplename) %>% 
    mutate(cmd = as.character(cmd))
  ret = list(flowmat = flowmat, 
             outfiles = recalibbam_o,
             outfiles2 = c(
               metricsfile = metricsfile,
               dedupbam = dedupbam,
               bqsr_report = bqsr_report_o,
               recalibbam = recalibbam_o,
               recalibbai = gsub(".bam", ".bai", recalibbam_o)
             ))
  return(ret)
}


# ${gatk_path} --java-options "-Dsamjdk.compression_level=${compression_level} ${java_opt}" 
# MarkDuplicates \
# --INPUT ${sep=' --INPUT ' input_bams} \
# --OUTPUT ${output_bam_basename}.bam \
# --METRICS_FILE ${metrics_filename} \
# --VALIDATION_STRINGENCY SILENT \
# --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
# --ASSUME_SORT_ORDER "queryname" \
# --CREATE_MD5_FILE true
# }
# 
# ${gatk_path} --java-options "${java_opt}" \
# BaseRecalibrator \
# -R ${ref_fasta} \
# -I ${input_bam} \
# --use-original-qualities \
# -O ${recalibration_report_filename} \
# --known-sites ${dbSNP_vcf} \
# --known-sites ${sep=" --known-sites " known_indels_sites_VCFs} \
# -L ${sep=" -L " sequence_group_interval}
# }
# 
# ${gatk_path} --java-options "${java_opt}" \
# GatherBQSRReports \
# -I ${sep=' -I ' input_bqsr_reports} \
# -O ${output_report_filename}
# }
# 
# ${gatk_path} --java-options "${java_opt}" \
# ApplyBQSR \
# -R ${ref_fasta} \
# -I ${input_bam} \
# -O ${output_bam_basename}.bam \
# -L ${sep=" -L " sequence_group_interval} \
# -bqsr ${recalibration_report} \
# --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 \
# --add-output-sam-program-record \
# --create-output-bam-md5 \
# --use-original-qualities
# }
# 
# ${gatk_path} --java-options "-Dsamjdk.compression_level=${compression_level} ${java_opt}" \
# GatherBamFiles \
# --INPUT ${sep=' --INPUT ' input_bams} \
# --OUTPUT ${output_bam_basename}.bam \
# --CREATE_INDEX true \
# --CREATE_MD5_FILE true
# }


