

fq_set <- function(fq1, fq2){
  
  
}



#' A wrapper around somatic mutation caller MuTect
#'
#' This generates a set of commandlines, per chromosome
#'
#' @param tumor_bam path to a tumor bam file
#' @param normal_bam path to a normal bam file
#' @param samplename name of the sample, to be added to the flowmat
#' @param outprefix output file name [optional].
#' By default, this is created using names of tumor and normal bam files.
#' @param is_merged specify if the input bam is already split by chromosomes (FALSE) 
#' or is a merged file with all chromosome (TRUE). [FALSE]
#' @param split_by_chr fork running mutect by chromosome to make it faster [TRUE].
#' Turning is option OFF is not fully supported.
#' @param java_exe path to java's executable [java]
#' @param java_tmp path to a temp folder to be used by java
#' @param mutect_jar path to mutect's jar file
#' @param cpu_mutect integer specifying number of thread to be used per mutect fork
#' @param mem_mutect amount of memory to be used by mutect [-Xmx12g]
#' @param ref_fasta_path path to the reference genome in fasta format
#' @param mutect_opts additional arguments passed onto mutect
#' 
#' @return The function returns a flowmat with all the commands. 
#' The final file is called \code{'outprefix'_merged.mutect.tsv}.
#'
#' @import flowr
#' 
#' @export
#'
#' @examples \dontrun{
#'
#' x = "tumor.bam"
#' y = "normal.bam"
#'
#' out = mutect(x, y, samplename="tumor_normal", is_merged = TRUE)
#'
#' }
mutect1.7 <- function(trk,
                      samplename,
                      
                      variant_calling_mode,
                      
                      java_exe = opts_flow$get("java17_exe"),
                      java_tmp = opts_flow$get("java17_tmp"),
                      java_mem = opts_flow$get("java17_mem_str"),
                      # for gatherVcfs
                      gatk4_exe = opts_flow$get("gatk4_exe"),
                      
                      ref_fasta = opts_flow$get('ref_fasta'),
                      
                      mutect1_jar = opts_flow$get('mutect1_jar'),
                      mutect1_cpu = opts_flow$get('mutect1_cpu'), # not used
                      mutect1_mem_str = opts_flow$get("mutect1_mem_str"),
                      
                      mutect1_opts = opts_flow$get('mutect1_opts'),
                      
                      mutect1_germline_vcf = opts_flow$get('mutect1_germline_vcf'),
                      mutect1_cosmic_vcf = opts_flow$get('mutect1_cosmic_vcf'),
                      mutect1_pon_vcf = opts_flow$get("mutect1_pon_vcf"),
                      mutect1_dir = pipe_str$mutect1$dir,

                      # previous option, we will ALWAYS assume this is merged,
                      # ONE bam per sample, not multiple!
                      # is_merged = TRUE,
                      # split_by = TRUE, # always split by pre-defined intervals
                      # default is WEX, can change this
                      interval_split = opts_flow$get("mutect1_intervals")    
){
  
  # expect ALL non null args
  check_args(ignore = "mutect1_pon_vcf")
  flog.debug(paste0("Generating a mutect1.7 flowmat for sample: ", samplename))

  # check generic columns
  trk = metadata_for_mutect1(trk)
  # file names are now, matched, OR pon
  # we are not keeping the ref info in file names now...
  expect_columns(trk, c("outprefix", "outprefix_paired"))
  # ** prep files ----------
  trk %<>% 
    mutate(mutect1_outprefix = file.path(mutect1_dir, outprefix), 
           mutect1_outprefix_paired = file.path(mutect1_dir, outprefix_paired),
           mutect1_vcf = paste0(mutect1_outprefix_paired, "_mutect.merged.vcf.gz"))
  
  trk_tum = filter(trk, normal == "NO")
  trk_norm = filter(trk, normal == "YES")

  
  # for now this will ONLY work for single tumor sample
  # this maybe needs to be a loop!
  i=1
  tmp = lapply(1:nrow(trk_tum), function(i){

    tumor_bam = trk_tum$bam[i]
    # assuming a SINGLE normal!!
    normal_bam = trk_norm$bam[1]
    outprefix = trk_tum$mutect1_outprefix_paired[i]
    
    bamset = bam_set(bam = tumor_bam, 
                      outprefix = outprefix,
                      ref_fasta = ref_fasta, 
                      interval_split = interval_split,
                      split_by = "interval_split")
    
    # pairedset = paired_bam_set(tumor_bam = tumor_bam, normal_bam = normal_bam, 
    #                            outprefix = outprefix, 
    #                            is_merged = is_merged, split_by_chr = split_by_chr)
    
    # pipename = match.call()[[1]]
    
    mutects <- paste0(bamset$outprefix_interval, ".mutect.txt")
    vcfs <- paste0(bamset$outprefix_interval, ".mutect.vcf")
    wigs <- paste0(bamset$outprefix_interval, ".wig")
    
    # lapply(list(tumor_bam, normal_bam, bamset$outprefix_chr, bamset$chrs_names), length)
    
    intervals = bamset$intervals
    if(!is.null(mutect1_pon_vcf) & mutect1_pon_vcf != "" )
      mutect1_opts = glue("{mutect1_opts} --normal_panel {mutect1_pon_vcf}")
    # hg19_cosmic_v54_120711.vcf and dbsnp_132_b37.leftAligned.vcf
    
    # https://gatkforums.broadinstitute.org/gatk/discussion/2226/cosmic-and-dbsnp-files-for-mutect
    # Sites that are in dbSNP and COSMIC do NOT use the prior as a site being germline during somatic classification. 
    # This is because dbSNP contains a number of sites that are common somatic events which were deposited into dbSNP in the past. 
    # We want to counteract this effect and not make these sites harder to call
    # Sites in COSMIC are exempt from the "Panel of Normals" filter -- again these are typically recurrent events and this is a 
    # mechanism to bypass this filter if necessary
    # Sites in the output call_stats file are annotated as "COSMIC"
    # germline_vcf = "/home/sseth/ref/human/b37/annotations/gatk_bundle/dbsnp_138.b37.excluding_sites_after_129.vcf.gz"
    # cosmic_vcf = "/home/sseth/ref/human/b37/annotations/gatk_bundle/hg19_cosmic_v54_120711.vcf"
    # For input string: "R", for input source: /rsrch3/home/iacs/sseth/ref/human/b37/annotations/broad-somatic-b37/af-only-gnomad.raw.sites.vcf.gz
    # The analysis MuTect currently does not support parallel execution with nt.  Please run your analysis without the nt option.                   
    cmd_mutect <- glue("{java_exe} {mutect1_mem_str} -Djava.io.tmpdir={java_tmp} -jar {mutect1_jar} -T MuTect --reference_sequence {ref_fasta} ",
                      "--input_file:tumor {tumor_bam} --input_file:normal {normal_bam} ",
                      "--out {mutects} --vcf {vcfs} --coverage_file {wigs} ",
                      "--dbsnp {mutect1_germline_vcf} ",
                      "--cosmic {mutect1_cosmic_vcf} ",
                      "{mutect1_opts} {intervals}")
    if(variant_calling_mode == "pon")
      cmd_mutect <- glue("{java_exe} {mutect1_mem_str} -Djava.io.tmpdir={java_tmp} -jar {mutect1_jar} -T MuTect --reference_sequence {ref_fasta} ",
                    "--input_file:tumor {tumor_bam} ",
                    "--out {mutects} --vcf {vcfs} --coverage_file {wigs} ",
                    "--dbsnp {mutect1_germline_vcf} ",
                    "--cosmic {mutect1_cosmic_vcf} ",
                    "{mutect1_opts} {intervals}")


    # test:
    # cmd_mutect[1:2]

    # .filter='judgement==KEEP'
    # in case of a single file, this mean a read and write operation
    mutect1_tsv = paste0(bamset$outprefix, "_mutect.merged.tsv.gz")
    mutect1_tsv_pass = paste0(bamset$outprefix, "_mutect.merged.pass.tsv.gz")
    mutect1_vcf = paste0(bamset$outprefix, "_mutect.merged.vcf.gz")
    cmd_mergetsv = sprintf("flowr ultraseq::merge_sheets x=%s outfile=%s",
                        paste(mutects, collapse = ","), mutect1_tsv)
    # ** merge VCFs ---------
    # mutect_vcf = paste0(outprefix, ".vcf.gz")
    mutect_vcfs_i = paste0(vcfs, collapse = " -I ")
    # GatherVcfs
    # https://software.broadinstitute.org/gatk/documentation/tooldocs/current/picard_vcf_GatherVcfs.php
    # gather should be faster
    cmd_mergevcf = glue("{gatk4_exe} --java-options {java_mem} GatherVcfs -I {mutect_vcfs_i} -O {mutect1_vcf}; ")
    # we dont have stats file for mutect1!
    #   "bash ~/Dropbox/public/flow-r/my.ultraseq/pipelines/bin/m2_merge_vcf_stats.sh {merged_mutect_vcf}.stats"
    
    # convert to vcf by default
    # cmd_vcf = glue("funr my.ultraseq:::to_vcf.mutect_call_stats x={merged_mutect_tsv} outfile={merged_mutect_vcf}")
    # --artifact_detection_mode: used when running the caller on a normal (as if it were a tumor) to detect artifacts
    
    # @sbamin, create a filtered file by default
    keep_mutect <- glue("zgrep -E 'contig|KEEP' {mutect1_tsv} > {mutect1_tsv_pass}")

    # cmd_merge = paste(cmd_mergetsv, keep_mutect, cmd_mergevcf, sep = ";")
    
    cmds_mrg = c(cmd_mergetsv, cmd_mergevcf, keep_mutect)
    
    #if(execute_cmds) sapply(cmds, system)
    
    flowmat = list(
      mutect1.splt = cmd_mutect,
      mutect1.mrg = cmds_mrg) %>% to_flowmat(samplename = samplename)
    lst = list(flowmat = flowmat, 
                outfiles = list(mutect1_tsv = mutect1_tsv, 
                                mutect1_vcf = mutect1_vcf))
  })
  lst = purrr::transpose(tmp)
  lst$flowmat %<>% bind_rows()
  lst$trk = trk

  return(lst)
  
}


# TEST:
# bsub.iacs
# cd /rsrch3/home/iacs/sseth/flows/SS/sarco/jco/wex/runs/sarco2/tmp

# run 2 mutect cmds (to test merging etc)
# https://gatkforums.broadinstitute.org/gatk/discussion/5810/error-message-unable-to-parse-header-with-error-for-input-string-r
"The Number entry is an Integer
that describes the number of values that can be included with the INFO field. For example, if the INFO field contains
a single number, then this value should be 1; if the INFO field describes a pair of numbers, then this value should
be 2 and so on. There are also certain special characters used to define special cases:"
"If the field has one value for each possible allele (including the reference), then this value should be `R'"
# ...so that is a valid value?
#  any idea how to make GATK accept my files, without significantly changing the meaning of their content?
#  if i change the "R"s to "4"s (for the 4 possible-possible alleles - i haven't called indels), it seems to work.
# do you think i've opened up any problems for myself in the future?
               
# The analysis MuTect currently does not support parallel execution with nt.  Please run your analysis without the nt option.                   
# module load jdk/1.7.0_79;java -Xmx16g -Djava.io.tmpdir=java_tmp -jar /risapps/noarch/mutect/1.1.7/mutect-1.1.7.jar -T MuTect --reference_sequence /rsrch3/home/iacs/sseth/ref/human/b37/fastas/Homo_sapiens_assembly19.fasta --input_file:tumor /rsrch3/home/iacs/sseth/flows/SS/sarco/jco/wex/bam/JCOxSARCO-sarco2-T_1208XX_ST1374_073_H09WGADXX.duplicates_marked.recalibrated.bam --input_file:normal /rsrch3/home/iacs/sseth/flows/SS/sarco/jco/wex/bam/JCOxSARCO-sarco2-N_1208XX_ST1374_073_H09WGADXX.duplicates_marked.recalibrated.bam --out WEX-sarco2-T___matched_001.mutect.txt --vcf WEX-sarco2-T___matched_001.mutect.vcf --coverage_file WEX-sarco2-T___matched_001.wig --dbsnp /home/sseth/ref/human/b37/annotations/gatk_bundle/dbsnp_138.b37.excluding_sites_after_129.vcf.gz --cosmic /rsrch3/home/iacs/sseth/ref/human/b37/annotations/cosmic/v88/CosmicMuts_sorted.vcf --enable_extended_output   -L /rsrch3/home/iacs/sseth/ref/human/b37/annotations/broad-somatic-b37/whole_exome_agilent_1.1_refseq_plus_3_boosters.Homo_sapiens_assembly19.targets.pad250.split200/0000-scattered.interval_list
# --normal_panel /rsrch3/home/iacs/sseth/flows/SS/sarco/jco/wex/ponm/ponm.vcf.gz 
# module load jdk/1.7.0_79;java -Xmx16g -Djava.io.tmpdir=java_tmp -jar /risapps/noarch/mutect/1.1.7/mutect-1.1.7.jar -T MuTect --reference_sequence /rsrch3/home/iacs/sseth/ref/human/b37/fastas/Homo_sapiens_assembly19.fasta --input_file:tumor /rsrch3/home/iacs/sseth/flows/SS/sarco/jco/wex/bam/JCOxSARCO-sarco2-T_1208XX_ST1374_073_H09WGADXX.duplicates_marked.recalibrated.bam --input_file:normal /rsrch3/home/iacs/sseth/flows/SS/sarco/jco/wex/bam/JCOxSARCO-sarco2-N_1208XX_ST1374_073_H09WGADXX.duplicates_marked.recalibrated.bam --out WEX-sarco2-T___matched_002.mutect.txt --vcf WEX-sarco2-T___matched_002.mutect.vcf --coverage_file WEX-sarco2-T___matched_002.wig --dbsnp /home/sseth/ref/human/b37/annotations/gatk_bundle/dbsnp_138.b37.excluding_sites_after_129.vcf.gz --cosmic /rsrch3/home/iacs/sseth/ref/human/b37/annotations/cosmic/v88/CosmicMuts_sorted.vcf --enable_extended_output -L /rsrch3/home/iacs/sseth/ref/human/b37/annotations/broad-somatic-b37/whole_exome_agilent_1.1_refseq_plus_3_boosters.Homo_sapiens_assembly19.targets.pad250.split200/0001-scattered.interval_list


# funr my.ultraseq:::to_vcf.mutect_call_stats x=WEX-sarco2-T___matched_200.mutect.txt outfile=WEX-sarco2-T___matched_200.mutect.vcf





# END

