

fq_set <- function(fq1, fq2){
  
  
}



#' A wrapper around somatic mutation caller MuTect
#'
#' This generates a set of commandlines, per chromosome
#'
#' @param tumor_bam path to a tumor bam file
#' @param normal_bam path to a normal bam file
#' @param samplename name of the sample, to be added to the flowmat
#' @param out_prefix output file name [optional].
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
#' The final file is called \code{'out_prefix'_merged.mutect.tsv}.
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
mutect2 <- function(tumor_bam,
                    normal_bam,
                    samplename = opts_flow$get("samplename"),
                    out_prefix,
                    
                    is_merged = TRUE,
                    split_by_chr = TRUE,
                    
                    java_exe = opts_flow$get("java_exe"),
                    java_tmp = opts_flow$get("java_tmp"),
                    
                    gatk_jar = opts_flow$get('gatk_jar'),
                    
                    cpu_mutect = opts_flow$get('cpu_mutect'), ## curently not suported
                    mem_mutect = opts_flow$get("java_mem"),
                    
                    ref_fasta_path = opts_flow$get('ref_fasta_path'),
                    
                    mutect_opts = opts_flow$get('mutect_opts')
                    
                    
){
  
  pairedset = paired_bam_set(tumor_bam = tumor_bam, normal_bam = normal_bam, 
                             out_prefix = out_prefix, 
                             is_merged = is_merged, split_by_chr = split_by_chr)
  
  check_args(ignore = "out_prefix")
  
  pipename = match.call()[[1]]
  message("Generating a ", pipename, " flowmat for sample: ", samplename)
  
  mutects <- paste0(pairedset$out_prefix_chr, ".mutect.txt")
  wigs <- paste0(pairedset$out_prefix_chr, ".wig")
  
  lapply(list(tumor_bam, normal_bam, pairedset$out_prefix_chr, pairedset$chrs_names), length)
  
  # mutect paired end mode
  mutect_opts = "--artifact_detection_mode -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS -L $regions_bed_fl"
  cmd_mutect <- sprintf("%s %s %s -T MuTect2 -R %s -I:tumor %s -I:normal %s -o %s -nct %s",
                        java_exe, mem_mutect, gatk_jar,
                        ref_fasta_path,
                        tumor_bam, normal_bam,
                        mutects, mutect_cores)
                        
                        
  #gatk -T MuTect2 -R $ref_fa -I:tumor $bam --artifact_detection_mode -L $regions_bed_fl -o mutect_59.vcf -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS -nct 20
  
  cmd_mutect <- sprintf("%s %s -Djava.io.tmpdir=%s -jar %s --T MuTect2 --reference_sequence %s --input_file:tumor %s --input_file:normal %s --out %s  --coverage_file %s %s %s",
                        java_exe, mem_mutect, java_tmp, gatk_jar, ref_fasta_path, 
                        tumor_bam, normal_bam,
                        mutects, wigs,
                        mutect_opts, pairedset$gatk_intervals)
  cmds <- list(mutect = cmd_mutect)
  
  # .filter='judgement==KEEP'
  # in case of a single file, this mean a read and write operation
  merged_mutect = paste0(pairedset$out_prefix, "_mutect.merged.tsv")
  cmd_merge1 = sprintf("flowr ultraseq::merge_sheets x=%s outfile=%s",
                       paste(mutects, collapse = ","), merged_mutect)
  
  # @sbamin, create a filtered file by default
  merged_filt = paste0(pairedset$out_prefix, "_mutect.merged.keep.tsv")
  keep_mutect <- sprintf("grep -E 'contig|KEEP' %s > %s_calls.KEEP.tsv",
                         merged_mutect, merged_filt)
  
  cmd_merge = paste(cmd_merge1, keep_mutect, sep = ";")
  
  cmds = c(cmds, mutect_merge = cmd_merge)
  
  #if(execute_cmds) sapply(cmds, system)
  
  flowmat = to_flowmat(cmds, samplename = samplename)
  return(list(flowmat=flowmat, outfiles=list(all = mutects, merged = merged_mutect)))
  
}


