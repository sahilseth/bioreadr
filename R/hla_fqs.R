#' hla_fqs
#'
#' @param bam 
#' @param samplename 
#' @param region 
#' @param execute 
#'
#' @import glue
#' @export
hla_fqs <- function(bam, samplename, 
                    
                    region = "6:29690552-33111102", 
                    odir = "hla_fqs",
                    samtools_exe = opts_flow$get("samtools_exe"),
                    java_exe = opts_flow$get("java_exe"),
                    java_mem = opts_flow$get("java_mem"),
                    picard_jar = opts_flow$get("picard_jar"),
                    
                    execute = TRUE # not used
){
  
  bamnm = basename(bam) %>% file.path(odir, .)
  bam_hla = gsub(".bam$", "_hla.bam", bamnm)
  bam_umap = gsub(".bam$", "_unmapped.bam", bamnm)
  
  fq_hla_1 = gsub("_hla.bam$", "_hla_1.fq", bam_hla)
  fq_hla_2 = gsub("_hla.bam$", "_hla_2.fq", bam_hla)
  fq_hla_3 = gsub("_hla.bam$", "_hla_3.fq", bam_hla)
  
  fq_umap_1 = gsub("_unmapped.bam$", "_umap_1.fq", bam_umap)
  fq_umap_2 = gsub("_unmapped.bam$", "_umap_2.fq", bam_umap)
  fq_umap_3 = gsub("_unmapped.bam$", "_umap_3.fq", bam_umap)
  
  fq1 = gsub("_hla.bam$", "_hla_umap_1.fq", bam_hla)
  fq2 = gsub("_hla.bam$", "_hla_umap_2.fq", bam_hla)
  fq3 = gsub("_hla.bam$", "_hla_umap_3.fq", bam_hla)
  
  # get bams  
  cmd_bam1 = glue("mkdir {odir};{samtools_exe} view -b {bam} {region} > {bam_hla};{samtools_exe} index {bam_hla}")
  cmd_bam2 = glue("{samtools_exe} view -b -f 12 {bam} > {bam_umap};{samtools_exe} index {bam_umap}")
  
  # get fqs
  cmd_fq1 <- glue("{java_exe} {java_mem} -jar {picard_jar} SamToFastq INPUT={bam_hla} ", 
                  "FASTQ={fq_hla_1} SECOND_END_FASTQ={fq_hla_2} UNPAIRED_FASTQ={fq_hla_3}")
  cmd_fq2 <- glue("{java_exe} {java_mem} -jar {picard_jar} SamToFastq INPUT={bam_umap} ", 
                  "FASTQ={fq_umap_1} SECOND_END_FASTQ={fq_umap_2} UNPAIRED_FASTQ={fq_umap_3}")
  
  cmd_fq = glue("cat {fq_hla_1} {fq_umap_1} > {fq1};", 
                "cat {fq_hla_2} {fq_umap_2} > {fq2};",
                "cat {fq_hla_3} {fq_umap_3} > {fq3}")
  
  cmds = list(hla_fq = c(cmd_bam1, cmd_bam2,
                         cmd_fq1, cmd_fq2, cmd_fq))
  flowmat = to_flowmat(cmds, samplename)
  
  list(flowmat = flowmat, outfiles = list(fq1 = fq1, fq2 = fq2, fq3 = fq3))
}
