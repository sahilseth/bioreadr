


bedtools_bam2fastq <- function(bam, samplename, fqgz,
                               bedtools_exe = opts_flow$get("bedtools_exe")){
  
  cmd = glue("{bedtools_exe} bamtofastq -i {bam} -fq /dev/stdout | gzip > {fqgz}") %>% as.character()
  
  
  flowmat = to_flowmat(list(bam2fq = cmd), samplename)

  list(flowmat = flowmat, outfiles = list(fqgz = fqgz))
  
}
