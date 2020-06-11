fastx_trimmer <- function(in.f1, 
                          in.f2,
                          trimmed_files1,
                          trimmed_files2,
                          
                          fastx_trimmer_exe,
                          
                          pe){
  
  trimmed_files1 = file.path(gsub(file_ext , ".trimmed1.fastq", basename(fqs1)))
  cmd_trim1 <- sprintf("%s/fastx_trimmer -i %s -o %s -f %s -l %s -Q %s",
                       fastx_dir, in.f1, trimmed_files1,
                       read1_start, read1_end, "Q33")
  cmd_trim <- cmd_trim1
  
  if(pe){
    trimmed_files2=file.path(gsub(".fastq.gz",".trimmed2.fastq", basename(fqs1)))
    cmd_trim2 <- sprintf("%s/fastx_trimmer -i %s -o %s -f %s -l %s -Q %s",
                         fastx_dir, in.f1, trimmed_files2,
                         ##from, to
                         read2_start,read2_end, "Q33")
    in.f2 = trimmed_files2
    cmd_trim <- paste(cmd_trim1, cmd_trim2, sep=";")
  }
  in.f1 = trimmed_files1
  cmds$trim = cmd_trim
}