extract_samplename <- function(bam){
  
  cmd = glue("module load samtools;samtools view -H {bam} | grep '^@RG'")
  tmp = system(cmd, intern = T)
  tmp %<>% strsplit(split = "\t") %>% unlist() %>% 
    grep("SM", ., value = TRUE) %>% 
    strsplit(split = "SM:") %>% unlist()
  tmp[tmp != ""] %>% unique()
}