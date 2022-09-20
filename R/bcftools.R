bcf_rename_samples <- function(infl, outfl, from, to){
  # create a tmp files with samples
  tmpfl = tempfile()
  data.frame(from, to) %>% 
    write.table(tmpfl, col.names = FALSE, sep = "\t", row.names=F, quote = F)
  cmd = glue("module load bcftools;bcftools reheader -s {tmpfl} {infl} > {outfl}")
  system(cmd)

}