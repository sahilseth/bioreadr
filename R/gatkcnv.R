


# read funcs -------
read_counts <- function(x){
  df <- data.table::fread(cmd = glue("grep -v '@' {x}"), data.table = F) %>%
      as_tibble() %>%
      clean_names()
}
read_counts.gatk = read_counts
read_gatk.counts = read_counts


 



# END
