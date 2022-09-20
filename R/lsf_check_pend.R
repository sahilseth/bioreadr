check_pend <- function(){
  npend = system("busers -w | awk '{print $5}' | tail -n1", intern = T) %>% as.integer()
  mp = system("busers -w | awk '{print $12}' | tail -n1", intern = T) %>% as.integer()
  if((mp + 1000 - npend) < 0){
    flog.info(glue("{npend} pending jobs, waiting for some of them to complete"))
    Sys.sleep(60)
    check_pend()
  }
  return("all set")
}