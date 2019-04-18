


#' superfreq
#' 
#' call superFreq from the pkg and store the results
#'
#' @param odir dir to store final results
#' @param ... all params passed onto superFreq::superFreq
#'
#' @export
superfreq <- function(odir = ".", ...){
  
  library(superFreq)
  
  # need to convert some params to numeric:
  # BQoffset=33 cpus=8 systematicVariance=0.02 maxCov=150
  
  # get all the args for this func
  args <- formals(superFreq::superFreq)
  
  # fetch all args passed from ...
  arguments <- list(...)
  
  out = try(withVisible(do.call(func, args = params)))
  
  out = superFreq::superFreq(...)
  write_rds(out, glue("{odir}/R/sfq_results.rds"))
  
}


read_sfq.trk.seg <- function(df_trk, 
                            col_fl, 
                            col_samplename){
  
  df_trk = data.frame(df_trk, stringsAsFactors = F)
  df_seg = lapply(1:nrow(df_trk), function(i){
    
    samplename = df_trk[i, col_samplename]
    seg_fl = df_trk[i, col_fl]
    message(samplename, " ", appendLF = F)
    
    df_seg = read_tsv(seg_fl, 
                      col_types = cols(
                        chr = col_character(),
                        start = col_double(),
                        end = col_double(),
                        x1 = col_double(),
                        x2 = col_double(),
                        M = col_double(),
                        width = col_double(),
                        df = col_double(),
                        var = col_double(),
                        cov = col_double(),
                        Nsnps = col_double(),
                        pHet = col_double(),
                        pAlt = col_double(),
                        odsHet = col_double(),
                        f = col_double(),
                        stat = col_double(),
                        nullStat = col_double(),
                        altStat = col_double(),
                        nullStatErr = col_double(),
                        altStatErr = col_double(),
                        postHet = col_double(),
                        ferr = col_double(),
                        call = col_character(),
                        clonality = col_double(),
                        clonalityError = col_double(),
                        sigma = col_double(),
                        pCall = col_double(),
                        subclonality = col_double(),
                        subclonalityError = col_double(),
                        genes = col_character()
                      )) %>% clean_names()
    df_seg %<>% mutate(samplename = samplename)
    df_seg
  }) %>% bind_rows()
  df_seg
  
}