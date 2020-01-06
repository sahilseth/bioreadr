

# https://bioconductor.org/packages/release/bioc/vignettes/GEOquery/inst/doc/GEOquery.html


if(FALSE){
  
  p_load(fs)
  # Single cell RNA sequencing of 1,534 cells in six fresh triple negative breast cancer tumors.
  gseid = "GSE118389"
  geodir = "~/rsrch1_data/public_data/geo"
  fs::dir_ls(geodir)
  
  # In total, 3,678 patients with BC were studied. For 405 tumors, 
  # a comprehensive multi-rater histopathologic evaluation was performed.
  # This SuperSeries is composed of the SubSeries GSE81538 [cohort 405] and GSE96058 [cohort 3273] linked to below.
  gseid = "GSE81540" # super series
  gseid = "GSE81538"
  gseid = "GSE96058"
  
}

gse_summary <- function(gseid){
  
  message("there are ", length(gse@header$sample_id), " samples")
  
    
}

geo_read_supp_mats <- function(fls_supp, 
                               rownames = T, 
                               sep = "\t"){
  
  # fls = fls_supp %>% rownames()
  x = fls_supp[1]
  lst_mats = lapply(fls_supp, function(x){
    
    message("reading ", x)
    # this would not work, its a matrix!
    mat = data.table::fread(x) %>% data.frame() %>% wranglr::to_mat()
    # if(!rownames){
    #   mat = read.delim(x, sep = sep)
    # }else{
    #   mat = read.table(x, row.names = 1, sep = sep) %>% data.matrix()
    # }
    
    })
  # check
  # mat2[1:5, 1:5]
  
  # name the lists
  names(lst_mats) = basename(fls_supp)
  
  # lets see the dimensions
  lapply(lst_mats, dim)
  
  lst_mats
}



geo_download <- function(gseid, geodir){
  
  library(pacman)
  p_load(GEOquery, dplyr, magrittr, readr, glue)
  
  # get the series matrix
  gse_mat <- getGEO(gseid, destdir = geodir, GSEMatrix = T, AnnotGPL = FALSE, getGPL = TRUE)

  # extract phenodata
  df_trk = phenoData(gse_mat[[1]])@data# %>% tbl_df()
  # df_trk = lapply(gse_mat, function(x){ phenoData(x)@data }) %>% 
  #   do.call(rbind, .)
  df_trk = lapply(gse_mat, function(x){ phenoData(x)@data }) %>% 
    bind_rows() %>% data.frame()
  rownames(df_trk) = df_trk$title
  
  # change the rowname
  rownames(df_trk) = df_trk$title
  
  # to download suppl files
  gse <- getGEO(gseid, destdir = geodir, GSEMatrix = F, AnnotGPL = FALSE, getGPL = TRUE)
  # class(gse[[1]])
  length(gse)
  # platform:
  gse@gpls
  
  # download cel files
  fls_supp = getGEOSuppFiles(gseid, baseDir = geodir)
  
  # read supp fls
  lst_mats = geo_read_supp_mats(fls_supp)
  
  # create a mae
  p_load(MultiAssayExperiment, SummarizedExperiment)
  # confirm the order is correct:
  # as.character(df_trk$title) == colnames(lst_mats[[1]])
  
  mae = MultiAssayExperiment(experiments = lst_mats, 
                             colData = df_trk)
  write_rds(mae, glue("{geodir}/{gseid}/mae.rds"))
  
  # working with expression sets
  # show(gse)
  # bad idea for single cell data etc!
  
  
  # example of one sample:
  class(gse@gsms)
  # mat = exprs(gse)

  # extract the matrix
  # mat = exprs(gse$)
  
  
  # boxplot(mat)
  
  
}

