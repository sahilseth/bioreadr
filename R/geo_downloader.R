

# https://bioconductor.org/packages/release/bioc/vignettes/GEOquery/inst/doc/GEOquery.html
# http://einstein:2500/notebooks/projects/packs_my.db/misc/geo/download_process_geo.ipynb#

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
  
  gseid = "GSE87455"
  
  
  
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

if(FALSE){
  geodir = "~/rsrch1_data/public_data/geo"
  source('~/Dropbox/public/flowr/my.ultraseq/my.ultraseq/R/geo_downloader.R')
  
  gseid = "GSE87455"
}

geo_download.illumina <- function(gseid, geodir, 
                                  download_suppl = T,
                                  read_suppl = F
                                  ){
  library(pacman)
  p_load(GEOquery, dplyr, magrittr, readr, glue)
  
  # get the series matrix
  message("download matrix")
  # https://bioconductor.org/packages/release/bioc/vignettes/GEOquery/inst/doc/GEOquery.html
  gse_mat <- getGEO(gseid, destdir = geodir, GSEMatrix = T, AnnotGPL = FALSE, getGPL = TRUE)

  message("extract phenodata")
  df_trk = phenoData(gse_mat[[1]])@data# %>% tbl_df()
  # df_trk = lapply(gse_mat, function(x){ phenoData(x)@data }) %>% 
  #   do.call(rbind, .)
  df_trk = lapply(gse_mat, function(x){ phenoData(x)@data }) %>% 
    bind_rows() %>% data.frame()
  rownames(df_trk) = df_trk$title
  
  # to download mats
  message("download ")
  gse <- getGEO(gseid, destdir = geodir, GSEMatrix = F, AnnotGPL = FALSE, getGPL = TRUE)
  # class(gse[[1]])
  length(gse)
  # platform:
  gse@gpls
  
  message("download cell files")
  fls_supp = getGEOSuppFiles(gseid, baseDir = geodir)
  
  message("read supp fls")
  if(read_suppl)
    lst_mats = geo_read_supp_mats(fls_supp)
  
  # create a mae
  p_load(MultiAssayExperiment, SummarizedExperiment)
  # confirm the order is correct:
  # as.character(df_trk$title) == colnames(lst_mats[[1]])
  
  mae = MultiAssayExperiment(experiments = lst_mats, 
                             colData = df_trk)
  write_rds(mae, glue("{geodir}/{gseid}/mae.rds"))
  
    
}

geo_download <- function(gseid, geodir){
  
  library(pacman)
  p_load(GEOquery, dplyr, magrittr, readr, glue)
  
  # get the series matrix
  message("download matrix")
  # https://bioconductor.org/packages/release/bioc/vignettes/GEOquery/inst/doc/GEOquery.html
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
  
  # to download mats
  message("download ")
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

#' @param fls is directly from `getGEOSuppFiles`; df with 
untar_raw_fls <- function(df_fl){
  fl = rownames(df_fl)
  message("files (will use only 1st)", fl)
  fl = fl[1]
  fl = tools::file_path_as_absolute(fl)
  #getOption("unzip")
  # requires SINGLE file
  tmp = untar(fl, exdir = dirname(fl))
  fls_unzip = list.files(path = dirname(fl), full.names = T)
  fls_unzip = setdiff(fl)
  fls_unzip
}

# convert the eset into SummarizedExperiment
to_se.eset <- function(es){
  # add to metadata
  value_description = gse@gsms[[1]]@dataTable@columns %>% filter(Column == "VALUE") %>% 
    pull(Description) %>% as.character()
  
  mat = exprs(out$es[[1]])
  # extract coldata
  phen = phenoData(out$es[[1]])
  df_coldata = phen@data  %>% clean_names()
  df_coldata_desc = phen@varMetadata
  # extract rowdata
  feat = featureData(out$es[[1]])
  df_rowdata = feat@data  %>% clean_names()
  df_rowdata_desc = feat@varMetadata
  # create SummarizedExperiment
  se = SummarizedExperiment(list(mat = mat), 
                            rowData = df_rowdata, colData = df_coldata, 
                            metadata = list(df_coldata_desc = df_coldata_desc,
                                            df_rowdata_desc = df_rowdata_desc,
                                            value_description = value_description))
  se
}






# END
