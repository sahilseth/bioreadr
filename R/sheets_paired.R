#' @name read_paired_samplesheet
#' 
#' @description read_paired_samplesheet
#'
#' @param x x
#' @param ... ignored
#' 
#' @export
#' 
read_paired_samplesheet <- function(x, ...){
    ext <- file_ext(x)
    if(ext=="tsv"){
        mat <- read.table(x, as.is=TRUE, sep="\t", header=TRUE)
    }else if(ext=="csv"){
        mat <- read.csv2(x, as.is=TRUE, comment.char = '#', strip.white=TRUE,
                         blank.lines.skip=TRUE, sep=",", header=TRUE)
    }
    mat <- mat[! (is.na(mat$samplename) | mat$samplename %in% c("NA", "NULL", "")), ]
    ## convert paired_mat for xenome, calling snps etc
    single_mat <- rbind(cbind(mat$project, mat$samplename, mat$sampbam, mat$db_sampleid),
                        cbind(mat$project, mat$refname, mat$refbam, mat$db_refid))
    colnames(single_mat) <- c("project", "samplename", "bam", "db_id")
    single_mat <- data.frame(unique(single_mat),stringsAsFactors=FALSE)
    return(list(single_mat=single_mat, paired_mat=mat))
}

read_paired_bam_trk <- function(trk, test_file_exists = T){
  
  pacman::p_load(testit, magrittr)
  
  if(!is.data.frame(trk))
    trk = readr::read_tsv(trk) %>% data.frame(stringsAsFactors = F)
  
  # TESTS --------
  # make sure trk has the reqd columns
  cols_df = colnames(df_paired_bam)
  cols_reqd = c("patient_id", "samp1_bam", "samp2_bam", "samp1", "samp2", "pairnm")
  # testthat::expect_(cols_expected %in% colnames(trk))
  testit::assert("we have reqd columns in trk", {
    cols_reqd %in% cols_df
  })
  
  # make sure files exist
  bams = c(trk$samp1_bam, trk$samp2_bam)
  bais = gsub(".bam$", ".bam.bai", bams)
  
  if(test_file_exists)
    testit::assert("check if all files exists", {
      file.exists(bams) &
        file.exists(bais)
    })
  trk
}

## x is a object from seqan_tooling
#' @name create_tooling_paired_samplesheet
#' @description create_tooling_paired_samplesheet
#' @param x x
#' @param outfile outfile
#' @param tumor.only tumor.only
#' @param normal.only normal.only
#' 
#' @export
#' 
create_tooling_paired_samplesheet <- function(x, outfile, tumor.only=FALSE, normal.only=FALSE){
    colnames(x)  <- tolower(colnames(x))
    project=apply(x[,c("project","subproject")], 1, paste, collapse="-")
    if(tumor.only){
        out_mat <- data.frame(project, samplename=x$tumorsampleid, refname=NA,
                              sampbam=x$tumorsamplepreprocbam, refbam=NA,
                              db_sampleid=x$tumorsamplerunitemid, db_refid=0)
    }
    if(!missing(outfile)){
        dir.create(dirname(outfile),recursive=TRUE)
        write.table(out_mat, file=outfile, sep="\t", quote=FALSE, row.names=FALSE)
    }
    return(out_mat)
}

## x is a matrix with two columns sampname refname
## refname can be common_normal_01

#' @name create_paired_samplesheet
#' @description create_paired_samplesheet
#' @param x x
#' @param fqmat fqmat
#' @param sampbam sampbam
#' @param refbam refbam
#' @param outfile outfile
#' @param db_sampleid db_sampleid
#' @param db_refid db_refid
#' @param outprefix outprefix
#' @param bampath bampath
#' @param project project
#' @param normal_bam_01 normal_bam_01
#' 
#' @importFrom tools file_path_sans_ext
#' 
#' @export
#' 
create_paired_samplesheet <- function(x, fqmat, 
                                      sampbam, refbam, ## full paths to bam files
                                      outfile, db_sampleid = 0, db_refid = 0, outprefix, 
                                      bampath, ## to be supplied if fastq sheet is supplied
                                      project = 'project',  subproject = 'subproject',
                                      normal_bam_01){
  ## if fqmat is available we will get things from here
  if(!missing(fqmat)){
      tmp <- unlist(sapply(x[,1], function(s) fqmat[fqmat$samplename == s, "recal_bam"][1]))
      if(!is.null(tmp)) tmp_sampbam <- tmp
      if(!missing(bampath)) tmp_sampbam = file.path(bampath, tmp_sampbam)
      refbam <- unlist(sapply(x[,2], function(s){  
        if(s == "common_normal_01") return(normal_bam_01)
        fqmat[fqmat$samplename == s, "recal_bam"][1]
      }))    
      project <- unlist(sapply(x[,1], function(s){  
        paste(fqmat[fqmat$samplename == s, c("project","subproject"), drop = FALSE][1,], collapse="-")}))        
  }else{ ## if no fastq sheet is supplied
    refbam <- ifelse(x[,2] == "common_normal_01", normal_bam_01, refbam)
  }
  if(missing(sampbam)) sampbam <- tmp_sampbam ## if missing them replace from fqmat
  if(missing(outprefix))
    outprefix = sprintf("%s___%s", basename(file_path_sans_ext(sampbam)), basename(file_path_sans_ext(refbam)))
  out_mat = cbind(project = project, samplename = x[,1], refname = x[,2], sampbam = sampbam, refbam = refbam, 
                  db_sampleid = db_sampleid, db_refid = db_refid, outprefix = outprefix)
  write.table(out_mat, file=outfile, sep="\t", quote=FALSE, row.names=FALSE)
  return(out_mat)
}



# read_clone_arch_input_trk()



