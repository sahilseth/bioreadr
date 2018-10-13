## currently this supports mutect, (should work for platypus and other SNP callers toos)
## this can use a direct VCF input as well
## sseth@mda




if(FALSE){

  ## pindel ex
  pindel="/rsrch2/iacs/ngs_runs/150513_SN1120_0334_BC76R7ACXX/pindel/GumbsC-LungCancerMOON-Mulflungcancer-Jay-07-T_C76R7ACXX-1-ATCACG_60x--GumbsC-LungCancerMOON-Mulflungcancer-Jay-09-N_C76R7ACXX-7-TGACCA_60x"
  tumorBAM="/rsrch2/iacs/ngs_runs/150513_SN1120_0334_BC76R7ACXX/bams/GumbsC-LungCancerMOON-Mulflungcancer-Jay-07-T_C76R7ACXX-1-ATCACG_60x.bwa_recalibed.bam"
  normalBAM="/rsrch2/iacs/ngs_runs/150513_SN1120_0334_BC76R7ACXX/bams/GumbsC-LungCancerMOON-Mulflungcancer-Jay-09-N_C76R7ACXX-7-TGACCA_60x.bwa_recalibed.bam"
  annotatePindel(pindel = pindel, tumorBAM = tumorBAM, normalBAM = normalBAM)


}



#' annotate variants
#'
#' annotate variants
#'
#' @param x A data.frame or a file with an extension recognized by read_sheet
#' @param out_path output path
#' @param out_file output file
#' @param execute logical, run annovar dry-run
#' @param annovar_path Path to annovar executable
#' @param db_path path to annovar database
#' @param build build version of the input files
#' @param protocols_g region-type annotation databases, see Annovar's website for details
#' @param protocols_f filter-type annotation databases, see Annovar's website for details
#' @param chr column which represents the chromosome
#' @param start column representing start location
#' @param end column representing end location
#' @param ref column representing reference allele
#' @param alt column representing alternate allele
#'
#' @importFrom stringr str_c
#' @importFrom params read_sheet
#' @import dplyr
#'
#' @details
#'
#' Important default arguments:
#'
#' \code{-otherinfo -nastring . -remove -buildver}
#'
#'
#' @export
annotate_variants <- function(x, outpath, outfile,
                              execute = TRUE,
                              annovar_path = "/scratch/iacs/apps/annovar/latest",
                              db_path = "/scratch/iacs/apps/annovar/humandb",
                              annovar_opts = "-otherinfo -nastring . -remove -buildver",

                              build = "hg19",
                              protocols_g,
                              protocols_f,

                              # specifying input columns:
                              # these defaults work well for mutect
                              chr = "contig", start = "position", end = "position",
                              ref = "ref_allele", alt = "alt_allele",

                              use_uuid = TRUE
){

  if(build == "hg19"){
    if(missing(protocols_g))
      protocols_g = "refGene,knownGene,ensGene"
    if(missing(protocols_f))
      protocols_f = "avsift,ljb26_all,cg69,cosmic72,cosmic72conf,esp6500siv2_all,snp138,nci60,clinvar_20150330,cytoBand,tcga_nb_v1"
  }
  
  # if(length(protocols_f) > 1) protocols_f = paste(protocols_f, collapse = ",")
  protocols = str_c(c(protocols_g, protocols_f), collapse = ",")
  # 	operations = sprintf("%s,%s",
  # 	paste(rep("g", sum(gregexpr(",", protocols_g)[[1]] > 0) + 1), collapse = ","),
  # 	paste(rep("f", sum(gregexpr(",", protocols_f)[[1]] > 0) + 1), collapse = ","))
  op_g = NULL; op_f = NULL
  if(length(protocols_g) > 0)
    op_g = rep("g", sum(gregexpr(",", protocols_g)[[1]] > 0) + 1)
  if(length(protocols_f) > 0)
    op_f = rep("f", sum(gregexpr(",", protocols_f)[[1]] > 0) + 1)
  operations = str_c(c(op_g, op_f), collapse = ",")
  message("Using build: ", build)
  message("Using protocols: ", protocols)
  message("Using operations: ", operations)
  #return()

  ## --- Detecting if x is a data.frame or not
  if(!is.data.frame(x)){
    if(file.exists(x)){
      message("Looks like x is a file. Lets read it...")
      ## --- used for creating output file name
      fl = x
      #x = data.table::fread(x, data.table = FALSE)
      x = read_sheet(x, quote = "")
      if(missing(outpath))
        outpath = dirname(fl)
      # annovar input file
      if(use_uuid)
        input = file.path(outpath, flowr::get_unique_id(tools::file_path_sans_ext(basename(fl)),
                                                        suffix = ".tsv"))
      else
        input = sprintf("%s/%s.tsv", outpath, tools::file_path_sans_ext(basename(fl)))

    }
  }else{
    ## --- if x is DF, we need out_path
    if(missing(outpath))
      stop("Please specify outpath")
    if(missing(outfile))
      outfile = "ann"

    ## --- in outpath create this file
    if(use_uuid)
      input = file.path(outpath, flowr::get_unique_id(outfile, suffix = ".temp"))
    else
      input = sprintf("%s/%s.temp", outpath, tools::file_path_sans_ext(basename(outfile)))
    outfile = file.path(outpath, outfile)
  }
  ## --- outfile is w/o extension (annovar adds multianno extension on its own later)
  if(missing(outfile))
    outfile = sprintf("%s/%s.tsv", outpath, tools::file_path_sans_ext(basename(input)))

    #outfile = file.path(outpath, tools::file_path_sans_ext(basename(input)))

  message("Using ", outfile, " as the final output file")

  ## --- get out_path
  if(!file.exists(outpath))
    dir.create(outpath)

  mutMat = tbl_df(x)
  ## rename the columns so that there are not duplicates
  colnames(mutMat) <- make.names(tolower(colnames(mutMat)),
                                 unique=TRUE, allow_ = TRUE)
  ## --- check if all cols are present
  cols = c(chr, start, end, ref, alt)
  if(sum(!cols %in% colnames(mutMat)) > 0)
    stop("Some columns are missing: ",
         paste(cols[!cols %in% colnames(mutMat)], collapse = " "),
         "\nFew columns of matrix are: ",
         paste(colnames(mutMat)[1:10], collapse = " "))
  ## --- rename columns if necessary
  ## --- Adding a key to mutmat to be able to merge data later...
  ## --- use mutate instead of rename, to make sure all is good !
  #mutMat$annovarkey = paste("key", 1:nrow(mutMat), sep = "-")
  mutMat <- mutate_(mutMat,
                    chr = chr, start = start, end = end, ref = ref, alt = alt) %>%
    mutate(key=paste(chr, start, end,ref,alt, sep=":"))
  message("Fixing chr names, X-->22, Y-->23") ## annovar compatibility
  mutMat <- fix_chrnames(mutMat, build)
  mutMat2 = dplyr::select(mutMat,
                          chr, start, end, ref, alt, key) %>% unique()
  message("Working on: ", nrow(mutMat2), " unique variations")
  #on.exit(file.remove(input))
  write.table(mutMat2, file = input, sep = "\t",
              col.names = FALSE, row.names = FALSE, quote = FALSE)

  #outfile = flowr:::get_unique_id(tools::file_path_sans_ext(basename(fl)))
  cmd = paste("perl", file.path(annovar_path, "table_annovar.pl"),
              ## input, db_path, get other columsns, put . for NA, remove tmp files, build
              input, file.path(db_path), annovar_opts, build,
              "-protocol", protocols, "-operation", operations,
              "-out", outfile)
  message(cmd)
  if(execute) system(cmd)

  ann_out = sprintf("%s.%s_multianno.txt", outfile, build)
  #final_out = sprintf("%s.tsv", outfile)
  final_out = outfile
  message("Merging with input file")
  ann = tbl_df(read_sheet(ann_out, quote = ""))
  #ann = fread(outfile)
  #print(head(ann))
  #print(head(mutMat))
  ann = rename_(ann, key = 'Otherinfo')
  colnames(ann) = tolower(colnames(ann))
  ann2 = merge(mutMat, ann, by = c("chr", "start", "end", "ref", "alt", "key"))
  message("Writing out final output: ", final_out)
  write.table(ann2, file = final_out, sep = "\t", quote = FALSE, row.names = FALSE)
  warnings()

  ## remove temp files
  unlink(c(ann_out, input))

  invisible(list(cmd = cmd, outfile = final_out, data = ann2))
}





fix_chrnames <- function(mutMat, build){
  mutMat$chr = gsub("X", "23", mutMat$chr)
  mutMat$chr = gsub("Y", "24", mutMat$chr)
  return(mutMat)
}







#' annotate_mutect
#'
#' annotate_mutect
#'
#' @param x a mutect file
#'
#' @export
#'
annotate_mutect <- function(x, odir, oprefix,
                            tumorBAM, normalBAM, ## optional

                            ## --- annovar defaults
                            organism = "human", build,
                            annovar_path = "/scratch/iacs/apps/annovar/latest",

                            ## annovar's additional options.
                            protocols_g, protocols_f,

                            ## annotate ALL the variants (even those rejected by mutect)
                            full = FALSE,

                            ## --- input file details, fixed for mutect
                            chrom = "contig", start = "position",
                            end = "position", ref = "ref_allele", alt = "alt_allele",

                            use_uuid = TRUE){

  require(parallel) ## used is split by chrs
  ## -- get missing values if missing
  if(missing(odir)){
    odir <- dirname(x)
  }
  if(missing(oprefix)){
    #oprefix <- gsub("\\.txt$", "", basename(x))
    oprefix = tools::file_path_sans_ext(basename(x))
  }
  if(missing(tumorBAM) | missing(normalBAM)){
    message("bams are missing, skipping recounts")
    reCount=FALSE
  }

  ## --- setting defauls for human
  if(organism == "human"){
    if(missing(build)){
      build = "hg19"; message("Using build ", build)}
    if(missing(protocols_g)) protocols_g = "refGene,knownGene"
    if(missing(protocols_f)) protocols_f = "avsift,ljb26_all,cg69,cosmic70,cosmic72,esp6500siv2_all,snp138,nci60,clinvar_20150330,cytoBand,tcga_nb_v1"
  }

  ## --- setting defaults for mouse
  if(organism == "mouse"){
    if(missing(build)) build = "mm10"
    if(missing(protocols_g)) protocols_g = "refGene,knownGene"
    if(missing(protocols_f)) protocols_f = "avsift,ljb26_all,cg69,cosmic70,cosmic72,esp6500siv2_all,snp138,nci60,clinvar_20150330,cytoBand,tcga_nb_v1"
  }


  if(FALSE){
    dbtype = match.arg(dbtype)
    if(dbtype == "refgene"){
      filters = "refseq_mrna"
      att = c("uniprot_swissprot_accession", "interpro", "refseq_peptide", "ucsc")
    }else{
      filters = "ucsc"
      att = c("uniprot_swissprot_accession", "interpro", "refseq_peptide", "refseq_mrna")
    }
  }

  tab = readMutect(x, full = full)
  ## --- ANNOTATE the file
  message("annotate_mutect, using oprefix: ", oprefix)
  tmp <- annotate_variants(tab,
                           outpath = odir,
                           outfile = oprefix,
                           annovar_path = annovar_path,
                           build = build, use_uuid = use_uuid,
                           protocols_g = protocols_g, protocols_f = protocols_f,
                           ## --- details of how to use the columns
                           chr = "contig", start = "position", end = "position",
                           ref = "ref_allele",alt= "alt_allele")


  ann = tmp$data
  annfile <- tmp$outfile

  write_sheet(ann, file = annfile, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  return(annfile)
}

#' Not tested
#'
#' @export
annotate_vcf <- function(x, oprefix, odir,
                         build, use_uuid,
                         annovar_path, protocols_f, protocols_g, ...){

  tab = parse_vcf(x)

  if(missing(odir))
    odir = dirname(x)
  message("annotate_mutect, using oprefix: ", oprefix)
  tmp <- annotate_variants(tab,
                           outfile = oprefix,
                           outpath = odir,
                           annovar_path = annovar_path,
                           build = build, use_uuid = use_uuid,
                           protocols_g = protocols_g, protocols_f = protocols_f,
                           ## --- details of how to use the columns
                           chr = "chrom", start = "pos", end = "pos",
                           ref = "ref", alt= "alt")
}



annotate_pindel <- function(){

}


#' Generates a cmdline, when invoked annotates files using annovar
#'
#' @param x files to be annotates
#' @param samplename name of the sample, to be added to the flowmat
#' @param outfile name of the final merged file.
#'
#' @export
annotate <- function(files, samplename = opts_flow$get("samplename"),
                     outfile,
                     build = "hg19",
                     annovar_dir = opts_flow$get("annovar_dir"),
                     annotate_func = "ultraseq::annotate_mutect"){

  if(missing(outfile))
    outfile <- paste0(tools::file_path_sans_ext(basename(files[1])),  "_merged.annovar.tsv")

  check_args(ignore = "outfile")

  # output file names
  anns = paste0(tools::file_path_sans_ext(basename(files)), "_annovar.txt")

  cmd_ann = sprintf("flowr %s x=%s oprefix=%s build=%s full=TRUE use_uuid=FALSE annovar_path=%s",
                    annotate_func, x, anns, build, annovar_dir)

  cmd_merge = sprintf("flowr ultraseq::merge_sheets x=%s outfile=%s",
                      paste(anns, collapse = ","), outfile)
  cmds = list(annotate = cmd_ann, annotate_merge = cmd_merge)

  flowmat = to_flowmat(cmds, samplename = samplename)

  return(list(flowmat = flowmat, outfiles = list(merged = outfile)))
}
attr(annotate, "type") <-"module"


## --- testing annotate_mutect
if(FALSE){

  debug(annotate_mutect)
  mutect = "/rsrch1/iacs/ngs_runs/Sarcomatoid/mutect/WISTUBA-lung_sarco-WEX-973-T__WISTUBA-lung_sarco-WEX-973-N_merged.muTect_call_stats.txt"
  tmp <- annotate_mutect(mutect)

  x="/rsrch2/iacs/flowr/runs/sarco/rna_mutect/rna_mutect-610-20150915-04-08-15-V41B7wF2/tmp/RNA_Barcode_None_001.R_2014_08_22_14_58_37_user_PRO-65-Sarcomatoid_RNA_Seq_610-T2_trimmed_trimmed.fastq_rg_reorder_chr20.mutect.txt"
  require(stringr);require(dplyr)
  #source('~/iacsSVN/RPacks/IACSUtil/R/annotate_variants.R')
  #devtools::load_all("~/Dropbox/iacsSVN/RPacks/IACSUtil")
  #debug(annotate_mutect)
  #debug(formatMutect)
  annotate_mutect(x, protocols_f = "cosmic70")


}
