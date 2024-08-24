
# GOAL:
# read CNV data
# conv to GR; then to gene level
# GENExSAMPLE matrix
# ragged experiment


if(FALSE){
  
  library(pacman)
  p_load(tidyverse, janitor, glue)
  p_load(cowplot, egg)
  
  # this will be wrapped in a function
  gencode_fl = "~/rsrch2_home/ref/human/b37/annotations/gencode/v19/gencode.v19.annotation_gene.bed"
  col_sample = "samplename"
  col_chr = "chrom"
  col_start = "start"
  col_end = "end"
  col_num_mark = "num_mark"
  col_seg_mean = "cnlr_median"
  
  expand_using = "genes"
  # ** test one ------
  odir = "~/projects2/ss_tnbc/analysis/art/wex/2019Q1_b2/facets"
  df_seg = read_rds("~/projects2/ss_tnbc/analysis/art/wex/2019Q1_b2/facets/df_seg.rds")
  df_facet_summ = read_rds("~/projects2/ss_tnbc/analysis/art/wex/2019Q1_b2/facets/df_facet_summ.rds")
  df_seg = left_join(df_seg, df_facet_summ, by = c("samplename" = "t_path_sample_id"))
  # get chr lens
  source('~/projects/packs_dnaseq/R/plot_heat_cn_v1.R')
  df_chr_len = get_hg19_chr_lens()
  # get data heatmap
  
  # ** exomecn seg --------
  df_seg = read_rds("~/projects2/ss_tnbc/analysis/art/wex/2019Q1_b2/facets/df_seg_exomecn.rds");df_seg
  # to igv style & annotate with genes
  source('~/Dropbox/projects/packs_dnaseq/R/to_mae.seg.R')
  lst_seg = to_igv_seg(df_seg,
                       col_sample = "samplename",
                       col_start = "loc_start",
                       col_end = "loc_end",
                       col_seg = "seg_mean")
  
  lst_seg_ann = annotate_cnv.seg(lst_seg$gr_seg)
  
  
  # create a ht of all 20K genes!
  source('~/projects/packs_dnaseq/R/to_mae.seg.R')
  lst_seg_ann = annotate_cnv.seg(lst_seg$gr_seg)
  
  chr_splits = gsub("chr", "", rowData(lst_seg_ann$se_cnv)$chrom)
  chr_splits = factor(chr_splits, levels = c(1:23, "X", "Y"))
  p_load(ComplexHeatmap, SummarizedExperiment)
  packageVersion("ComplexHeatmap")
  ht = Heatmap(assay(lst_seg_ann$se_cnv),
               name = "CNV",
               row_gap = unit(0.5, "mm"),
               row_split = gsub("chr", "", rowData(lst_seg_ann$se_cnv)$chrom),
               raster_device = "CairoPNG", column_names_gp = gpar(cex = 0.3),
               use_raster = T, cluster_rows = F, show_row_names = F, show_row_dend = F)
  # plot_size(12, 8)
  wd = "~/projects2/ss_tnbc/analysis/art/wex/2019Q1_b2"
  setwd(wd);odir="facets";wranglr::mkdir(odir)
  pdf(glue("{odir}/p_ht_exomecn.pdf"), width = 8, height = 7)
  draw(ht)
  dev.off()
  
  
  # ** test gatkcnv sarco -------
  df_seg = read_tsv('/rsrch3/home/iacs/sseth/flows/SS/sarco/mda/wex/gatkcnv/merged.matched.called.igv.tsv') %>%
    clean_names()
  str(df_seg)
  # debug(to_igv_seg)
  lst_seg = to_igv_seg(df_seg,
                       col_sample = "sample",
                       col_start = "start",
                       col_end = "end",
                       col_seg = "segment_mean",
                       col_chr = "chromosome",
                       col_num_mark = "num_probes")
  
  # create a ht of all 20K genes!
  source('~/projects/packs_dnaseq/R/to_mae.seg.R')
  lst_seg_ann = annotate_cnv.seg(lst_seg$gr_seg)
  
  p_load(ComplexHeatmap, MultiAssayExperiment, SummarizedExperiment)
  chr_splits = gsub("chr", "", rowData(lst_seg_ann$se_cnv)$chrom)
  chr_splits = factor(chr_splits, levels = c(1:23, "X", "Y"))
  p_load(ComplexHeatmap, SummarizedExperiment)
  packageVersion("ComplexHeatmap")
  ht = Heatmap(assay(lst_seg_ann$se_cnv),
               name = "CNV",
               row_gap = unit(0.5, "mm"),
               row_split = gsub("chr", "", rowData(lst_seg_ann$se_cnv)$chrom),
               raster_device = "CairoPNG", column_names_gp = gpar(cex = 0.3),
               use_raster = T, cluster_rows = F, show_row_names = F, show_row_dend = F)
  ht
  # plot_size(12, 8)
  wd = "~/projects2/ss_tnbc/analysis/art/wex/2019Q1_b2"
  setwd(wd);odir="facets";wranglr::mkdir(odir)
  pdf(glue("{odir}/p_ht_exomecn.pdf"), width = 8, height = 7)
  draw(ht)
  dev.off()
  
  
  
  
  
}




# https://software.broadinstitute.org/software/igv/SEG
# The segmented data file format is the output of the Circular Binary Segmentation algorithm (Olshen et al., 2004).
# 'ID	chrom	loc.start	loc.end	num.mark	seg.mean
# GenomeWideSNP_416532	1	51598	76187	14	-0.7116
# GenomeWideSNP_416532	1	76204	16022502	8510	-0.029
# GenomeWideSNP_416532	1	16026084	16026512	6	-2.0424
# GenomeWideSNP_416532	1	16026788	17063449	424	-0.1024
to_gr_seg.seg <- function(df_seg,
                          col_sample = "samplename",
                          col_chr = "chrom",
                          col_start = "start",
                          col_end = "end",
                          col_num_mark = "num_mark",
                          col_seg_mean = "cnlr_median",
                          
                          save_fls = FALSE
                          
                          
){
  p_load(futile.logger)
  source('~/Dropbox/public/github_wranglr/R/expect_columns.R')
  col_igv_default = c(col_sample, col_chr,
                      col_start, col_end,
                      col_num_mark, col_seg_mean)
  # check column names
  expect_columns(df_seg, columns = col_igv_default)
  
  col_others = setdiff(colnames(df_seg), col_igv_default)
  col_igv_new = c(".sample_id", "chrom", "start", "end", "num_mark", col_others, col_seg_mean)
  
  flog.debug("create a IGV table")
  df_igv_seg <- df_seg[, c(col_sample, col_chr, col_start, col_end, col_num_mark,
                           # other columns
                           col_others,
                           # SEG is LAST column, and should be last
                           col_seg_mean)]
  colnames(df_igv_seg) = col_igv_new
  
  if(save_fls){
    write_tsv(df_igv_seg, glue("{odir}/df_igv_seg.seg"))
    write_rds(df_igv_seg, glue("{odir}/df_igv_seg.rds"))
  }
  
  flog.debug("create genomic ranges object")
  gr_seg = data.frame(df_igv_seg) %>%
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = T)
  
  # class(df_igv_seg) <- c("df_igv_seg", "as_tibble", "tbl", "data.frame")
  list(df_igv_seg = df_igv_seg,
       gr_seg = gr_seg)
}
to_igv.seg = to_gr_seg.seg





#' to_long.igv_seg
#'
#' convert the IGV segments to genes (expanded version)
#'
#' @param gr_seg
#' @param gencode_fl
#' @param expand_using
#'
#' @export
to_long.gr_seg <- function(gr_seg,
                           # gene positions for hg19
                           expand_using = c("marks", "genes")){
  
  expand_using = match.arg(expand_using)
  
  # testit::assert("input is of class df_igv_seg", {
  #   class(df_igv_seg)[1] ==  "df_igv_seg"
  # })
  
  # gr_seg
  # message("annotate with genes")
  
  
  message("conv to bins for clustering")
  df_igv_seg_mrk = expand_bins(df_igv_seg)
}


to_gr.gencode <- function(gencode_fl = "~/.rsrch3/home/iacs/sseth/ref/human/b37/annotations/gencode/v19/gencode.v19.annotation_gene_hgnc.bed"){
  pacman::p_load(GenomicRanges, dplyr, janitor, glue, readr)
  df_gencode = data.table::fread(gencode_fl, data.table = F) %>% as_tibble() %>% clean_names()
  message("keep PC genes from gencode...")
  df_gencode_pc = dplyr::filter(df_gencode, gene_type == "protein_coding");dim(df_gencode_pc)
  df_gencode_pc_f1 = dplyr::select(df_gencode_pc, number_seqnames:strand, gene_id, gene_name) %>%
    group_by(number_seqnames, start, end, width, strand, gene_name) %>%
    summarise(gene_id = paste(gene_id, collapse = ";"))
  df_gencode_pc_f1 %>% dplyr::filter(grepl(";", gene_id))
  # cases where we have same start, end, gene_name, merge gene_id
  gr_genes = GenomicRanges::makeGRangesFromDataFrame(df_gencode_pc,
                                                     seqnames.field = "number_seqnames",
                                                     keep.extra.columns = T,
                                                     starts.in.df.are.0based = FALSE)
  seqlevelsStyle(gr_genes) <- "UCSC"
  
  list(df_gencode = df_gencode,
       df_gencode_pc = df_gencode,
       gr_genes = gr_genes)
  
}

to_gr.ens91 <- function(ens91_fl = "~/stash_home/ref/annotations/annotation_hub/ensembl91_hg38_AH60773.rds"){
  pacman::p_load(GenomicRanges, dplyr, janitor, glue, readr)

  ens91 = read_rds(ens91_fl)
  # we will keep all genes, and look at protein coding genes, only later
  # there might be issue of same region represented twice, but will figure that out later
  # cases where we have same start, end, gene_name, merge gene_id
  gr_genes = GenomicRanges::makeGRangesFromDataFrame(ens91$df_ens)

  seqlevelsStyle(gr_genes) <- "UCSC"
  
  list(df_gencode = df_gencode,
       df_gencode_pc = df_gencode,
       gr_genes = gr_genes)
  
}

if(FALSE){
  gr_seg = read_rds("/stash/results/dev/seths3/P04492-20220524-0006-fl-rel/data/ngs360/wes/cnv/merged/sclust/gr_seg.rds")
  do.call( )
  gr_genes = ens91$hs.genes

  df_cnv_ann_lng %>% 
    write_rds("/stash/results/dev/seths3/P04492-20220524-0006-fl-rel/data/ngs360/wes/cnv/merged/sclust/df_cnv_ann_lng.rds")


}


# annotate CNV segment file
# also create a matrix of genes X samples
# https://support.bioconductor.org/p/97647/
#' annotate_cnv.seg
#'
#' incase of multiple values per gene, get the more extreme value
#'
#' @param gr_seg
#' @param gencode_fl
#'
#' @export
annotate.gr_seg <- function(gr_seg,
                            gr_genes = lst_gen$gr_genes
){
  # print(class(gr_seg) == "GRanges")
  p_load(assertthat)
  assert_that(class(gr_seg) == "GRanges", msg = "gr_seg is not of class GRanges")
  assert_that(class(gr_genes) == "GRanges", msg = "gr_genes is not of class GRanges")
  
  library(GenomicRanges)
  pacman::p_load(GenomicRanges, dplyr, janitor, glue, readr, futile.logger)
  
  # gr_genes = lst_gen$gr_genes
  
  # length(gr_genes)
  flog.debug("make sure both use NCBI style")
  genomeStyles("Homo sapiens") %>% head(2)
  
  seqlevelsStyle(gr_seg) <- "NCBI"
  seqlevelsStyle(gr_genes) <- "NCBI"
  
  GenomicRanges::seqnames(gr_seg)
  GenomicRanges::seqnames(gr_genes)
  
  flog.debug("intersect seg & genes")
  # moe on this here: https://support.bioconductor.org/p/67118/#67148
  # ann with only PC genes (for matrix) - easier to understand
  # for gene with no PC regions what should be done??
  # this is whole-exon data (so having data on PC regions is a given!!)
  olaps = GenomicRanges::findOverlaps(gr_seg, gr_genes)
  long_annotated = gr_seg[queryHits(olaps)]
  # there would certainly be multiple gene per seg
  long_annotated$gene_name = gr_genes[subjectHits(olaps)]$gene_name
  long_annotated$gene_id = gr_genes[subjectHits(olaps)]$gene_id
  long_annotated$gene_biotype = gr_genes[subjectHits(olaps)]$gene_biotype
  # df_ann = ann[subjectHits(olaps)]
  # df_cnv_ann = cbind(long_annotated, df_ann)
  
  df_cnv_ann_lng = as.data.frame(long_annotated) %>% as_tibble()
  df_cnv_ann_lng
}


to_se.ann_seg.exomecn <- function(df_cnv_ann_lng, lst_gen){
  p_load(futile.logger)
  
  message("check dups")
  # --------+++++++++------------ (seg mean)
  # ----GGGGGGGGG---------------- (gene annotation)
  # gene has two different values (which one to use?)
  # tmp = dplyr::count(df_cnv_ann_lng, .sample_id, gene_name) %>% filter(n>1)
  # head(tmp)
  
  message("convert into a matrix")
  # two segments intersect with the SAME gene!
  # filter(df_cnv_ann_lng, gene_name %in% tmp$gene_name[2], .sample_id %in% tmp$.sample_id[2])
  # filter(df_gencode, gene_name %in% tmp$gene_id[1])
  # idea:
  get_extreme <- function(x){
    if(length(x) == 1) return(x)
    # mean(x, na.rm = T)
    # if(is.na(x)) return(x)
    x[which(abs(x) == max(abs(x)))][1]
  }
  # x = c(0.2, 0.2);get_extreme(x)
  # x = c(0.2, 1.3);get_extreme(x)
  # x = c(0.2, -1.3);get_extreme(x)
  
  # df_tmp = sample_frac(df_cnv_ann_lng, size = 0.1);dim(df_tmp)
  mat_cnv_ann_wd = reshape2::dcast(df_cnv_ann_lng,
                                   gene_name ~ .sample_id,
                                   value.var = "seg_mean",
                                   fun.aggregate = get_extreme,
                                   fill = -100)
  dim(mat_cnv_ann_wd)
  # head(df_cnv_ann_wd)
  p_load(magrittr)
  mat_cnv_ann_wd %<>% wranglr::to_mat()
  # replace -100 with NA
  mat_cnv_ann_wd[mat_cnv_ann_wd == -100] <- NA
  
  flog.info("prepare rowdata")
  rowdata = filter(lst_gen$df_gencode_pc,
                   gene_name %in% rownames(mat_cnv_ann_wd)) %>%
    mutate(index = 1:n()) %>%
    group_by(gene_name) %>% top_n(1, index) %>%
    ungroup() %>%
    # getting chr from gencode
    dplyr::rename(chrom = number_seqnames) %>%
    mutate(chrom = factor(chrom, levels = seqlevels(lst_gen$gr_genes))) %>%
    data.frame(row.names = .$gene_name)
  # head(rowdata)
  rowdata %<>% arrange(chrom, start, end)
  
  se_cnv <- SummarizedExperiment::SummarizedExperiment(
    list(cnv_log2 = mat_cnv_ann_wd[rowdata$gene_name, ]),
    rowData = rowdata)
  se_cnv
  
  list(mat_cnv_ann_wd = mat_cnv_ann_wd,
       df_cnv_ann_lng = df_cnv_ann_lng,
       se_cnv = se_cnv)
}

# conv to integer CNV
# https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/CNV_Pipeline/
# Genes with focal CNV values smaller than -0.3 are categorized as a "loss" (-1)
# Genes with focal CNV values larger than 0.3 are categorized as a "gain" (+1)
# Genes with focal CNV values between and including -0.3 and 0.3 are categorized as "neutral" (0).


to_mae.titan <- function(df_ann, trk, gencode_fl){
  
  # ** create matrix -----
  message("convert into a matrix")
  resolve_state <- function(x){
    if(length(x) == 1) return(x)
    # mean(x, na.rm = T)
    # if(is.na(x)) return(x)
    x[which(abs(x) == max(abs(x)))][1]
  }
  # x=c("NEUT", "HETD")
  resolve_corrected_call <- function(x){
    if(length(x) == 1) return(x)
    df_titan_ref = titan_cna_states()
    indxs = which(df_titan_ref$titan_corrected_call2 %in% x)
    df_titan_ref$titan_corrected_call2[max(indxs)]
  }
  # x=c("DLOH", "DHET")
  resolve_call <- function(x){
    if(length(x) == 1) return(x)
    df_titan_ref = titan_cna_states()
    indxs = which(df_titan_ref$titan_call %in% x)
    df_titan_ref$titan_state[max(indxs)]
  }
  
  mat_ann_state_wd = df_ann %>%
    pivot_wider(id_cols = "gene_name", names_from = ".sample_id", values_from = "titan_state",
                values_fn = resolve_state) %>% to_mat()
  mat_ann_call_wd = df_ann %>%
    pivot_wider(id_cols = "gene_name", names_from = ".sample_id", values_from = "titan_call",
                values_fn = resolve_call) %>% to_mat()
  mat_ann_corrected_call_wd = df_ann %>%
    pivot_wider(id_cols = "gene_name", names_from = ".sample_id",
                values_from = "corrected_call", values_fn = resolve_corrected_call) %>%
    to_mat()
  
  lst_gencode = to_grseg.gencode(gencode_fl = gencode_fl)
  
  # prepare rowdata
  gr_genes - lst2$gr_seg[[1]]
  rowdata = filter(lst_gencode$df_gencode, gene_name %in% rownames(mat_ann_corrected_call_wd)) %>%
    mutate(index = 1:n()) %>%
    group_by(gene_name) %>% top_n(1, index) %>%
    ungroup() %>%
    # getting chr from gencode
    dplyr::rename(chrom = number_seqnames) %>%
    mutate(chrom = factor(chrom, levels = seqlevels(lst_gencode$gr_genes))) %>%
    data.frame(row.names = .$gene_name)
  # head(rowdata)
  rowdata %<>% arrange(chrom, start, end)
  
  se_corrected_call <- SummarizedExperiment::SummarizedExperiment(mat_ann_corrected_call_wd[rowdata$gene_name, ],
                                                                  rowData = rowdata)
  se_call <- SummarizedExperiment::SummarizedExperiment(mat_ann_call_wd[rowdata$gene_name, ],
                                                        rowData = rowdata)
  se_state <- SummarizedExperiment::SummarizedExperiment(mat_ann_state_wd[rowdata$gene_name, ],
                                                         rowData = rowdata)
  
  p_load(MultiAssayExperiment)
  coldata = trk %>% data.frame(row.names = .$name)
  mae = MultiAssayExperiment(experiments = list(se_corrected_call = se_corrected_call,
                                                se_call = se_call,
                                                se_state = se_state), colData = coldata)
  
  mae
  
}



to_mae.sclust <- function(df_ann_icn, df_ann_as, trk, gr_genes){
  
  p_load(MultiAssayExperiment, futile.logger)
  
  # ** create matrix -----
  # get the abs max, in case of duplicate values
  resolve_state <- function(x){
    if(length(x) == 1) return(x)
    # mean(x, na.rm = T)
    # if(is.na(x)) return(x)
    x[which(abs(x) == max(abs(x)))][1]
  }
  
  flog.info("> convert into a matrix")
  flog.info(">> resolve multiple states per gene, take abs max for duplicate values (un corr)")
  mat_cn_corr_wd = df_ann_icn %>%
    pivot_wider(id_cols = "gene_name", names_from = "bms_sample_id", values_from = "cnv_corrected",
                values_fn = resolve_state) %>% wranglr::to_mat()
  mat_cn_corr_wd[1:5, 1:5] 
  
  flog.info(">> resolve multiple states per gene, take abs max for duplicate values (corr)")
  mat_cnv_uncorr_wd = df_ann_icn %>%
    pivot_wider(id_cols = "gene_name", names_from = "bms_sample_id", values_from = "log2_copy_ratio",
                values_fn = resolve_state) %>% wranglr::to_mat()
  dim(mat_cnv_uncorr_wd)
  
  flog.info(">> resolve multiple states per gene, take abs max for duplicate values (corr)")
  mat_copy_nr_wd = df_ann_as %>%
    pivot_wider(id_cols = "gene_name", names_from = "bms_sample_id", values_from = "copy_nr",
                values_fn = resolve_state) %>% wranglr::to_mat()
  mat_copy_nr_wd[1:2,1:2];dim(mat_copy_nr_wd)
  
  flog.info("create rowdata/rowranges")
  # lst_gencode = to_grseg.gencode(gencode_fl = gencode_fl)
  rowdata = gr_genes[rownames(mat_cn_corr_wd), ]
  
  flog.info("create summrizedexperiment")
  se_cn_corr <- SummarizedExperiment::SummarizedExperiment(mat_cn_corr_wd,
                                                         rowRanges = rowdata)
  se_uncorr <- SummarizedExperiment::SummarizedExperiment(mat_cnv_uncorr_wd,
                                                              rowRanges = rowdata)
  se_copy_nr <- SummarizedExperiment::SummarizedExperiment(mat_copy_nr_wd,
                                                          rowRanges = gr_genes[rownames(mat_copy_nr_wd), ])
  
  flog.info("create coldata")
  coldata = trk %>% data.frame(row.names = .$name)
  
  flog.info("create MAE")
  mae = MultiAssayExperiment(experiments = list(se_cn_corr = se_cn_corr,
                                                se_uncorr = se_uncorr,
                                                se_copy_nr = se_copy_nr), colData = coldata)
  mae
}

to_mae.gatk <- function(df_ann, trk, gr_genes){
  
  p_load(MultiAssayExperiment, futile.logger)
  
  # ** create matrix -----
  # get the abs max, in case of duplicate values
  resolve_state <- function(x){
    if(length(x) == 1) return(x)
    # mean(x, na.rm = T)
    # if(is.na(x)) return(x)
    x[which(abs(x) == max(abs(x)))][1]
  }
  
  flog.info("> convert into a matrix")
  flog.info(">> resolve multiple logr per gene, take abs max for duplicate values")
  mat_logr_med_wd = df_ann %>%
    pivot_wider(id_cols = "gene_name", names_from = "sample_id", values_from = "log2_copy_ratio_posterior_50",
                values_fn = resolve_state) %>% wranglr::to_mat()
  dim(mat_logr_med_wd);mat_logr_med_wd[1:5, 1:5] 
  
  mat_segmean_wd = df_ann %>%
    pivot_wider(id_cols = "gene_name", names_from = "sample_id", values_from = "mean_log2_copy_ratio",
                values_fn = resolve_state) %>% wranglr::to_mat()
  dim(mat_segmean_wd);mat_segmean_wd[1:5, 1:5] 

  flog.info(">> resolve multiple baf per gene, take abs max for duplicate values (corr)")
  mat_baf_wd = df_ann %>%
    pivot_wider(id_cols = "gene_name", names_from = "sample_id", values_from = "minor_allele_fraction_posterior_50",
                values_fn = resolve_state) %>% wranglr::to_mat()
  dim(mat_baf_wd);mat_baf_wd[1:5, 1:5] 
  
  flog.info("create rowdata/rowranges")
  # lst_gencode = to_grseg.gencode(gencode_fl = gencode_fl)
  rowdata = gr_genes[rownames(mat_logr_med_wd), ]
  
  flog.info("create summrizedexperiment")
  se_logr_med <- SummarizedExperiment::SummarizedExperiment(mat_logr_med_wd,
                                                        rowRanges = rowdata)
  se_segmean <- SummarizedExperiment::SummarizedExperiment(mat_segmean_wd,
                                                            rowRanges = rowdata)
  se_baf <- SummarizedExperiment::SummarizedExperiment(mat_baf_wd,
                                                       rowRanges = rowdata)
  flog.info("create coldata")
  coldata = trk %>% data.frame(row.names = .$name)
  
  flog.info("create MAE")
  mae = MultiAssayExperiment(experiments = list(se_logr_med = se_logr_med,
                                                se_segmean = se_segmean,
                                                se_baf = se_baf), 
                             colData = coldata)
  mae
}

to_cnvranger <- function(gr_seg){
  # https://bioconductor.org/packages/devel/bioc/vignettes/CNVRanger/inst/doc/CNVRanger.html#input-data-format
  # create a list
  grl = split(gr_seg, gr_seg$.sample_id)
  
  ra = RaggedExperiment(grl)
  ra
  
  # assay(ra[1:5,1:5])
  # library(CNVRanger)
  # # gistic like
  # cnvrs <- populationRanges(grl, density=0.1, est.recur=TRUE)
  # cnvrs
  
  p_load(AnnotationHub)
  ah <- AnnotationHub()
  
  ahDb <- query(ah, pattern = c("Homo sapiens", "EnsDb"))
  ahDb
  
}





# COPIED from: https://rdrr.io/bioc/maftools/src/R/readSegs.R
#--- Change segment sizes into linear scale
transformSegments = function(segmentedData, build = 'hg19'){
  
  build.opts = c('hg19', 'hg18', 'hg38')
  
  if(!build %in% build.opts){
    stop('Available reference builds: hg18, hg19, hg38')
  }
  
  if(build == 'hg19'){
    chr.lens = c(249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663,
                 146364022, 141213431, 135534747, 135006516, 133851895, 115169878, 107349540,
                 102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566,
                 155270560, 59373566)
  } else if(build == 'hg18'){
    chr.lens = c(247249719, 242951149, 199501827, 191273063, 180857866, 170899992,
                 158821424, 146274826, 140273252, 135374737, 134452384, 132349534,
                 114142980, 106368585, 100338915, 88827254, 78774742, 76117153,
                 63811651, 62435964, 46944323, 49691432, 154913754, 57772954)
  } else if(build == 'hg38'){ #hg38
    chr.lens = c(248956422, 242193529, 198295559, 190214555, 181538259, 170805979,
                 159345973, 145138636, 138394717, 133797422, 135086622, 133275309,
                 114364328, 107043718, 101991189, 90338345, 83257441, 80373285,
                 58617616, 64444167, 46709983, 50818468, 156040895, 57227415)
  } else{
    stop('Available reference builds: hg18, hg19, hg38')
  }
  
  segmentedData[,Start_Position := as.numeric(as.character(Start_Position))]
  segmentedData[,End_Position := as.numeric(as.character(End_Position))]
  
  #Replace chr x and y with numeric value (23 and 24) for better ordering
  segmentedData$Chromosome = gsub(pattern = 'chr', replacement = '', x = segmentedData$Chromosome, fixed = TRUE)
  segmentedData$Chromosome = gsub(pattern = 'X', replacement = '23', x = segmentedData$Chromosome, fixed = TRUE)
  segmentedData$Chromosome = gsub(pattern = 'Y', replacement = '24', x = segmentedData$Chromosome, fixed = TRUE)
  
  segmentedData$Chromosome = factor(x = segmentedData$Chromosome, levels = 1:24, labels = 1:24)
  
  segmentedData = segmentedData[order(Chromosome, Start_Position, decreasing = FALSE)]
  
  seg.spl = split(segmentedData, segmentedData$Chromosome)
  
  seg.spl.transformed = seg.spl[[1]]
  if(nrow(seg.spl.transformed) > 0){
    seg.spl.transformed$Start_Position_updated = seg.spl.transformed$Start_Position
    seg.spl.transformed$End_Position_updated = seg.spl.transformed$End_Position
  }
  
  chr.lens.sumsum = cumsum(chr.lens)
  
  for(i in 2:length(seg.spl)){
    
    x.seg = seg.spl[[i]]
    if(nrow(x.seg) > 0){
      x.seg$Start_Position_updated = x.seg$Start_Position + chr.lens.sumsum[i-1]
      x.seg$End_Position_updated = x.seg$End_Position + chr.lens.sumsum[i-1]
    }
    seg.spl.transformed = rbind(seg.spl.transformed, x.seg, fill = TRUE)
  }
  
  return(seg.spl.transformed)
}




# END


