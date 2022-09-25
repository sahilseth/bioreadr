
# singularity pull docker://quay.io/bcbio/bcbio-vc


# https://github.com/gavinha/TitanCNA/tree/master/scripts/R_scripts
# numClusters=3
# numCores=4
# for ploidy in [2, 3, 4]:
# for num_clusters in [1, 2, 3]:
# ## run TITAN for each ploidy (2,3,4) and clusters (1 to numClusters)
# echo "Maximum number of clusters: $numClusters";
# for ploidy in $(seq 2 4)
# do
# echo "Running TITAN for $i clusters.";
# outDir=run_ploidy$ploidy
# mkdir $outDir
# for numClust in $(seq 1 $numClusters)
# do
# echo "Running for ploidy=$ploidy";
# Rscript titanCNA_v1.10.1.R --id test --hetFile test.het.txt --cnFile test.cn.txt \
# --numClusters $numClust --numCores $numCores --normal_0 0.5 --ploidy_0 $ploidy \
# --chrs "c(1:22, \"X\")" --estimatePloidy TRUE --outDir $outDir
# done
# echo "Completed job for $numClust clusters."
# done
# 
# ## select optimal solution
# Rscript selectSolution.R --ploidyRun2=run_ploidy2 --ploidyRun3=run_ploidy3 --ploidyRun4=run_ploidy4 --threshold=0.05 --outFile optimalClusters.txt

.run_titan_example <- function(){
  
  # test trk file:
  trk = "~/.rsrch3/scratch/iacs/sseth/flows/SS/tnbc/ms51_wex_b1/superfreq/df_sfq_metadata_ms51_dna_b1.tsv"
}


# https://github.com/broadinstitute/gatk/blob/master/scripts/unsupported/combine_tracks_postprocessing_cnv/combine_tracks.wdl
# https://github.com/bcbio/bcbio-nextgen/blob/520fad7da7fa46635e6514b8bb6b12b988cc20ac/bcbio/structural/titancna.py#L162

# https://github.com/gavinha/TitanCNA/issues/38
# floris/chapman posts here
# Convert GATK hets from ModelSegments to Titan input:
# cat SAMPLE.hets.tsv | awk '/^[^@]/ { print $1,$2,$5,$3,$6,$4 }' | tr ' ' '\t' > tmp.hets.tsv 


#' .titan_cn_file
#' 
#' Convert CNVkit or GATK4 normalized input into TitanCNA ready format.
#'
#' @param cnr_file 
#' 
#' @export
to_titan_cn_file <- function(cnr_file, 
                             titan_cnr_file,
                             file_type = c("gatkcnv", "cnvkit")){
  
  df_cnr = read_tsv(cnr_file, comment = "@", col_types = cols(.default = col_character()))
  
  file_type = match.arg(file_type)
  support_cols = list("cnvkit" = c("chromosome", "start", "end", "log2"),
                      "gatkcnv" = c("CONTIG", "START", "END", "LOG2_COPY_RATIO"))
  cnr_cols = support_cols[[file_type]]
  df_cnr = df_cnr[, support_cols[[file_type]]]
  # std col names
  # NAME BASED
  colnames(df_cnr) = c("chr", "start", "end", "log2_copy_ratio")
  # library(GenomicRanges)
  # library(TitanCNA)
  # df_cnr$chr <- setGenomeStyle(df_cnr$chr, genomeStyle = "NCBI")
  # gr_cnr <- data.frame(df_cnr, stringsAsFactors = F) %>% as("GRanges")		
  
  write_tsv(df_cnr, titan_cnr_file)
  
}


#' to_titan_hets_file
#'
#' @param hets_file 
#' @param titan_hets_file 
#'
#' @export
to_titan_hets_file <- function(hets_file, titan_hets_file){
  
  # cat SAMPLE.hets.tsv | awk '/^[^@]/ { print $1,$2,$5,$3,$6,$4 }' | tr ' ' '\t' > tmp.hets.tsv 
  
  # POSITON BASED
  df_hets = read_tsv(hets_file, comment = "@", col_types = cols(.default = col_character()))
  df_hets %<>% select(chr = CONTIG, pos = POSITION, 
                      ref_allele = REF_NUCLEOTIDE, 
                      ref_count = REF_COUNT,
                      alt_allele = ALT_NUCLEOTIDE, 
                      alt_count = ALT_COUNT)
  
  write_tsv(df_hets, titan_hets_file)
}




titancna_merge_opt_cluster <- function(df_trk){
  # df_trk =  df_trk
  i=1
  # df_trk$ttn_optimal_fl
  expect_columns(df_trk, "ttn_optimal_fl")
  tmp = lapply(1:nrow(df_trk), function(i){
    # being explicit is better for bind rows
    opt_clus = df_trk$ttn_optimal_fl[i] %>% 
      read_tsv(col_types = cols(
        Phi = col_double(),
        id = col_character(),
        barcode = col_character(),
        numClust = col_double(),
        cellPrev = col_character(),
        purity = col_double(),
        norm = col_double(),
        ploidy = col_double(),
        loglik = col_double(),
        sdbw = col_double(),
        path = col_character()
      )) %>% clean_names()
      opt_clus$individual = df_trk$individual[i]
      opt_clus
  }) %>% bind_rows()
  # tmp
  df_trk = left_join(df_trk, tmp, by = c("name" = "barcode", "individual"))
  # write_tsv(df_trk, "/rsrch3/home/iacs/sseth/flows/SS/sarco/mda/wex/titancna_v2/df_ttn_opt.tsv")
  df_trk
  }



# read titiancna -------
to_gr_seg.titancna <- function(segfl,
                               gencode_fl = "~/.rsrch3/home/iacs/sseth/ref/human/b37/annotations/gencode/v19/gencode.v19.annotation_gene_hgnc.bed",
                               lst_gen = NULL,
                               verbose = FALSE){
  
  seg = data.table::fread(segfl, data.table = F) %>% as_tibble() %>% 
    clean_names()
  seg
  seg$chromosome %>% table()
  
  # # http://seqanswers.com/forums/archive/index.php/t-32100.html
  # message("switch chr names")
  # seg <- seg %>% mutate(major_cn = tcn_em-lcn_em) %>% 
  #   dplyr::rename(chr = chrom, start = start, end = end, total_cn =  tcn_em, 
  #                 major_cn = major_cn, minor_cn = lcn_em, everything()) %>% 
  #   mutate(chr = gsub("23", "X", chr),
  #          chr = gsub("24", "Y", chr))
  # gr_seg <- makeGRangesFromDataFrame(seg, keep.extra.columns = T)
  
  # convert to GR
  lst = seg %>% to_gr_seg.seg(col_sample = "sample", 
                              col_chr = "chromosome", 
                              col_start = "start_position_bp", col_end = "end_position_bp", 
                              col_num_mark = "length_snp", col_seg_mean = "median_log_r")
  
  # read igvseg
  
  # read subclonal seg for igv
  segfl_sub = gsub(".segs.txt$", ".titan.txt", segfl)
  df_sub_seg = data.table::fread(segfl_sub, data.table = F)
  head(df_sub_seg)
  
  # convert to GR
  if(verbose)
    message("annotate")
  # use the supplied file
  if(is.null(lst_gen))
    lst_gen <- to_grseg.gencode(gencode_fl)

  source('~/Dropbox/projects/packs_dnaseq/R/to_mae.seg.R')
  df_ann = annotate.gr_seg(lst$gr_seg, lst_gen = lst_gen)
  dim(df_ann)
  # df_ann = df_ann %>% group_by(gene_name) %>% add_count()
  # df_ann %>% filter(n>1) %>% View()
  
  list(df_ann = df_ann, df_igv_seg = lst$df_igv_seg, gr_seg = lst$gr_seg)
  
}

# ** to_matrix -------
resolve_state <- function(x) {
  if (length(x) == 1) {
    return(x)
  }
  # mean(x, na.rm = T)
  # if(is.na(x)) return(x)
  x[which(abs(x) == max(abs(x)))][1]
}

# x=c("NEUT", "HETD")
resolve_corrected_call <- function(x) {
  if (length(x) == 1) {
    return(x)
  }
  indxs <- which(df_titan_ref$titan_corrected_call2 %in% x)
  df_titan_ref$titan_corrected_call2[max(indxs)]
}

# x=c("DLOH", "DHET")
resolve_call <- function(x) {
  if (length(x) == 1) {
    return(x)
  }
  indxs <- which(df_titan_ref$titan_call %in% x)
  df_titan_ref$titan_call[max(indxs)]
}


titan_cna_states <- function(){
  # based on discussions here
  # https://github.com/gavinha/TitanCNA/issues/8
  # read.delim(pipe("pbpaste")) %>% dput()

  df_titan_states <- structure(list(
    titan_state = -1:24,
    titan_genotype = c(
      "NULL",
      "NULL", "A", "AA", "AB", "AAA", "AAB", "AAAA", "AAAB", "AABB",
      "AAAAA", "AAAAB", "AAABB", "AAAAAA", "AAAAAB", "AAAABB", "AAABBB",
      "AAAAAAA", "AAAAAAB", "AAAAABB", "AAAABBB", "AAAAAAAA", "AAAAAAAB",
      "AAAAAABB", "AAAAABBB", "AAAABBBB"
    ),
    titan_total_copy_number = c(
      "NULL",
      "0", "1", "2", "2", "3", "3", "4", "4", "4", "5", "5", "5", "6",
      "6", "6", "6", "7", "7", "7", "7", "8", "8", "8", "8", "8"
    ),
    titan_call = c(
      "OUT", "HOMD", "DLOH", "NLOH", "HET", "ALOH",
      "GAIN", "ALOH", "ASCNA", "BCNA", "ALOH", "ASCNA", "UBCNA",
      "ALOH", "ASCNA", "UBCNA", "BCNA", "ALOH", "ASCNA", "UBCNA",
      "UBCNA", "ALOH", "ASCNA", "UBCNA", "UBCNA", "BCNA"
    ),
    titan_corrected_call = c(
      "OUT",
      "HOMD", "DLOH", "NLOH", "HET", "ALOH", "GAIN", "ALOH", "ASCNA",
      "BCNA", "ALOH", "ASCNA", "UBCNA", "ALOH", "ASCNA", "UBCNA",
      "BCNA", "ALOH", "ASCNA", "UBCNA", "UBCNA", "HLAMP", "HLAMP",
      "HLAMP", "HLAMP", "HLAMP"
    ),
    titan_corrected_call2 = c(
      "OUT",
      "HOMD", "HETD", "NLOH", "NEUT", "ALOH", "GAIN", "ALOH", "ASCNA",
      "BCNA", "ALOH", "ASCNA", "UBCNA", "ALOH", "ASCNA", "UBCNA",
      "BCNA", "ALOH", "AMP", "AMP", "AMP", "HLAMP", "HLAMP",
      "HLAMP", "HLAMP", "HLAMP"
    ),
    titan_call_description = c(
      "Outlier state",
      "Homozygous deletion", "Hemizygous deletion LOH", "Copy neutral LOH",
      "Diploid heterozygous", "Amplified LOH", "Gain/duplication of 1 allele",
      "Amplified LOH", "Allele-specific copy number amplification",
      "Balanced copy number amplification", "Amplified LOH", "Allele-specific copy number amplification",
      "Unbalanced copy number amplification", "Amplified LOH",
      "Allele-specific copy number amplification", "Unbalanced copy number amplification",
      "Balanced copy number amplification", "Amplified LOH", "Allele-specific copy number amplification",
      "Unbalanced copy number amplification", "Unbalanced copy number amplification",
      "Amplified LOH", "Allele-specific copy number amplification",
      "Unbalanced copy number amplification", "Unbalanced copy number amplification",
      "Balanced copy number amplification"
    ),
    titan_corrected_call_description = c(
      "Outlier state",
      "Homozygous deletion", "Hemizygous deletion LOH", "Copy neutral LOH",
      "Diploid heterozygous", "Amplified LOH", "Gain/duplication of 1 allele",
      "Amplified LOH", "Allele-specific copy number amplification",
      "Balanced copy number amplification", "Amplified LOH", "Allele-specific copy number amplification",
      "Unbalanced copy number amplification", "Amplified LOH",
      "Allele-specific copy number amplification", "Unbalanced copy number amplification",
      "Balanced copy number amplification", "Amplified LOH", "Allele-specific copy number amplification",
      "Unbalanced copy number amplification", "Unbalanced copy number amplification",
      "High-level amplification", "High-level amplification", "High-level amplification",
      "High-level amplification", "High-level amplification"
    )
  ), class = "data.frame", row.names = c(
    NA,
    -26L
  ))
  df_titan_states %>% dplyr::mutate(titan_total_copy_number = as.integer(titan_total_copy_number))
}
df_titan_ref <- titan_cna_states()


# plots -----
# compare purity, ploidy, and number of clusters

plot_purity_ploidy <- function(df_trk_ttn){
  # this is for 174 optcluster files
  p_load(ggplot2, cowplot, ggsci)
  p1 = df_trk_ttn %>% 
    ggplot(aes(purity, ploidy, color = as.factor(num_clust))) +
    geom_point() + theme_cowplot() +
    scale_color_npg(name = "# clusters")
  p2 = ggstatsplot::ggscatterstats(df_trk_ttn, 
    purity, ploidy, 
    type = "none", # type of test that needs to be run
    ggtheme = theme_cowplot(), # choosing a different theme

    formula = NA,
    smooth.line.args = list(size = 0, color = "black"),
    results.subtitle = FALSE,
    # points to label: 
    label.var = "individual", # variable for labeling data points
    label.expression = "ploidy > 5 | purity < 0.35", # expression that decides which points to label

    xfill = "pink", # color fill for x-axis marginal distribution
    yfill = "#009E73", # color fill for y-axis marginal distribution
    centrality.parameter = "median", # central tendency lines to be displayed
    centrality.label.args = list(size = 3, yintercept = 5.5),
    marginal.type = "density", # type of marginal distribution to be displayed

    messages = FALSE # turn off messages and notes
  )
  p = plot_grid(p1, p2, nrow=1)
  save_plot(glue("{odir}/p_purity_vs_ploidy.pdf"), p, base_width = 12, base_height = 6)

  # we have a few which have more than 4 ploidy, we did not even run
  # till that!
  filter(df_trk_ttn, ploidy > 4) %>% 
    select(phi, id, num_clust, purity, ploidy, loglik, sdbw, cell_prev) %>% 
    kable()


}