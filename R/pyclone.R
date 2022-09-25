# code:
# 
# instalation:
# https://bitbucket.org/aroth85/pyclone/wiki/Installation
# module load conda_/2.7
# conda create --name pyclone_py27 python=2.7
# source activate pyclone
# # now install it
# conda install pyclone -c aroth85
# 
# usage:
# 
# 
# issues:
# 
# - too many clusters:
# There is no parameter specify the number of clusters. If you are getting a large number of
# singleton clusters that is likely indicative of a problem with the analysis.
# You appear to be using an extremely large number of SNVs, PyClone is typically used with 100's not 1000's of SNVs.
# That is not necessarily an issue, but you will likely need to run the MCMC analysis for much longer for the chain to converge.
# I would guess 100,000-1,000,000 iterations would be required based on my experience.
# You should also be running multiple restarts to check the results match.
# 
# pyclone is NOT for WGS, but for targetted:
# For the purposes of PyClone analysis, less than 1000x is considered low coverage.
# For the 50x genomes the binomial may work better, and it will be slightly faster.
# I have to warn you that using PyClone on WGS data has not really been tested.
# We typically focus on targeted resequencing data.
# However, if you have multiple tumours per patient you may get reasonable results.
# Single tumour sample analysis will probably not work well though.
# 
# # Other ideas!
# 1) If you have single sample data, expands may work. It is designed to scale to exome/WGS data.
# 2) If you have multiple samples, but the copy number is relatively stable i.e. a
# large number of regions with no copy number change, then sciclone may work.
# 3) If neither of the above hold you could run something like a Gaussian or Binomial
# mixture model on the data to cluster the SNVs by variant allele frequency (VAF).
# Then you can randomly sample SNVs from each cluster as representatives of the cluster and use these for PyClone analysis.
# The reason 3) may work is as follows. SNVs with similar VAFs likely have the same mutational genotype
# and are present in the same clones. Thus SNVs in clusters from the GMM or BMM would effectively be treated the same by PyClone.
# Hence, we can just look at a subset of them to get similar performance from PyClone.
# 
# # solution: RESTART!
# Hi Andy, I have seen you mention "restart" in several posts. what do you mean by that?
# running the same exact command but specify different seed?
# 


# install pyclone-vi
# source activate pyclone-vi
# pip install git+https://github.com/Roth-Lab/pyclone-vi.git



# library(pacman)
# p_load(tidyverse, janitor)
# p_load(GenomicRanges)





# ** seg readers ---------
# **** sequenza ---------
seqz_read <- function(segfl){
  seg = read_tsv(segfl, 
                 col_types = c(chromosome = col_character(),
                               start.pos = col_integer(),
                               end.pos = col_integer(),
                               Bf = col_double(),
                               N.BAF = col_integer(),
                               sd.BAF = col_double(),
                               depth.ratio = col_double(),
                               N.ratio = col_integer(),
                               sd.ratio = col_double(),
                               CNt = col_integer(),
                               A = col_integer(),
                               B = col_integer(),
                               LPP = col_double())) %>% clean_names()
  seg <- seg %>% 
    dplyr::rename(chr = chromosome, start = start_pos, end = end_pos, major_cn = a, minor_cn = b, total_cn =  cnt)
  gr_seg <- makeGRangesFromDataFrame(seg, keep.extra.columns = T)
  
  list(gr_seg = gr_seg, df_seg = seg)
}


# **** facets ---------

# *** facets column description -----------
# data frame of segment summaries pre and post clustering of segments.
# The columns are:
# ‘chrom’ the chromosome to which the segment belongs;
# ‘seg’ the segment number;
# ‘num.mark’ the number of SNPs in the segment;
# ‘nhet’ the number of SNPs that are deemed heterozygous;
# ‘cnlr.median’ the median log-ratio of the segment;
# ‘mafR’ the log-odds-ratio summary for the segment;
# ‘segclust’ the segment cluster to which segment belongs;
# ‘cnlr.median.clust’ the median log-ratio of the segment cluster;
# ‘mafR.clust’ the log-odds-ratio summary for the segment cluster;
# ‘cf’ the cellular fraction of the segment;
# ‘tcn’ the total copy number of the segment;
# ‘lcn’ the minor copy number of the segment.
# 
# dataframe consisting of the columns of segmentation output as well as
# cellular fraction (cf), total (tcn) and lesser (lcn) copy number of each segment
# and their em counterpart (with .em suffix)
 
# reading primarily for pyclone
# https://github.com/mskcc/facets/issues/57
# Filtering may be a better idea. These are typically focal changes. So If tcn is large you can see if including them with lcn=1 will give you sensible answers.
# DF: chrom, 
#     start end,
#     major_cn minor_cn total_cn
#     clonality
# EXTRA
facets_read.segfl <- function(segfl){
  
  seg = read_tsv(segfl, 
                 col_types = c(chrom = col_integer(),
                               seg = col_integer(),
                               num.mark = col_integer(),
                               nhet = col_integer(),
                               cnlr.median = col_double(),
                               mafR = col_double(),
                               segclust = col_integer(),
                               cnlr.median.clust = col_double(),
                               mafR.clust = col_double(),
                               start = col_double(),
                               end = col_double(),
                               cf.em = col_double(),
                               tcn.em = col_integer(),
                               lcn.em = col_integer())) %>% clean_names()
  # http://seqanswers.com/forums/archive/index.php/t-32100.html
  # 
  message("switch chr names")
  seg <- seg %>% mutate(major_cn = tcn_em-lcn_em) %>% 
    dplyr::rename(chr = chrom, start = start, end = end, total_cn =  tcn_em, 
                  major_cn = major_cn, minor_cn = lcn_em, everything()) %>% 
    mutate(chr = gsub("23", "X", chr),
           chr = gsub("24", "Y", chr))
  gr_seg <- makeGRangesFromDataFrame(seg, keep.extra.columns = T)
  
  list(gr_seg = gr_seg, df_seg = seg)
}

# gm_read_facets.seg <- function(df_trk,
#                                col_fl = "exomcn_fl_full", 
#                                col_samplename = "t_path_sample_id"){
#     
#     lst_seg = lapply(1:nrow(df_trk), function(i) {
#       samplename = df_trk[i, col_samplename] %>% unlist()
#       seg_fl = df_trk[i, col_fl] %>% unlist()
#       # print(class(seg_fl))
#       message(samplename, appendLF = F)
#       lst_seg = facets_read.seg(seg_fl)
#       lst_seg$df_seg %<>% clean_names() %>% mutate(samplename = samplename)
#       lst_seg
#     }) 
#     lst_seg = transpose(lst_seg)
#     df_seg %>% bind_rows(lst_seg$df_seg)
#     df_seg
# }

facets_segfl_fitfl <- function(seg_fl){
  fit_fl = gsub("_cncf.tsv", ".rds", seg_fl)
  
}

# **** purecn -----
# segfl = "/rsrch3/scratch/iacs/sseth/flows/SS/tnbc/ms51_wex/ss_cl_het_gb/purecn_v2/185_011/IPCT-S4012-MOON0051-Cap2515-4-HTID331___matched.csv"
# segfl = "/rsrch3/scratch/iacs/sseth/flows/SS/tnbc/ms51_wex/ss_cl_het_gb/purecn_v2/185_011/IPCT-S4012-MOON0051-Cap2515-4-HTID331___matched.rds"
# segfl = "/rsrch3/scratch/iacs/sseth/flows/SS/tnbc/ms51_wex/ss_cl_het_merged/cnv/purecn_v2/seg_mrg.rds"
# Rscript ~/Dropbox/public/flowr/my.ultraseq/pipelines/dnaseq/purecn/PureCN_v2.R --out purecn2/IPCT-S4012-MOON0051-Cap2515-4-HTID331___matched --sampleid IPCT-S4012-MOON0051-Cap2515-4-HTID331 --tumor /rsrch3/scratch/iacs/sseth/flows/SS/tnbc/ms51_wex/purecn/coverage/IPCT-S4012-MOON0051-Cap2515-4-HTID331_ngs-pipe-3-TCCGGAGCAGAA.bwa_recalibed_coverage_loess.txt.gz --normaldb /rsrch3/home/iacs/sseth/ref/human/b37/common_normal/ms51_v1/purecn_pon/pon_db/normalDB_seqez_v3_hg19.rds --mappingbiasfile /rsrch3/home/iacs/sseth/ref/human/b37/common_normal/ms51_v1/purecn_pon/pon_db/mapping_bias_seqez_v3_hg19.rds --vcf mutect2/IPCT-S4012-MOON0051-Cap2515-4-HTID331___matched.vcf.gz --outvcf --cosmic.vcf.file /rsrch3/home/iacs/sseth/ref/human/b37/annotations/cosmic/v88/CosmicCodingMuts.vcf.gz --intervals /rsrch3/home/iacs/sseth/ref/az_ref_beds/ss_downloads/SeqCapEZ_Exome_v3.0_Design_Annotation_files/purecn/SeqCap_EZ_Exome_v3_hg19_capture_targets_nochr_purecnv_interval.txt --snpblacklist /rsrch3/home/iacs/sseth/ref/human/b37/annotations/purecn/hg19_simple_repeats.bed.gz --parallel --cores 24 --postoptimize --seed 123 --minaf 0.03 --genome hg19 --funsegmentation PSCBS --centromere_seq_style NCBI   --additionaltumors /rsrch3/scratch/iacs/sseth/flows/SS/tnbc/ms51_wex/purecn/coverage/IPCT-S4012-MOON0051-Cap2516-4-HTID343_ngs-pipe-3-GATATCTAATTT.bwa_recalibed_coverage_loess.txt.gz
# pureCN was run with information from multiple tumors, where available
purecn_read.segfl <- function(segfl){
  seg = read_tsv(segfl)
  seg = read_rds(segfl)
  seg$results %>% length()
  x1 = seg$results[[1]]$seg
  x2 = seg$results[[5]]$seg
  dim(x1);dim(x2)
  table(x1$C, x2$C)
  
  res <- PureCN::readCurationFile(segfl)
  
}

# ** mut readers ---------
# **** mutect ---------

# mutect_read.tsv
# mutect_read.vcf
# mutect_read.annovar
mutect_read.annovar <- function(mutfl){
  mut <- read.table(mutfl, header = T, stringsAsFactors = F, sep = "\t") %>% clean_names()
  mut <- mut %>% 
    dplyr::select(chr = chr, start = start, end = end, 
                  ref_counts = t_ref_count, var_counts = t_alt_count, 
                  tumor_name = tumor_name, 
                  ref_allele = ref_allele, alt_allele = alt_allele)
  gr_mut <- makeGRangesFromDataFrame(mut, keep.extra.columns = T)
  
  list(gr_mut = gr_mut, df_mut = mut)
}

# mutect_read.readcount
mutect_read.bamreadcount <- function(mutfl){
  mut <- read.table(mutfl, header = T, stringsAsFactors = F, sep = "\t") %>% clean_names()
  mut <- mut %>% dplyr::select(chr = chr, start = start, end = end, 
                  ref_counts = ref_count, var_counts = alt_count, 
                  tumor_name = samplename, 
                  ref_allele = ref_allele, alt_allele = alt_allele)
  gr_mut <- makeGRangesFromDataFrame(mut, keep.extra.columns = T)
  
  list(gr_mut = gr_mut, df_mut = mut)
}

# **** strelka ----




# END
