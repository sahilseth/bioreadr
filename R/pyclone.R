#' code:
#' 
#' instalation:
#' https://bitbucket.org/aroth85/pyclone/wiki/Installation
#' module load conda_/2.7
#' conda create --name pyclone #python=2
#' source activate pyclone
#' # now install it
#' conda install pyclone -c aroth85
#' 
#' usage:
#' 
#' 
#' issues:
#' 
#' - too many clusters: 
#' There is no parameter specify the number of clusters. If you are getting a large number of 
#' singleton clusters that is likely indicative of a problem with the analysis. 
#' You appear to be using an extremely large number of SNVs, PyClone is typically used with 100's not 1000's of SNVs.
#' That is not necessarily an issue, but you will likely need to run the MCMC analysis for much longer for the chain to converge. 
#' I would guess 100,000-1,000,000 iterations would be required based on my experience. 
#' You should also be running multiple restarts to check the results match.
#' 
#' pyclone is NOT for WGS, but for targetted:
#' For the purposes of PyClone analysis, less than 1000x is considered low coverage. 
#' For the 50x genomes the binomial may work better, and it will be slightly faster. 
#' I have to warn you that using PyClone on WGS data has not really been tested. 
#' We typically focus on targeted resequencing data. 
#' However, if you have multiple tumours per patient you may get reasonable results. 
#' Single tumour sample analysis will probably not work well though.
#' 
#' # Other ideas!
#' 1) If you have single sample data, expands may work. It is designed to scale to exome/WGS data.
#' 2) If you have multiple samples, but the copy number is relatively stable i.e. a 
#' large number of regions with no copy number change, then sciclone may work.
#' 3) If neither of the above hold you could run something like a Gaussian or Binomial 
#' mixture model on the data to cluster the SNVs by variant allele frequency (VAF). 
#' Then you can randomly sample SNVs from each cluster as representatives of the cluster and use these for PyClone analysis.
#' The reason 3) may work is as follows. SNVs with similar VAFs likely have the same mutational genotype 
#' and are present in the same clones. Thus SNVs in clusters from the GMM or BMM would effectively be treated the same by PyClone.
#' Hence, we can just look at a subset of them to get similar performance from PyClone.
#' 
#' # solution: RESTART!
#' Hi Andy, I have seen you mention "restart" in several posts. what do you mean by that?
#' running the same exact command but specify different seed?
#'
#'


# load ---------
library(pacman)
p_load(tidyverse, janitor)
p_load(GenomicRanges)



# pyclone TSV inputs ----------


#' get TSV files
#' focussing on somatic mutations
#' The required fields in this file are:
#' mutation_id - A unique ID to identify the mutation. Good names are thing such a the genomic co-ordinates of the mutation i.e. chr22:12345. Gene names are not good IDs because one gene may have multiple mutations, in which case the ID is not unique and PyClone will fail to run or worse give unexpected results. If you want to include the gene name I suggest adding the genomic coordinates i.e. TP53_chr17:753342.
#' ref_counts - The number of reads covering the mutation which contain the reference (genome) allele.
#' var_counts - The number of reads covering the mutation which contain the variant allele.
#' normal_cn - The copy number of the cells in the normal population. For autosomal chromosomes this will be 2 and for sex chromosomes it could be either 1 or 2. For species besides human other values are possible.
#' minor_cn - The minor copy number of the cancer cells. Usually this value will be predicted from WGSS or array data.
#' major_cn - The major copy number of the cancer cells. Usually this value will be predicted from WGSS or array data.
#' If you do not major and minor copy number information you should set the minor copy number to 0, and the major copy number to the predicted total copy number. If you do this make sure to use the total_copy_number for the --prior flag of the build_mutations_file, setup_analysis and run_analysis_pipeline commands. DO NOT use the parental copy number or major_copy_number information method as it assumes you have knowledge of the minor and major copy number.
#' Any additional columns in the tsv file will be ignored so feel free to add additional annotation fields.

# adapted from:
# https://gitlab.com/tangming2005/snakemake_DNAseq_pipeline/blob/lancet/scripts/generate_pyclone_input.R


if(FALSE){
  library(pacman)
  p_load(tidyverse, GenomicRanges, janitor)
  
  segfl = "~/.rsrch2/iacs/iacs_dep/sseth/tmp/test_clone_conv/185_006/185_006_v1/185_006_v1_031420181521069040_segments.txt"  
  # better since, possibility of aggressive filtering
  mutfl = "~/.rsrch2/iacs/ngs_runs/sahil_170713_MS_BREAST_PDX_SPANCER/mutect/MS-PDX-300x-1728D-CB223_sahil-170713-MS-BREAST-PDX-SPANCER-1--100x-1728D-CB226_sahil-170713-MS-BREAST-PDX-SPANCER-1.030220181519980509.annotated.tsv"
  #mutfl = "~/.rsrch2/iacs/ngs_runs/sahil_170713_MS_BREAST_PDX_SPANCER/mutect/MS-PDX-300x-1728D-CB223_sahil-170713-MS-BREAST-PDX-SPANCER-1--100x-1728D-CB226_sahil-170713-MS-BREAST-PDX-SPANCER-1_merged.muTect_call_stats.txt"
  
  pyc_inp = pyclone_prep_input.sequenza(segfl, mutfl)
  head(pyc_inp)
  # temp step, need to figure out why chrX is NA!
  pyc_inp = pyc_inp[complete.cases(pyc_inp), ]
  write_tsv(pyc_inp, "~/.rsrch2/iacs/iacs_dep/sseth/tmp/test_clone_conv/185_006/185_006_v1/pyclone_inp.tsv")
  
}

# https://gitlab.com/tangming2005/snakemake_DNAseq_pipeline/blob/lancet/scripts/sequenza.R


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

if(FALSE){
  segfl = "~/.rsrch2/iacs/iacs_dep/sseth/tmp/test_clone_conv/185_006/facets/facets_bam0-20180317-22-09-54-KyqNQY1j_fit_cncf.tsv"
}

# **** facets ---------

# *** facets column description -----------
#' data frame of segment summaries pre and post clustering of segments. 
#' The columns are: 
#' ‘chrom’ the chromosome to which the segment belongs; 
#' ‘seg’ the segment number; 
#' ‘num.mark’ the number of SNPs in the segment; 
#' ‘nhet’ the number of SNPs that are deemed heterozygous; 
#' ‘cnlr.median’ the median log-ratio of the segment; 
#' ‘mafR’ the log-odds-ratio summary for the segment; 
#' ‘segclust’ the segment cluster to which segment belongs; 
#' ‘cnlr.median.clust’ the median log-ratio of the segment cluster; 
#' ‘mafR.clust’ the log-odds-ratio summary for the segment cluster;
#' ‘cf’ the cellular fraction of the segment; 
#' ‘tcn’ the total copy number of the segment; 
#' ‘lcn’ the minor copy number of the segment.
#' 
#' dataframe consisting of the columns of segmentation output as well as 
#' cellular fraction (cf), total (tcn) and lesser (lcn) copy number of each segment 
#' and their em counterpart (with .em suffix)
 
# reading primarily for pyclone
# https://github.com/mskcc/facets/issues/57
# Filtering may be a better idea. These are typically focal changes. So If tcn is large you can see if including them with lcn=1 will give you sensible answers.
facets_read <- function(segfl){
  
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
  
  seg <- seg %>% mutate(major_cn = tcn_em-lcn_em) %>% 
    dplyr::rename(chr = chrom, start = start, end = end, major_cn = major_cn, minor_cn = lcn_em, total_cn =  tcn_em)
  gr_seg <- makeGRangesFromDataFrame(seg, keep.extra.columns = T)
  
  list(gr_seg = gr_seg, df_seg = seg)
                 
  
}

# ** mut readers ---------
# **** mutect ---------

# mutect_read.tsv
# mutect_read.vcf
# mutect_read.readcount
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

mutect_read.bamreadcount <- function(mutfl){
  mut <- read.table(mutfl, header = T, stringsAsFactors = F, sep = "\t") %>% clean_names()
  mut <- mut %>% dplyr::select(chr = chr, start = start, end = end, 
                  ref_counts = ref_count, var_counts = alt_count, 
                  tumor_name = samplename, 
                  ref_allele = ref_allele, alt_allele = alt_allele)
  gr_mut <- makeGRangesFromDataFrame(mut, keep.extra.columns = T)
  
  list(gr_mut = gr_mut, df_mut = mut)
}

# **** strelka ---------

# setup pyclone ---------
# **** prep TSVs ---------

#' prepare input for pyclone
#'
#' @param seg sequenza segment file
#' @param mut mutect mutations
#'
#' @details 
#' need to expand to support facets input
#' @export
#'
pyclone_prep_input <- function(segfl, mutfl, 
                               lst_bams = list(),
                               segfl_type = c("sequenza", "facets"), 
                               mutfl_type = c("mutect_ann", "mutect_bamreadcount"), 
                               outfile){
  
  if(segfl_type == "sequenza")
    gr_seg = seqz_read(segfl)$gr_seg
  if(segfl_type == "facets")
    gr_seg = facets_read(segfl)$gr_seg
  
  #mut<- read_tsv(mutfl, header =T, stringsAsFactors = F, sep = "\t")
  if(mutfl_type == "mutect_ann")
    gr_mut <- mutect_read.annovar(mutfl)$gr_mut
  if(mutfl_type == "mutect_bamreadcount")
    gr_mut <- mutect_read.bamreadcount(mutfl)$gr_mut
  
  # convert X Y to 23; TODO
  mut_seg_hits <- findOverlaps(gr_mut, gr_seg)
  
  df_merge = bind_cols(as.data.frame(gr_mut[queryHits(mut_seg_hits)]),
                       as.data.frame(gr_seg[subjectHits(mut_seg_hits)]) %>% 
                         dplyr::select(-seqnames, -start, -end, -width, -strand))
  df_merge = mutate(df_merge, 
                    mutation_id = paste("chr", seqnames, start, sep = ":"), # not including patient ID
                    #ref_counts = t_ref_AD, var_counts = t_alt_AD,
                    normal_cn = 2, variant_case = tumor_name, 
                    #variant_freq = t_alt_AD/(t_ref_AD + t_alt_AD),
                    variant_freq = var_counts/(ref_counts + var_counts),
                    genotype = paste(ref_allele, alt_allele, sep = ">")) %>%
    dplyr::select(mutation_id, ref_counts, var_counts, normal_cn, minor_cn, major_cn, 
                  variant_case, variant_freq, genotype)
  head(df_merge)
  
  # write out?'
  # pyc_inp = 
  df_merge[complete.cases(df_merge), ] %>% 
    filter(major_cn > 0) %>% 
    write_tsv(outfile)
  
  
  return(df_merge)
}


# https://gitlab.com/tangming2005/snakemake_DNAseq_pipeline/blob/lancet/Snakefile
### sometimes if one has a lot of variants (>1000), plotting function of PyClone will complain image size too big.
### if I run `PyClone run_analysis_pipeline`
### I have to run each step separately, and leave the plotting by R myself.
## https://groups.google.com/forum/#!topic/pyclone-user-group/UM9keeluDIQ
# this set_up step takes very short time, make it a local rule
## also make the backend of matplotlib to Agg if you are using pyclone on a cluster without X forwarding.

# PyClone build_table --config_file 16pyclone_output/Pa01_pyclone_analysis/config.yaml --out_file 16pyclone_output/Pa01_pyclone_analysis/tables/loci.tsv --table_type loci;
# PyClone build_table --config_file 16pyclone_output/Pa01_pyclone_analysis/config.yaml --out_file 16pyclone_output/Pa01_pyclone_analysis/tables/cluster.tsv --table_type cluster


if(FALSE){
  
  pyclone_exe = "module load conda_/2.7;source activate pyclone;PyClone"
  
  
}


# run pyclone ---------

# module load conda_/2.7;source activate pyclone
# PyClone 

# https://bitbucket.org/aroth85/pyclone/wiki/Usage
#' funcs, ending with fl, denote output is a flowmat type object
#' 
pyclone_flo <- function(pyclone_input, 
                        pyclone_exe = "module load conda_/2.7;source activate pyclone;PyClone",
                        samplename,
                        out_prefix){
  
  # this can be done, for pyclone_exe
  #cmd0 = "module load pyclone"
  cmd_build_mut = glue("{pyclone_exe} build_mutations_file --in_files {pyclone_input}")
  cmd1 = glue("{pyclone_exe} setup_analysis --in_file {pyclone_input} --tumour_contents {tumour_contents} --samples {params.samples} --density pyclone_binomial --working_dir {params.working_dir} --prior major_copy_number")
  cmd2 = glue("{pyclone_exe} run_analysis --config_file {input} --seed 1000")
  
}



#' run pyclone cmd
#' 
#'
#' @param df 
#' @param odir 
#' @param pyclone_params 
#' 
#' @details --num_iters 10000, usually, increase to 50000
#'
pyclone.run_analysis_pipeline <- function(df, 
                      odir, 
                      pyclone_params = "--num_iters 10000 --density pyclone_beta_binomial --prior parental_copy_number --burnin 2500"){
  
  pyclone_exe = "module load conda_/2.7;source activate pyclone;PyClone"
  
  # run pyclone
  mutfls = paste(df$pyclone_inp_tsv, collapse = " ")
  tum_contents = paste(round(df$purity, 2), collapse = " ")
  cmd_pyc = glue("{pyclone_exe} run_analysis_pipeline --in_files {mutfls} --working_dir {odir} --tumour_contents {tum_contents} {pyclone_params}")
  #system(cmd_pyc)
  
  # parse loci
  #pyclone_parse_cluster(df, odir)
  # run clonal evol
  
  # get trees
  
  # get sig muts
  
  # run pathway analysis
  
  return(list(cmd = cmd_pyc))
  
}


# df = df_pyclone_inp
# odir = 
pyclone_parse_cluster <- function(df, odir = "."){
  
  # read clusters:
  clust_fl = glue("{odir}/tables/cluster.tsv");
  loci_fl = glue("{odir}/tables/loci.tsv");
  lociann_fl = glue("{odir}/tables/loci_ann.tsv");
  
  # conf and size of different clusters
  df_clust = read_tsv(clust_fl)
  df_loci = read_tsv(loci_fl)

  df_mut = read_tsv(df$mutfl[1]) %>% 
    mutate(key_pyclone = paste0("chr", ":", chr, ":", start, ""))
  df_loci_ann = left_join(df_loci, df_mut, by = c("mutation_id" = "key_pyclone")) %>% 
    left_join(df_clust, by = c("sample_id", "cluster_id"))
  df_loci_ann
  
  write_tsv(df_loci_ann, lociann_fl)
  
  invisible(list(df_loci_ann = df_loci_ann))
}



patid = "002"
pyclone_pipe_r <- function(trk, 
                           segfl_type = "facets", 
                           mutfl_type = "mutect_bamreadcount"){
  #patid = paste(df$trialid[1], "_", df$path_patient_id[1], sep = "")
  patid = trk$participantid[1]
  message(patid)
  
  i = 2
  df_pyclone_inp = lapply(1:nrow(trk), function(i){
    sampid = trk$samplelbl1[i]
    outfile = glue("{sampid}_pyclone_inp.tsv");outfile
    fitfl = gsub("_cncf.tsv", ".rds", trk$segfl[i])
    segfl = trk$segfl[i]
    mutfl = trk$mutfl[i]
    source('~/Dropbox/public/flow-r/my.ultraseq/my.ultraseq/R/pyclone.R')
    pyc_inp = pyclone_prep_input(segfl, mutfl, 
                                 segfl_type = segfl_type, 
                                 mutfl_type = mutfl_type, 
                                 outfile = outfile)
    fit = read_rds(fitfl)
    
    data.frame(pyclone_inp_tsv = outfile, 
               purity = fit$purity, 
               ploidy = fit$ploidy, stringsAsFactors = F,
               mutfl = mutfl, segfl = segfl, 
               pyclone_inp_fl = outfile)
  }) %>% bind_rows()
  df_pyclone_inp
  
  # run pyclone cmds
  source('~/Dropbox/public/flow-r/my.ultraseq/my.ultraseq/R/pyclone.R')
  #py_run_path = glue("{pyclonepath}/{run_ver}/{patid}");py_run_path
  #if(!dir.exists(py_run_path)) dir.create(py_run_path, recursive = T)
  lst = pyclone.run_analysis_pipeline(df_pyclone_inp, odir = ".")
  paste0(lst$cmd, " > pyclone.log 2>&1") %>% system()
  
  # annotate pyclone output
  source('~/Dropbox/public/flow-r/my.ultraseq/my.ultraseq/R/pyclone.R')
  out = pyclone_parse_cluster(df_pyclone_inp)
  
  # pyclone_parse_cloneevol
  
  
}











# END
