# code: https://github.com/mskcc/facets
# 
# ** install -------
# You can install the current version (along with the vignette) using the command
# devtools::install_github("mskcc/facets", build_vignettes = TRUE)
# pctGCdata is a required package. So install that also (needs to be done only once)
# dont install this: devtools::install_github("mskcc/pctGCdata")
# devtools::install_github("veseshan/pctGCdata")
#     this is v0.3
# 
# usage (brief):
# The new version estimates the log-ratio level corresponding to the diploid state.
# It is embedded into the procSample call. In terms of using the package you can now do:
# rcmat <- readSnpMatrix(filename, ...)
# xx <- preProcSample(rcmat, ...)
# # specify cval you like
# oo <- procSample(xx, cval = 300)
# 
# snp-pileup -g -q15 -Q20 -P100 -r25,0 vcffile outputfile normalbam tumorbam
# snp-pileupto  generate  the  read  countmatrix are given in the included README.txt
# pileup out:
# File#R, File#A, File#E and File#D giving the counts of number of reads with theref allele, alt allele,
# errors (neither ref nor alt) and deletions in that position.
# 
# 
# 
# notes/comments on the method:
# ndpeth: 35, and upper: 1000 (remove); ONLY in matched normal [not tumor]
# 
# they use a 150-250bp window, to space out the SNPs (preventing hypersegmentation) and
# local patterns of serial dependencies
# 
# focus on SNPs with allelic imbalance: VAF>.25 and <0.75
# usually yields: 250K snps from TCGA WEX
# 
# then calc, logR and logOR (not sure what OR is!)
# normalize for library depth + GC-bias (loess regression, over GC content)
# 
# IMP: preProcSample: samples SNPs to be USED; not sure on the impact!!
# The SNPs in a genome are not evenly spaced. Some regions have multiple SNPs in a small neighborhood.
# Thus using all loci will induce serial correlation in the data.
# To avoid it we sample loci such that only a single locus is used in an interval of length ‘snp.nbhd’.
# So in order to get reproducible results use ‘set.seed’ to fix the random number generator seed.
# 


# install in R 3.5.2
# module load conda_/3.6
# conda install htslib
# cd /rsrch2/iacs/iacs_dep/sseth/R/x86_64-pc-linux-gnu-library/3.5/facets/extcode
# g++ -std=c++11 snp-pileup.cpp -lhts -o snp-pileup


# installation:
# #module load conda_/2.7
# #conda install htslib # does not work!!
# #cd /rsrch2/iacs/iacs_dep/sseth/R/x86_64-pc-linux-gnu-library/3.4/facets/extcode
# #g++ -std=c++11 snp-pileup.cpp -lhts -o snp-pileup

# installing htslib: CRAM support needs newer libs [bgzip is built in]
# ./configure --prefix=/rsrch2/iacs/iacs_dep/sseth/apps/htslib/1.7 --disable-lzma
# make;make install


# using local htslib! only way this works!
# hts=/rsrch2/iacs/iacs_dep/sseth/apps/htslib/1.7
# g++ -std=c++11 -I$hts/include snp-pileup.cpp -L$hts/lib -lhts -Wl,-rpath=$hts/lib -o snp-pileup 

# snp-pileup=/rsrch2/iacs/iacs_dep/sseth/R/x86_64-pc-linux-gnu-library/3.4/facets/extcode

# JEFF
# I usually vary the critical value, nbhd (window), and nhet (het snps).
# nhet doesn’t make much of a difference, the window size cleans up a bit.  
# The critical value probably has the strongest influence.  
# Usually, we end up with values around critical value 300, nbhd 250, nhet 30.

# >    seems to be a good default:
# >    critical value 300
# >    ndepth 35
# >    nbhd 250
# >    nhet 30


# setup -------
# library(pacman)
# p_load(funr)
# p_load(facets, readr, dplyr, glue)



if(FALSE){
  
  rcfl = "~/.rsrch2/iacs/iacs_dep/sseth/tmp/test_clone_conv/185_006/facets/MS-PDX-300x-1728D-CB223_sahil-170713-MS-BREAST-PDX-SPANCER-1-NoIndex.bwa_recalibed.bampileup.tmp.gz"
  #rcmat2 = data.table::fread(rcfl, data.table = F, sep = ",")
  rcmat = readSnpMatrix(rcfl)
  set.seed(1234)
  
  #     tmp <- segsnps(dmat, cval, hetscale, deltaCN)
  # https://github.com/mskcc/facets/blob/bc1aa9a076192a2e4d6346a9fb74d2a3bfe6d619/R/facets-wrapper.R
  xx <- preProcSample(rcmat, 
                      snp.nbhd = 250, # window size, default: 250
                      cval = 300, # critical value for segmentation; default: 25
                      hetscale = T, # logical variable to indicate if logOR should get more weight
                      # in the test statistics for segmentation and clustering.
                      # Usually only 10% of snps are hets and hetscale gives the
                      # logOR contribution to T-square as 0.25/proportion of hets.
                      gbuild = "hg19"
  )
  
  # check if critical value makes a diff  
  xx2 <- preProcSample(rcmat, 
                       snp.nbhd = 250, # window size, default: 250
                       cval = 25, # critical value for segmentation; default: 25
                       hetscale = T, # logical variable to indicate if logOR should get more weight
                       # in the test statistics for segmentation and clustering.
                       # Usually only 10% of snps are hets and hetscale gives the
                       # logOR contribution to T-square as 0.25/proportion of hets.
                       gbuild = "hg19"
  )
  
  # values which depend on cval:
  #  "seg.tree"    "jointseg"  "hscl"        "chromlevels"
  
  #  We call the logR for the2-copy statelogR0which is output below.
  
  # preProcSample does a fine segmentation which procSample coarsens. 
  # This way focal changes can be identified and added back to broad changes and avoid hyper-segmentation due to wave artifacts.
  
  
  # specify cval you like
  oo <- procSample(xx, cval = 300, 
                   min.nhet = 30 # minimum number of heterozygote snps in a segment used for
                   # bivariate t-statistic during clustering of segments
  )
  
  # prevent hypersegmentation
  # Can you run preProcSample with a higher cval (say 50) and see if all the local changes you see are suppressed
  # I suggest you use 50 or 75 for preProcSample and bump cval to 300 for procSample.
  # preProcSample does a fine segmentation which procSample coarsens. 
  # This way focal changes can be identified and added back to broad changes and avoid hyper-segmentation due to wave artifacts.
  oo2 <- procSample(xx2, cval = 300, 
                    min.nhet = 30 # minimum number of heterozygote snps in a segment used for
                    # bivariate t-statistic during clustering of segments
  )
  
  
  oo$dipLogR
  # [1] -0.4337544
  
  # Call allele-speci c copy number and associated cellular fraction, estimate tumor purity and ploidy
  fit  = emcncf(oo)
  fit1.2  = emcncf(oo2)
  
  # Once  the  logR  value  for  the  diploid  state  is  obtained  we  calculate  the  observed  copynumber for each cluster as exp(logRclogR0)
  # where logRcis the logR summary for thecluster and logR0is the diploid state level.
  # Once the observed total number is obtainedwe obtain the allele speci c copy numbers m and p and the cellular fractionusing thelogOR data.
  # The cellular fraction is associated with the aberrant genotype.
  # For clonal copynumber alterations,equals tumor purity.  For subclonal events,will be lower than theoverall sample purity.

  pdf("genomewide.pdf")
  plotSample(x = oo, emfit = fit)
  logRlogORspider(oo$out, oo$dipLogR)
  dev.off()
  
  
  fit2  = emcncf2(oo)
  pdf("genomewide_fit2.pdf")
  plotSample(x = oo, emfit = fit2)
  dev.off()
  
  
  fit$purity
  #[1] 0.661544
  fit$ploidy
  #[1] 3.06038
  
  write_tsv(fit$cncf, "fit_cncf.tsv")
  write_tsv(fit2$cncf, "fit2_cncf.tsv")
  
  write_rds(xx, "xx.rds")
  write_rds(oo, "oo.rds")
  write_rds(fit, "fit.rds")
  write_rds(fit2, "fit2.rds")
  
}

# facets ART WEX B1 --------
if(FALSE){
  # funr devtools::install pkg=~/Dropbox/public/flow-r/my.ultraseq/my.ultraseq/
  
  library(pacman)
  p_load(flowr, tidyverse, glue, janitor)
  source('~/Dropbox/public/flow-r/my.ultraseq/my.ultraseq/R/facets.R')
  
  df_pairs = params::read_sheet("~/projects/samplesheets/tnbc/seq_gm_full_details.xlsx", sheet = "pt_pairs")
  df_pairs = df_pairs[complete.cases(df_pairs),]
  
  bampath = "/rsrch2/iacs/ngs_runs/sahil_170713_MS_BREAST_PDX_SPANCER/bams"
  bampath = "/rsrch2/iacs/ngs_runs/sahil_180301_pre16_art09/bams"
  odir = "~/.rsrch2/iacs/iacs_dep/sseth/tmp/test_clone_conv/facets"
  
  
  
  i = 1
  tmp = lapply(1:nrow(df_pairs), function(i){
    facets_flo(bam_t = file.path(bampath, df_pairs$bam1[i]),
               bam_n = file.path(bampath, df_pairs$bam2[i]),
               samplename = df_pairs$participantid[i],
               oprefix = df_pairs$samplelbl1[i])
  })
  flowmat = lapply(tmp, "[[", "flowmat") %>% bind_rows()
  write_tsv(flowmat, glue("{odir}/flowmat_b6.tsv"))
  setwd(odir)
  to_flow(x='flowmat_b6.tsv', def='flowdef.tsv', flowname='facets', flow_run_path=odir, execute=TRUE)
  #write_tsv(to_flowdef(flowmat), glue("{odir}/flowdef0.tsv"))
  # flowr to_flow x=flowmat.tsv def=flowdef.tsv flowname=facets flow_run_path=. execute=TRUE
  
  # add to flow:
  facetfls = lapply(tmp, "[[", "oprefix") %>% unlist()
  facetfits = glue("{facetfls}_fit.rds")
  
  # get pyclone flow:
  df = filter(df_pairs, participantid == "x185_002")
  
  purities = sapply(facetfits, facets_purity)
  
  # get facet files
  # get mut files
  
  #  
  
}


if(FALSE){
  library(facets)
  library(glue)
  setwd(dirname(rcfl))
  
  # b0
  rcfl = "~/.rsrch2/iacs/iacs_dep/sseth/tmp/test_clone_conv/185_006/facets/MS-PDX-300x-1728D-CB223_sahil-170713-MS-BREAST-PDX-SPANCER-1-NoIndex.bwa_recalibed.bampileup.tmp.gz"
  facets_r(rcfl,flowr::get_unique_id("facets_bam0"),pre_snp.nbhd=250,pre_cval=50,proc_cval=300,proc_min.nhet=30)
  facets_r(rcfl,flowr::get_unique_id("facets_bam0"),pre_snp.nbhd=250,pre_cval=100,proc_cval=300,proc_min.nhet=30)
  facets_r(rcfl,flowr::get_unique_id("facets_bam0"),pre_snp.nbhd=250,pre_cval=300,proc_cval=300,proc_min.nhet=30)
  
  # b1
  rcfl = "~/.rsrch2/iacs/iacs_dep/sseth/tmp/test_clone_conv/185_006/facets/MS-PDX-300x-1728D-CB224_sahil-170713-MS-BREAST-PDX-SPANCER-1-NoIndex.bwa_recalibed.bampileup.tmp.gz"
  facets_r(rcfl,flowr::get_unique_id("facets_bam1"),pre_snp.nbhd=250,pre_cval=50,proc_cval=300,proc_min.nhet=30)
  facets_r(rcfl,flowr::get_unique_id("facets_bam1"),pre_snp.nbhd=250,pre_cval=100,proc_cval=300,proc_min.nhet=30)
  facets_r(rcfl,flowr::get_unique_id("facets_bam1"),pre_snp.nbhd=250,pre_cval=25,proc_cval=300,proc_min.nhet=30)
  
  # bs
  rcfl = "~/.rsrch2/iacs/iacs_dep/sseth/tmp/test_clone_conv/185_006/facets/MS-PDX-300x-1728D-CB225_sahil-170713-MS-BREAST-PDX-SPANCER-1-NoIndex.bwa_recalibed.bampileup.tmp.gz"
  facets_r(rcfl,flowr::get_unique_id("facets_bams"),pre_snp.nbhd=250,pre_cval=50,proc_cval=300,proc_min.nhet=30)
  facets_r(rcfl,flowr::get_unique_id("facets_bams"),pre_snp.nbhd=250,pre_cval=100,proc_cval=300,proc_min.nhet=30)
  facets_r(rcfl,flowr::get_unique_id("facets_bams"),pre_snp.nbhd=250,pre_cval=300,proc_cval=300,proc_min.nhet=30)
  
}

#fitrds = "~/.rsrch2/iacs/iacs_dep/sseth/tmp/test_clone_conv/185_006/facets/x185_006_T0_fit.rds"
facets_purity <- function(fitrds){
  fit = read_rds(fitrds)
  round(fit$purity, 2)
}


facets_flo <- function(bam_t, bam_n, 
                       samplename,
                       oprefix,
                       
                       snp_pileup_exe = opts_flow$get("snp_pileup_exe"),
                       snp_pipeup_param = "-g -q15 -Q20 -P100 -r25,0",
                       vcffl = opts_flow$get("vcffl"),
                       rscript_exe = opts_flow$get("rscript_exe"),
                       funr = "funr",

                       # pre_snp.nbhd = 250,
                       # pre_cval = 50,
                       # proc_cval = 300,
                       # proc_min.nhet = 30,
                       # pre_hetscale = TRUE,
                       # pre_gbuild = "hg19",
                       # pre_unmatched = F,
                       # pre_het.thresh = 0.25,
                       # pre_ndepth = 35, # default
                       # 
                       # proc_cval = 300, 
                       # proc_min.nhet = 30,

                       # pre_snp.nbhd = opts_flow$get("pre_snp.nbhd"),
                       # pre_cval = opts_flow$get("pre_cval"),
                       # pre_hetscale = opts_flow$get("pre_hetscale"),
                       # pre_gbuild = opts_flow$get("pre_gbuild"),
                       # pre_unmatched = opts_flow$get("pre_unmatched"),
                       # pre_het.thresh = opts_flow$get("pre_het.thresh"),
                       # pre_ndepth = opts_flow$get("pre_ndepth"), # default
                       # 
                       # proc_cval = opts_flow$get("proc_cval"), 
                       # proc_min.nhet = opts_flow$get("proc_min.nhet"),
                       facets_params = "pre_snp.nbhd=250 pre_cval=50 pre_hetscale=TRUE pre_gbuild=hg19 pre_unmatched=FALSE pre_het.thresh=0.25 pre_ndepth=35 proc_min.nhet=30"
                       ){
  # >    seems to be a good default:
  # >    critical value 300
  # >    ndepth 35
  # >    nbhd 250
  # >    nhet 30
  
  # my params:
  # pre_snp.nbhd	250
  # pre_cval	50
  # proc_min.nhet	30
  # pre_hetscale	TRUE
  # pre_gbuild	"hg19"
  # pre_unmatched	FALSE
  # pre_het.thresh	0.25
  # pre_ndepth	35
  # 
  # proc_cval	300
  # proc_min.nhet	30
  
  pileupfl = glue("{oprefix}_pileup.tmp.gz")
  snp_pipeup_param = "-g -q15 -Q20 -P100 -r25,0"
  cmd0 = glue("{snp_pileup_exe} -g -q15 -Q20 -P100 -r25,0 {vcffl} {pileupfl} {bam_n} {bam_t}")
  cmd0 = as.character(cmd0)
  #cmd1 = glue("funr my.ultraseq::facets_r rcfl={pileupfl} pre_snp.nbhd=250 pre_cval=50 proc_cval=300 proc_min.nhet=30")
  
  # unique ID can be good, but can cause issues as well
  oprefix = flowr::get_unique_id(oprefix)
  # facets = "~/Dropbox/public/flow-r/my.ultraseq/my.ultraseq/R/facets.R"
  # cmd1 = glue("{rscript_exe} {facets} facets_r rcfl={pileupfl} oprefix={oprefix} pre_snp.nbhd=250 pre_cval=50 proc_cval=300 proc_min.nhet=30")
  cmd1 = glue("{funr} my.ultraseq::facets_r rcfl={pileupfl} oprefix={oprefix} {facets_params}")
  cmd1 = as.character(cmd1)
  
  lst = list(facets_snp = cmd0, facets_r = cmd1)
  flowmat = to_flowmat(lst, samplename = samplename)
  
  # flowr to_flow x=flowmat.tsv def=flowdef.tsv flowname=facets flow_run_path=. execute=TRUE
  outfile = glue("{oprefix}_cncf.tsv")
  
  list(flowmat = flowmat, oprefix = oprefix, outfile = outfile)
}


#' @name facets_r
#' 
#' @param rcfl snp matrix fl
#' @param oprefix something
#' @param pre_snp.nbhd something
#' @param pre_cval something
#' @param pre_hetscale something
#' @param pre_gbuild something
#' @param pre_unmatched indicator of whether the normal sample is unmatched. When this is TRUE hets are called using tumor reads only and logOR calculations are different. Use het.thresh = 0.1 or lower when TRUE.
#' @param het.thresh vaf threshold to call a SNP heterozygous
#' @param proc_cval something
#' @param proc_min.nhet something
#' 
#' @export
#' 
facets_r <- function(rcfl,
                     oprefix,
                     
                     pre_snp.nbhd = 250, 
                     pre_cval = 50,
                     pre_hetscale = T,
                     pre_gbuild = "hg19",
                     pre_unmatched = F,
                     pre_het.thresh = 0.25,
                     pre_ndepth = 35, # default
                     
                     proc_cval = 300, 
                     proc_min.nhet = 30){
  
  library(pctGCdata)
  
  # read it
  set.seed(1234)
  message("reading SNP matrix")
  rcmat = facets::readSnpMatrix(rcfl)
  
  # check if critical value makes a diff  
  message("preProcSample")
  xx <- facets::preProcSample(rcmat, 
                      snp.nbhd = pre_snp.nbhd, # window size, default: 250
                      cval = pre_cval, # critical value for segmentation; default: 25
                      hetscale = pre_hetscale, # logical variable to indicate if logOR should get more weight
                      # in the test statistics for segmentation and clustering.
                      # Usually only 10% of snps are hets and hetscale gives the
                      # logOR contribution to T-square as 0.25/proportion of hets.
                      # ndepth: minimum normal sample depth to keep
                      ndepth = pre_ndepth,
                      gbuild = pre_gbuild,
                      unmatched = pre_unmatched,
                      het.thresh = pre_het.thresh
  )
  
  # values which depend on cval:
  #  "seg.tree"    "jointseg"  "hscl"        "chromlevels"
  #  We call the logR for the2-copy statelogR0which is output below.
  # preProcSample does a fine segmentation which procSample coarsens. 
  # This way focal changes can be identified and added back to broad changes and avoid hyper-segmentation due to wave artifacts.
  message("procSample")
  oo <- facets::procSample(xx, 
                   cval = proc_cval, 
                   min.nhet = proc_min.nhet # minimum number of heterozygote snps in a segment used for
                   # bivariate t-statistic during clustering of segments
                   
                  # dipLogR: diploid level obtained from a fit, typically using a higher
                  # cval, can be used with lower cval to recover focal changes
  )
  
  # oo$dipLogR
  # [1] -0.4337544
  
  # Call allele-speci c copy number and associated cellular fraction, estimate tumor purity and ploidy
  fit  = emcncf(oo)
  
  # Once  the  logR  value  for  the  diploid  state  is  obtained  we  calculate  the  observed  copynumber for each cluster as exp(logRclogR0)
  # where logRcis the logR summary for thecluster and logR0is the diploid state level.
  # Once the observed total number is obtainedwe obtain the allele speci c copy numbers m and p and the cellular fractionusing thelogOR data.
  # The cellular fraction is associated with the aberrant genotype.
  # For clonal copynumber alterations,equals tumor purity.  For subclonal events,will be lower than theoverall sample purity.
  purity = round(fit$purity, 2)
  ploidy = round(fit$ploidy, 2)
  sname = glue("win:{pre_snp.nbhd},cval:{pre_cval},{proc_cval},nhet:{proc_min.nhet},purity:{purity},ploidy:{ploidy} ")
  sname
  
  pdf(glue("{oprefix}_fit_genomewide.pdf"))
  facets::plotSample(x = oo, emfit = fit, sname = sname)
  facets::plotSample(x = oo, emfit = fit, clustered = T, sname = sname)
  facets::logRlogORspider(oo$out, oo$dipLogR)
  dev.off()
  
  # fit2  = emcncf2(oo)
  # purity = round(fit2$purity, 2)
  # ploidy = round(fit2$ploidy, 2)
  # sname2 = glue("win:{pre_snp.nbhd},cval:{pre_cval},{proc_cval},nhet:{proc_min.nhet},purity:{purity},ploidy:{ploidy} ")
  # sname2
  # 
  # pdf(glue("{oprefix}_fit2_genomewide.pdf"))
  # plotSample(x = oo, emfit = fit2, sname = sname2)
  # plotSample(x = oo, emfit = fit2, clustered = T, sname = sname2)
  # logRlogORspider(oo$out, oo$dipLogR)
  # dev.off()
  
  
  readr::write_tsv(fit$cncf, glue("{oprefix}_fit_cncf.tsv"))
  #write_tsv(fit2$cncf, glue("{oprefix}_fit2_cncf.tsv"))
  
  readr::write_rds(xx, glue("{oprefix}_xx.rds"))
  readr::write_rds(oo, glue("{oprefix}_oo.rds"))
  readr::write_rds(fit, glue("{oprefix}_fit.rds"))
  #write_rds(fit2, glue("{oprefix}_fit2.rds"))
  
  invisible(list(xx = xx, oo = oo, fit = fit))
  # fit2 = fit2
  
}



facets_read_to_granges <- function(...){
  facets_read.seg(...)
}


# out: data frame of segment summaries pre and post clustering of
#       segments. The columns are: 
# ‘chrom’ the chromosome to which the segment belongs; ‘seg’ the segment number; 
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
# https://user-images.githubusercontent.com/1105264/53222479-4ff3cd80-36b1-11e9-8137-6839027c4364.png
# cncf: dataframe consisting of the columns of segmentation output as
#       well as cellular fraction (cf), total (tcn) and lesser (lcn)
#       copy number of each segment and their em counterpart (with
#       .em suffix)
read_facets.trk.seg <- function(df_trk, 
                                col_fl, 
                                col_samplename){
  
  df_trk = data.frame(df_trk, stringsAsFactors = F)
  df_seg = lapply(1:nrow(df_trk), function(i){
    
    samplename = df_trk[i, col_samplename]
    seg_fl = df_trk[i, col_fl]
    message(samplename, " ", appendLF = F)
    
    df_seg = read_tsv(seg_fl, 
                      col_types = cols(chrom = col_double(),
                                       seg = col_double(),
                                       num.mark = col_double(),
                                       nhet = col_double(),
                                       cnlr.median = col_double(),
                                       mafR = col_double(),
                                       segclust = col_double(),
                                       cnlr.median.clust = col_double(),
                                       mafR.clust = col_double(),
                                       start = col_double(),
                                       end = col_double(),
                                       cf.em = col_double(),
                                       tcn.em = col_double(),
                                       lcn.em = col_double())) %>% clean_names()
    df_seg %<>% mutate(samplename = samplename) %>% 
      mutate(major_cn = tcn_em - lcn_em) %>% 
      # this file is perfect for viewing in IGV
      dplyr::select(samplename, chrom, start, end,
                    num_mark = num_mark, 
                    seg_num = seg, 
                    
                    everything(), 
                    
                    log2_rd_r = cnlr_median, # log-ratio      (read depth ratio)      log2
                    log_baf_r = maf_r,       # log-odds-ratio (ratio of allele freqs) log_e  
                    
                    cluster = segclust,
                    minor_cn = lcn_em,
                    major_cn,
                    total_cn = tcn_em)
    df_seg
  }) %>% bind_rows()
  df_seg
  
}

# A SEG file (segmented data; .seg or .cbs) is a tab-delimited text file that lists loci and associated numeric values. 
# The first row contains column headings and each subsequent row contains a locus and an associated numeric value. 
# IGV ignores the column headings. 
# It reads the first four columns as track name, chromosome, start location, and end location. 
# It reads the last column as the numeric value for that locus (if the value is non-numeric, IGV ignores the row). 
# IGV ignores all other columns.
# 
# The segmented data file format is the output of the Circular Binary Segmentation algorithm (Olshen et al., 2004).


read_facets.trk.fits <- function(df_trk){
  
  df_facet_summ <- lapply(1:nrow(df_trk), function(i){
    message(i, appendLF = F)
    seg_fl = df_trk$facets_fl_full[i]
    fit_fl = gsub("_cncf.tsv", ".rds", seg_fl)
    
    fit = read_rds(fit_fl)
    data.frame(
      
      t_path_sample_id = df_trk$t_path_sample_id[i],
      n_path_sample_id = df_trk$n_path_sample_id[i],
      t_sample_type = df_trk$t_sample_type[i],
      path_patient_id = df_trk$path_patient_id[i],
      
      purity = fit$purity, 
      ploidy = fit$ploidy, 
      # if null, or absent, return NA
      loglik = ifelse(is.null(fit$loglik), NA, fit$loglik),
      dipLogR = fit$dipLogR)
  }) %>% bind_rows()
  df_facet_summ
}


#out = funr(commandArgs(trailingOnly = TRUE))


# END
