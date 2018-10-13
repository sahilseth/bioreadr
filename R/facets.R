#' code: https://github.com/mskcc/facets
#' 
#' install:
#' You can install the current version (along with the vignette) using the command
#' devtools::install_github("mskcc/facets", build_vignettes = TRUE)
#' pctGCdata is a required package. So install that also (needs to be done only once)
#' devtools::install_github("mskcc/pctGCdata")
#' 
#' usage (brief):
#' The new version estimates the log-ratio level corresponding to the diploid state. 
#' It is embedded into the procSample call. In terms of using the package you can now do:
#' rcmat <- readSnpMatrix(filename, ...)
#' xx <- preProcSample(rcmat, ...)
#' # specify cval you like
#' oo <- procSample(xx, cval = 300)
#' 
#' snp-pileup -g -q15 -Q20 -P100 -r25,0 vcffile outputfile normalbam tumorbam
#' snp-pileupto  generate  the  read  countmatrix are given in the included README.txt
#' pileup out: 
#' File#R, File#A, File#E and File#D giving the counts of number of reads with theref allele, alt allele, 
#' errors (neither ref nor alt) and deletions in that position. 
#'
#'
#'
#' notes/comments on the method:
#' ndpeth: 35, and upper: 1000 (remove); ONLY in matched normal [not tumor]
#' 
#' they use a 150-250bp window, to space out the SNPs (preventing hypersegmentation) and
#' local patterns of serial dependencies
#' 
#' focus on SNPs with allelic imbalance: VAF>.25 and <0.75
#' usually yields: 250K snps from TCGA WEX
#' 
#' then calc, logR and logOR (not sure what OR is!)
#' normalize for library depth + GC-bias (loess regression, over GC content)
#' 
#' IMP: preProcSample: samples SNPs to be USED; not sure on the impact!!
#' The SNPs in a genome are not evenly spaced. Some regions have multiple SNPs in a small neighborhood.
#' Thus using all loci will induce serial correlation in the data. 
#' To avoid it we sample loci such that only a single locus is used in an interval of length ‘snp.nbhd’. 
#' So in order to get reproducible results use ‘set.seed’ to fix the random number generator seed.
#' 
#' 

# installation:
# module load conda_/2.7
# conda install htslib
# cd /rsrch2/iacs/iacs_dep/sseth/R/x86_64-pc-linux-gnu-library/3.4/facets/extcode
# g++ -std=c++11 snp-pileup.cpp -lhts -o snp-pileup

# using local htslib! only way this works!
# hts=/rsrch2/iacs/iacs_dep/sseth/apps/htslib/1.7
# g++ -std=c++11 -I$hts/include snp-pileup.cpp -L$hts/lib -lhts -Wl,-rpath=$hts/lib -o snp-pileup 



# installing htslib: CRAM support needs newer libs [bgzip is built in]
# ./configure --prefix=/rsrch2/iacs/iacs_dep/sseth/apps/htslib/1.7 --disable-lzma
# make;make install

# snp-pileup=/rsrch2/iacs/iacs_dep/sseth/R/x86_64-pc-linux-gnu-library/3.4/facets/extcode

# JEFF
# I usually vary the critical value, nbhd (window), and nhet (het snps).
# nhet doesn’t make much of a difference, the window size cleans up a bit.  
# The critical value probably has the strongest influence.  
# Usually, we end up with values around critical value 300, nbhd 250, nhet 30.


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
  
  #' values which depend on cval:
  #'  "seg.tree"    "jointseg"  "hscl"        "chromlevels"
  
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
  
  #' Once  the  logR  value  for  the  diploid  state  is  obtained  we  calculate  the  observed  copynumber for each cluster as exp(logRclogR0) 
  #' where logRcis the logR summary for thecluster and logR0is the diploid state level.  
  #' Once the observed total number is obtainedwe obtain the allele speci c copy numbers m and p and the cellular fractionusing thelogOR data. 
  #' The cellular fraction is associated with the aberrant genotype.  
  #' For clonal copynumber alterations,equals tumor purity.  For subclonal events,will be lower than theoverall sample purity.
  
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

#fitrds = "~/.rsrch2/iacs/iacs_dep/sseth/tmp/test_clone_conv/185_006/facets/x185_006_T0_fit.rds"
facets_purity <- function(fitrds){
  fit = read_rds(fitrds)
  round(fit$purity, 2)
}


facets_flo <- function(bam_t, bam_n, 
                       samplename,
                       oprefix,
                       snp_pileup_exe = "/rsrch2/iacs/iacs_dep/sseth/R/x86_64-pc-linux-gnu-library/3.4/facets/extcode/snp-pileup",
                       vcffl = '/rsrch2/iacs/iacs_dep/sseth/ref/human/b37/annotations/ncbi/00-common_all.vcf'){
  
  pileupfl = glue("{oprefix}_pileup.tmp.gz")
  cmd0 = glue("{snp_pileup_exe} -g -q15 -Q20 -P100 -r25,0 {vcffl} {pileupfl} {bam_n} {bam_t}")
  #cmd1 = glue("funr my.ultraseq::facets_r rcfl={pileupfl} pre_snp.nbhd=250 pre_cval=50 proc_cval=300 proc_min.nhet=30")
  oprefix = flowr::get_unique_id(oprefix)
  facets = "~/Dropbox/public/flow-r/my.ultraseq/my.ultraseq/R/facets.R"
  cmd1 = glue("Rscript {facets} facets_r rcfl={pileupfl} oprefix={oprefix} pre_snp.nbhd=250 pre_cval=50 proc_cval=300 proc_min.nhet=30")
  
  lst = list(facets_snp = cmd0, facets_r = cmd1)
  flowmat = to_flowmat(lst, samplename = samplename)
  
  # flowr to_flow x=flowmat.tsv def=flowdef.tsv flowname=facets flow_run_path=. execute=TRUE
  outfile = glue("{oprefix}_cncf.tsv")
  
  list(flowmat = flowmat, oprefix = oprefix, outfile = outfile)
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

#' facets_r
#'
#' @param rcfl 
#' @param oprefix 
#' @param pre_snp.nbhd 
#' @param pre_cval 
#' @param pre_hetscale 
#' @param pre_gbuild 
#' @param pre_unmatched indicator of whether the normal sample is unmatched.  
#' When this is TRUE hets are called using tumor reads only and logOR calculations are different. 
#' Use het.thresh = 0.1 or lower when TRUE.
#' @param het.thresh vaf threshold to call a SNP heterozygous
#' @param proc_cval 
#' @param proc_min.nhet 
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
                     
                     proc_cval = 300, 
                     proc_min.nhet = 30){
  
  # read it
  set.seed(1234)
  rcmat = readSnpMatrix(rcfl)
  
  # check if critical value makes a diff  
  xx <- preProcSample(rcmat, 
                      snp.nbhd = pre_snp.nbhd, # window size, default: 250
                      cval = pre_cval, # critical value for segmentation; default: 25
                      hetscale = pre_hetscale, # logical variable to indicate if logOR should get more weight
                      # in the test statistics for segmentation and clustering.
                      # Usually only 10% of snps are hets and hetscale gives the
                      # logOR contribution to T-square as 0.25/proportion of hets.
                      gbuild = pre_gbuild,
                      unmatched = pre_unmatched,
                      het.thresh = pre_het.thresh
  )
  
  #' values which depend on cval:
  #'  "seg.tree"    "jointseg"  "hscl"        "chromlevels"
  #  We call the logR for the2-copy statelogR0which is output below.
  # preProcSample does a fine segmentation which procSample coarsens. 
  # This way focal changes can be identified and added back to broad changes and avoid hyper-segmentation due to wave artifacts.
  oo <- procSample(xx, 
                   cval = proc_cval, 
                   min.nhet = proc_min.nhet # minimum number of heterozygote snps in a segment used for
                   # bivariate t-statistic during clustering of segments
  )
  oo$dipLogR
  # [1] -0.4337544
  
  # Call allele-speci c copy number and associated cellular fraction, estimate tumor purity and ploidy
  fit  = emcncf(oo)
  
  #' Once  the  logR  value  for  the  diploid  state  is  obtained  we  calculate  the  observed  copynumber for each cluster as exp(logRclogR0) 
  #' where logRcis the logR summary for thecluster and logR0is the diploid state level.  
  #' Once the observed total number is obtainedwe obtain the allele speci c copy numbers m and p and the cellular fractionusing thelogOR data. 
  #' The cellular fraction is associated with the aberrant genotype.  
  #' For clonal copynumber alterations,equals tumor purity.  For subclonal events,will be lower than theoverall sample purity.
  purity = round(fit$purity, 2)
  ploidy = round(fit$ploidy, 2)
  sname = glue("win:{pre_snp.nbhd},cval:{pre_cval},{proc_cval},nhet:{proc_min.nhet},purity:{purity},ploidy:{ploidy} ")
  sname
  
  pdf(glue("{oprefix}_fit_genomewide.pdf"))
  plotSample(x = oo, emfit = fit, sname = sname)
  plotSample(x = oo, emfit = fit, clustered = T, sname = sname)
  logRlogORspider(oo$out, oo$dipLogR)
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
  
  
  write_tsv(fit$cncf, glue("{oprefix}_fit_cncf.tsv"))
  #write_tsv(fit2$cncf, glue("{oprefix}_fit2_cncf.tsv"))
  
  write_rds(xx, glue("{oprefix}_xx.rds"))
  write_rds(oo, glue("{oprefix}_oo.rds"))
  write_rds(fit, glue("{oprefix}_fit.rds"))
  #write_rds(fit2, glue("{oprefix}_fit2.rds"))
  
  invisible(list(xx = xx, oo = oo, fit = fit))
  #, fit2 = fit2
  
}



#out = funr(commandArgs(trailingOnly = TRUE))


# END