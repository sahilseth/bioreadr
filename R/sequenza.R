
# help:
# http://www.cbs.dtu.dk/biotools/sequenza/

# mpileup:
# "/scratch/rists/hpcapps/x86_64/samtools/0.1.19/samtools mpileup -l /scratch/iacs/ngs/captureBeds/E1_sureSelectV4_hg19.bed -f ~/.rsrch2/iacs/iacs_dep/sseth/ref/human/b37/fastas/Homo_sapiens_assembly19.fasta -Q 20 ~/.rsrch1/iacs/ngs_runs/160502_SN1222_0345_BC8P7HACXX/bams/GDraetts-CSMS-G001-PATC53R15-T_C8P7HACXX-1-AGCCATGC.bwa_recalibed.bam | gzip > ~/.rsrch1/iacs/ngs_runs/160502_SN1222_0345_BC8P7HACXX/sequenza/PATC53R15---commN1_lung_blood_11_sample/GDraetts-CSMS-G001-PATC53R15-T_090520171504652725_pileup.gz"

# https://cran.r-project.org/web/packages/sequenza/vignettes/sequenza.pdf
# 1 samtools mpileup −f hg19.fasta −Q 20 normal.bam | gzip > normal.pileup.gz
# 2 samtools mpileup −f hg19.fasta −Q 20 tumor.bam | gzip > tumor.pileup.gz

# need to do ONCE
# sequenza−utils.py GC−windows −w 50 hg19.fa | gzip > hg19.gc50Base.txt.gz


# sequenza−utils.py bam2seqz −gc hg19.gc50Base.txt.gz
# −−fasta hg19.fasta \
#  −n normal.fifo \
# −t tumor.fifo | gzip > out.seqz.gz

# sequenza−utils.py seqz−binning −w 50 \
# −s out.seqz.gz | gzip > out small.seqz.gz

# 4.3
# gc.stats <- gc.sample.stats(data.file)

# test <- sequenza.extract(data.file)

# CP.example <- sequenza.fit(test)


# issues  ---------
# ** v high ploidy
# https://groups.google.com/forum/#!topic/sequenza-user-group/LpqHxZnLkks

# ** vs facets --------
# https://groups.google.com/forum/#!searchin/sequenza-user-group/segmentation%7Csort:date/sequenza-user-group/pDzhj440nUA/qf-_p1JVAAAJ
# facets gives: provide sublclonal heterogeneity information
#
# Regarding the pre-processing of sequencing data, 
# sequenza genotypes the samples by looking at the base content of the normal sample, 
# while FACETS can select SNP from dbSNP and the 1k genome projects.
#
# So the segments resulting from the two tools may be completely different 
# (the segmentation is the most challenging task, in my opinion).
# 
# both uses depth ratio and the “B allele” frequency, although 
# FACETS uses the log-transformed data for both.
# 
# ploidy/cellularity model seems similar, The main differences are in the 
# implementations (obviously) and also a quite different approach on the inference 
# of the parameters.
# 
# https://groups.google.com/forum/#!searchin/sequenza-user-group/Questions$20on$20BAF$20and$20LogR%7Csort:date/sequenza-user-group/LAeG7BPeHEo/F81bDDk1AgAJ




# run all binning and preproc steps
sequenza_cmds <- function(){
  
}

# a set of R funcs, callable by sequenza_cmds
if(FALSE){
  seqz = ''
  isMale = TRUE
  chroms = c(1:22, "X")
  cores = 10
  odir
  
  
  # test a extract:
  library(pacman)
  p_load(tidyverse, readr, sequenza)
  load("~/.rsrch2/iacs/iacs_dep/sseth/tmp/test_clone_conv/185_006/185_006_v1/185_006_v1_031420181521069040_sequenza_extract.RData")
  # v2 is with female
  odir = "~/.rsrch2/iacs/iacs_dep/sseth/tmp/test_clone_conv/185_006/185_006_v2"
  seqz_ext = `185_006_v1_031420181521069040_sequenza_extract`
  seqz_fit <- sequenza.fit(seqz_ext, female = T, chromosome.list = chroms, mc.cores = cores)
  seqz_res = sequenza.results(sequenza.extract = seqz_ext, 
                              cp.table = seqz_fit, 
                              sample.id = "185_006_v1", 
                              out.dir = odir, 
                              female = TRUE, chromosome.list = chroms)
  
  
}


# seqz ex --------
if(FALSE){
  library(pacman)
  p_load(sequenza)
  seqz0 = "~/.rsrch2/iacs/iacs_dep/sseth/tmp/test_clone_conv/185_006/mpileup_sequenza_temp/MS-PDX-300x-1728D-CB223_sahil-170713-MS-BREAST-PDX-SPANCER-1-NoIndex.bwa_recalibed.bam_seqz.gz"
  seqz1 = "~/.rsrch2/iacs/iacs_dep/sseth/tmp/test_clone_conv/185_006/mpileup_sequenza_temp/MS-PDX-300x-1728D-CB224_sahil-170713-MS-BREAST-PDX-SPANCER-1-NoIndex.bwa_recalibed.bam_seqz.gz"
  seqzs = "~/.rsrch2/iacs/iacs_dep/sseth/tmp/test_clone_conv/185_006/mpileup_sequenza_temp/MS-PDX-300x-1728D-CB225_sahil-170713-MS-BREAST-PDX-SPANCER-1-NoIndex.bwa_recalibed.bam_seqz.gz"
  
  chroms = c(1:22, "X")
  odir = "~/.rsrch2/iacs/iacs_dep/sseth/tmp/test_clone_conv/185_006/185_006_v2"
  seqz = seqz1;oprefix = "bam1"
  
}


# ** gc ------
# python_exe="/rsrch2/iacs/iacs_dep/sseth/apps/conda/2.7/bin/python"
# ref_fasta="/rsrch2/iacs/iacs_dep/sseth/ref/human/b37/fastas/Homo_sapiens_assembly19.fasta"
# sequenza_utils_exe="/rsrch2/iacs/iacs_dep/sseth/R/x86_64-pc-linux-gnu-library/3.5/sequenza/exec/sequenza-utils.py"
# need to do ONCE
# ${python_exe} ${sequenza_utils_exe} GC-windows -w 50 ${ref_fasta} | gzip > Homo_sapiens_assembly19.gc50Base.txt.gz

sequenza_seqz <- function(bam_t, 
                          bam_n, 
                          samplename,
                          sequenza_utils_exe = "/rsrch2/iacs/iacs_dep/sseth/R/x86_64-pc-linux-gnu-library/3.5/sequenza/exec/sequenza-utils.py",
                          oprefix){
  
  # facets uses the following:
  # pileupfl = glue("{oprefix}_pileup.tmp.gz")
  # snp_pipeup_param = "-g -q15 -Q20 -P100 -r25,0"
  # cmd0 = glue("{snp_pileup_exe} -g -q15 -Q20 -P100 -r25,0 {vcffl} {pileupfl} {bam_n} {bam_t}")
  # cmd0 = as.character(cmd0)
  
  # https://cran.r-project.org/web/packages/sequenza/vignettes/sequenza.pdf
  # 1 samtools mpileup −f hg19.fasta −Q 20 normal.bam | gzip > normal.pileup.gz
  # 2 samtools mpileup −f hg19.fasta −Q 20 tumor.bam | gzip > tumor.pileup.gz
  
  # need to do ONCE
  # sequenza−utils.py GC−windows −w 50 hg19.fa | gzip > hg19.gc50Base.txt.gz
  
  # sequenza−utils.py bam2seqz −gc hg19.gc50Base.txt.gz
  # −−fasta hg19.fasta \
  #  −n normal.fifo \
  # −t tumor.fifo | gzip > out.seqz.gz
  
  # sequenza−utils.py seqz−binning −w 50 \
  # −s out.seqz.gz | gzip > out small.seqz.gz
}

sequenza_r_pipe <- function(seqz, 
                            oprefix,
                            odir,
                            is_male = TRUE, 
                            chroms = c(1:22, "X"), 
                            cores = 10,
                            ....){
  
  # window: ONLY for plotting, window size controlled by breaks.method ONLY
  # https://groups.google.com/forum/#!searchin/sequenza-user-group/too$20many$20segments%7Csort:date/sequenza-user-group/Hnzswcq-34s/nZ94ehwAGQAJ
  # min.type.freq=0.9: to remove freq type (in case of multiple muts, at location)
  #     var ignored at: 5/6(0.83) 0/6(0) and 1/6(0.17)
  # weighted.mean=TRUE; eighted mean, positions with more reads have more weight
  #     in the calculation of the mean value for the segment
  #     (a segment value is an average of all the positions included in it)...
  #     In some case is better with the feature turned off in presence of bias
  #     (positions with huge amount of reads)
  seqz_ex <- sequenza.extract(seqz, chromosome.list = chroms, ...)

  # 
  #   N.ratio.filter=10; minimum no of obs, required to consider a segmetn
  #   The fit step is particularly important to filter out unreliable segments.
  #   This will produce the estimation, we are not throwing out segments for
  #   downstream analysis.
  # 
  #   N.BAF.filter=1; ONLY HETEROZYGOUS SNPS ARE USED
  #     So a segment could have 100 "ratio" pos, but only 2 "heterozygous" pos.
  #     I think you are confusing yourself with the variants in the segment in this case:
  #     the only variant that we consider for the segment are the germlines heterozygous,
  #     so position that are detected het already in the normal.
  # 
  #   ratio.priority: FALSE
  #    This is advisable if the BAF profile looks very noisy.
  #    In WES it could be the case, in WGS the heterozygous position are much more and more
  #    reliable, but in exome it depends a lot on the specific dataset.
  seqz_fit = sequenza.fit(seqz_ex, female = !is_male, 
                          chromosome.list = chroms,
                          ratio.priority = T,
                          N.ratio.filter = 20, # ask for more signal
                          mc.cores = cores)
  sequenza::cp.plot.contours(seqz_fit)
  plot(density(seqz_fit$ploidy))
  seqz_fit$cellularity
  
  seqz_res = sequenza.results(sequenza.extract = seqz_ex, 
                              cp.table = seqz_fit, 
                              sample.id = oprefix, 
                              out.dir = odir, 
                              female = !is_male, 
                              chromosome.list = chroms)
  
  # get segs
  ADRatio <- mean(extracted$gc$adj[, 2])
  seg_mat <- na.exclude(do.call("rbind", extracted$segments))
  cns <- baf.bays(BF = segMat$Bf, depth.ration = segMat$depth.ratio, 
                  cellularity = cellularity, ploid = ploidy, 
                  avg.depth.ratio = ADRatio)
  
  # plots
  
  
  # write out the stuff
  save_rds(fit, glue("{oprefix}_fit.rds"))
  save_rds(res, glue("{oprefix}_res.rds"))
  save_rds(extract, glue("{oprefix}_extract.rds"))
  save_rds(seg_mat, glue("{oprefix}_seg_mat.rds"))
  save_rds(cns, glue("{oprefix}_cns.rds"))
  
  
}


# getSeqzSegs <- function(seqz, cellularity, ploidy, odir, oprefix){
#   require(sequenza)
#   if(missing(odir)){
#     odir <- dirname(seqz)
#   }
#   if(missing(oprefix)){
#     oprefix <- gsub("_.*", "", basename(oprefix))
#   }
#   extracted <- sequenza.extract(seqz)
#   ADRatio <- mean(extracted$gc$adj[, 2])
#   segMat <- na.exclude(do.call("rbind", extracted$segments))
#   cns <- baf.bays(BF = segMat$Bf, depth.ration = segMat$depth.ratio, 
#                   cellularity = cellularity, ploid = ploidy, 
#                   avg.depth.ratio = ADRatio)
#   
#   return(cbind(segMat, cns))
# }




## segMat - returned by getSeqzSegs
plotABCopy <- function(segMat){
  genome.view(seg.cn = segMat, info.type = "CNt")
  legend("bottomright", bty = "n", "Overall copy number", col = "red", 
         inset = c(0, -0.4), pch = 15, xpd = TRUE)
  return(invisble())	
}

plotASCopy <- function(segMat){
  genome.view(seg.cn = segMat, info.type = "AB")
  legend("bottomright", bty = "n", c("A-allele", "B-allele"), 
         col = c("red", "blue"), inset = c(0, -0.45), pch = 15, xpd = TRUE)
  return(invisble())		
}

plotQA <- function(gcStats){ 
  gcStats <- gc.norm(seqzData$depth.ratio, gc = seqzData$GC.percent)	
  gcVect <- setname(gcStats$raw.mean, gc.stats$gc.value)
  seqzData$adjusted.ratio <- seqzData$depth.ratio/gc.vect[as.character(seqzData$GC.percent)]
  
}














# END

