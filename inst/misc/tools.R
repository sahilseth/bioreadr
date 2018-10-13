#### PICARD, logFile sample place as out
## tool="EstimateLibraryComplexity.jar";params="R=/IACS1/NGS/hg19BWAIndex/Homo_sapiens_assembly19.fasta"
runPicard <- function(bamFile,tool,params="",outsuf="picard",outFile,picardPath=getOption("ngs.picardPath"),
	  outPath=getOption("ngs.sampleQAPath"),java=getOption("ngs.java"),
	  javaMem=getOption("ngs.javaMem"),insuf=".bam|.sam",force=FALSE,verbose=TRUE,execute=TRUE){
    if(missing(outFile))	outFile <- file.path(outPath,gsub(insuf,sprintf(".%s.out",outsuf),basename(bamFile)))
    logFile <- file.path(outPath,gsub(insuf,sprintf(".%s.log",outsuf),basename(bamFile)))
    log <- paste(as.character(Sys.time()),"\t",bamFile,"\t","BGN",tool,"run\n")
    cat(log,file=logFile,append=TRUE); if(verbose){cat(log)}
    if(!force & file.exists(outFile))     return(outFile)
    cmd <- sprintf("time %s %s -Djava.io.tmpdir=/IACS1/tmp -jar %s/%s INPUT=%s OUTPUT=%s VALIDATION_STRINGENCY=LENIENT %s 2>>%s",
                   java,javaMem,picardPath,tool,bamFile,outFile,params,logFile)
    if(verbose) cat(cmd,"\n")
    if(execute) system(cmd)
    else return(cmd)
    log <- paste(as.character(Sys.time()),"\t",bamFile,"\t","END",tool,"run\n")
    cat(log,file=logFile,append=TRUE); if(verbose){cat(log)}
    return(outFile)
}

runPicard.BamIndexStats <- function(bamFile,tool="BamIndexStats.jar",
                                    params="",outsuf,picardPath=getOption("ngs.picardPath"),outPath=getOption("ngs.sampleQAPath"),
                                    java=getOption("ngs.java"),javaMem=getOption("ngs.javaMem"),insuf=".bam|.sam",
                                    prefix=gsub(".bam|.sam","",basename(bamFile)),
                                    force=FALSE,
                                    verbose=TRUE){
  outFile <- file.path(outPath,gsub(insuf,sprintf(".%s.out",outsuf),basename(bamFile)))
  logFile <- file.path(outPath,gsub(insuf,sprintf(".%s.log",outsuf),basename(bamFile)))
  log <- paste(as.character(Sys.time()),"\t",bamFile,"\t","BGN",tool,"run\n")
  cat(log,file=logFile,append=TRUE); if(verbose){cat(log)}
  if(!force & file.exists(outFile))     return(outFile)
  cmd <- sprintf("time %s %s -Djava.io.tmpdir=/IACS1/tmp -jar %s/%s INPUT=%s VALIDATION_STRINGENCY=LENIENT %s 1>>%s 2>>%s",
                 java,javaMem,picardPath,tool,bamFile,params,outFile,logFile)
  if(verbose) cat(cmd,"\n")
  system(cmd)
  log <- paste(as.character(Sys.time()),"\t",bamFile,"\t","END",tool,"run\n")
  cat(log,file=logFile,append=TRUE); if(verbose){cat(log)}
  return(outFile)
}

### http://picard.sourceforge.net/command-line-overview.shtml picard is designed for 2gb ram
runPicard.fixRGTags <- function(bam,tool="AddOrReplaceReadGroups.jar",
                                picardPath=getOption("ngs.picardPath"),
                                outBam,rgid,rgsm,rglb,rgpu,rgpl="Illumina_HiSeq2000",rgcn='IACS',
                                javaMem=getOption("ngs.javaMem"),
                                javaTemp="/IACS1/tmp",
                                verbose=FALSE, execute=FALSE){
    ## picard get the correct name
    cat("Working on ",bam,"\n")
    log=gsub(".bam",".log",outBam)
    cmd=sprintf("%s %s -Djava.io.tmpdir=%s -jar %s/%s INPUT=%s OUTPUT=%s SORT_ORDER=coordinate RGID=%s RGLB=%s RGPL=%s RGPU=%s RGSM=%s RGCN=%s VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true",
    getOption("ngs.java"),javaMem,javaTemp,picardPath,tool,bam,outBam,rgid,rglb,rgpl,rgpu,rgsm,rgcn)
    if(verbose) cat(cmd,"\n");
    if(execute) system(cmd)
}

##debug(runSamFlagStat)
runSamFlagStat <- function(bamFile,samtools=file.path(getOption("ngs.samtoolPath"),"samtools"),outsuf="flagstat",
                           outPath=getOption("ngs.sampleQAPath"),insuf=".bam|.sam",verbose=TRUE,force=FALSE){
    dir.create(outPath,recursive=TRUE)
    outFile <- file.path(outPath,gsub(insuf,sprintf(".%s",outsuf),basename(bamFile)))
    logFile <- file.path(outPath,gsub(insuf,sprintf(".%s.log",outsuf),basename(bamFile)))
    log <- paste(as.character(Sys.time()),"\t",bamFile,"\t","BGN FLAGSTAT run\n")
    if(!force & file.exists(outFile))     return(outFile)
    cat(log,file=logFile,append=TRUE); if(verbose){cat(log)}
    cmd <- paste("time ",samtools," flagstat ",bamFile," > ",outFile,sep="")
    if(verbose) cat(cmd,"\n")
    try(system(cmd))                  #spits out too much stuff...
    log <- paste(as.character(Sys.time()),"\t",bamFile,"\t","END FLAGSTAT run\n")
    cat(log,file=logFile,append=TRUE); if(verbose){cat(log)}
    return(outFile)
}


## tool="DepthOfCoverage";params="-L MT";outsuf="dopMT";outPath="/IACS1/home/sseth/projects/moon/dopMT"
## tools="UnifiedGenotyper"; params="-L MT -glm BOTH -dfrac 1 --sample_ploidy 2 --dbsnp /IACS1/NGS/bundle1.2_b37/dbsnp_132.b37.vcf"
## type="gatk.SomaticIndel";tool <- "SomaticIndelDetector";logPath="/IACS1/GCC/log/indel";opts <- ""
runGatk <- function(bamFile,tool,params="",outsuf,gatkPath=getOption("ngs.gatkPath"),outPath=getOption("ngs.sampleQAPath"),
                    logPath=getOption("ngs.sampleQAPath"),
                    java=getOption("ngs.java"),javaMem=getOption("ngs.javaMem"),insuf=".bam|.sam",verbose=TRUE,force=FALSE){
    ## out <- as.list(match.call(expand.dots = TRUE)[-1])
    ## print(out);out2=sapply(out,function(x) try(get(as.character(x))));print(out2)
    outFile <- file.path(outPath,gsub(insuf,sprintf(".%s.out",outsuf),basename(bamFile)))
    logFile <- file.path(logPath,gsub(insuf,sprintf(".%s.log",outsuf),basename(bamFile)))
    log <- paste(as.character(Sys.time()),"\t",bamFile,"\t","BGN",tool,"run\n")
    cat(log,file=logFile,append=TRUE); if(verbose){cat(log)}
    if(!force & file.exists(outFile)){
        if(verbose) cat(sprintf("The outfile already exists: %s Use some force!",outFile));
        return(outFile)}
    cmd <- sprintf("time %s %s -Djava.io.tmpdir=/IACS1/tmp -jar %s/GenomeAnalysisTK.jar -T %s -I %s -o %s %s 1>> %s 2>>%s",
                   java,javaMem,gatkPath,tool,bamFile,outFile,params,logFile,logFile)
    if(verbose) cat(cmd,"\n")
    try(system(cmd))
    log <- paste(as.character(Sys.time()),"\t",bamFile,"\t","END",tool,"run\n")
    cat(log,file=logFile,append=TRUE); if(verbose){cat(log)}
    return(outFile)
}


ngs.splitBamFileName <- function(bamFile){
    tmp <- strsplit(gsub("(.*)_([0-9]{6}.*)_s_([1-8]{1})_([TCGA]*).(.*)","\\1,\\2,\\3,\\4,\\5",basename(bamFile)),",")[[1]]
    sample <- tmp[1];runid <- tmp[2];lane <- tmp[3];index <- tmp[4]; ext <- tmp[5];merge.version <- "NA";
    if(length(grep("merged",bamFile))>0){
        ## runid <- basename(dirname(bamFile))
        tmp <- strsplit(gsub("(.*)_([0-9]{6}.*)_(.*)_(.*)","\\1,\\2,\\3,\\4",basename(bamFile)),",")[[1]]
        sample <- tmp[1];runid <- tmp[2];merge.version <- tmp[3];ext <- tmp[4]
        lane <- "NA";index="NA"
    }
    return(list(runid=runid,sample=sample,lane=lane,index=index,merge.version=merge.version,ext=ext))
}


if(FALSE){
    ## insSize
    pdfFile <- file.path(outPath,gsub(suf,sprintf("%s.pdf",outsuf),basename(bamFile)))
    param <- sprintf("HISTOGRAM_FILE=%s",pdfFile)
    outFile <- runPicard(bamFile,tool="CollectInsertSizeMetrics.jar",params=params,outsuf="isize",
                         picardPath=getOption("ngs.picardPath"),outPath=getOption("ngs.sampleQAPath"))
}


fastmean <- function(dat) {
    with(dat, sum(freq*value)/sum(freq) )
}
fastRMSE <- function(dat) {
    mu <- fastmean(dat)
    with(dat, sqrt(sum(freq*(value-mu)^2)/(sum(freq)-1) ) )
}
