## outFile <- "/IACS1/GCC/QA/130306_SN208_0456_AC1REFACXX/MOON-LAL-AL603A666M8_1_TATCAG/MOON-LAL-AL603A666M8_130306_SN208_0456_AC1REFACXX_s_1_TATCAG.rg.sorted.estLib.out"
picard.parse.Libsize <- function(outFile){
  tmp=try(scan(outFile,what="character",sep="\n")) #read the output
  if(class(tmp)=="try-error"){return(list(libsize="NA"))} #at times when file is missing
  tab=try(as.matrix(do.call(rbind,strsplit(tmp[-(1:grep("duplication_group_count\t",tmp))],"\t"))))
  if(class(tab)=="try-error"){return(list(mrq="NA"))} #at times the table is empty
  tab=apply(tab,2,as.numeric);colnames(tab)=c("count","freq")
  return(tab)
}

## outFile <- "/IACS1/GCC/QA/130306_SN208_0456_AC1REFACXX/MOON-LAL-AL603A666M8_1_TATCAG/MOON-LAL-AL603A666M8_130306_SN208_0456_AC1REFACXX_s_1_TATCAG.rg.sorted.isize.out"
picard.parse.InsSize <- function(outFile){
  tmp <- try(scan(outFile,what="character",sep="\n"))
  if(class(tmp)=="try-error"){return(list(m="NA",sd="NA",sd2="NA"))}
  tab=as.matrix(do.call(rbind,strsplit(tmp[-(1:grep("insert_size",tmp))],"\t")))## read the file, and get rows after insert_size
  tab=apply(tab,2,as.numeric);tab2 <- tab[,1:2];colnames(tab2)=c("iSize","freq")
  ## get the summary:
  m=tab[which(sum(tab2[,"freq"])*0.5 <= cumsum(tab2[,"freq"]))[1],1] #median
  mu <- fastmean(data.frame(freq=tab2[,"freq"],value=tab2[,"iSize"]))
  s <- round(fastRMSE(data.frame(freq=tab2[,"freq"],value=tab2[,"iSize"])),2)
  n <- which(tab2[,"iSize"]<=500 )
  s1 <- round(fastRMSE(data.frame(freq=tab2[n,"freq"],value=tab2[n,"iSize"])),2)
  if(ncol(tab)>2){m <- paste(m,"E",sep="")}
  return(list(m=m,sd=s,sd2=s1))                   #sd excluding the outliers
}

## outFile <- "/IACS1/GCC/LevelII/130308_SN1222_0177_BD1UT7ACXX/MOON-LAL-AL639A739M9_130308_SN1222_0177_BD1UT7ACXX_s_6_AGTCAC.rg.sorted.marked.metrics"
## picard.parse.MarkDup(outFile)
picard.parse.MarkDup <- function(outFile,names.prefix=""){
  tmp <- try(scan(outFile,what="character",sep="\n"))
  names <- c('unpaired_reads_examined','read_pairs_examined','unmapped_reads',
             'unpaired_read_duplicates','read_pair_duplicates','read_pair_optical_duplicates',
             'percent_duplication','estimated_library_size')
  names <- paste(names.prefix,names,sep="")
  if(class(tmp)=="try-error"){
    values <- rep('NA',length(names));names(values) <- names
    return(values)
  }
  len <- grep("LIBRARY\tUNPAIRED_READS_EXAMINED",tmp)
  cols=strsplit(tmp[len],"\t")[[1]]
  vals=strsplit(tmp[len+1],"\t")[[1]][-1]
  pct.cols <- grep("percent",names);vals[pct.cols] <- as.numeric(vals[pct.cols])*100
  names(vals) <- names
  return(vals)
}

## outFile <- "/IACS1/GCC/QA/130426_SN1120_0250_BD24WYACXX/MOON-LAL-AL745A953S5_8_GCCAAG/MOON-LAL-AL745A953S5_130426_SN1120_0250_BD24WYACXX_s_8_GCCAAG.rg.sorted.recalibed.hsMetrics.out"
picard.parse.hsMetric <- function(outFile,names.prefix=""){
  as.n=as.numeric
  tmp <- try(scan(outFile,what="character",sep="\n"))
  names <- c('genome_size','bait_territory','target_territory','bait_design_efficiency','total_reads','pf_reads',
             'pf_unique_reads','pct_pf_reads','pct_pf_uq_reads','pf_uq_reads_aligned','pct_pf_uq_reads_aligned',
             'pf_uq_bases_aligned','on_bait_bases','near_bait_bases','off_bait_bases','on_target_bases',
             'pct_selected_bases','pct_off_bait','on_bait_vs_selected','mean_bait_coverage','mean_target_coverage',
             'pct_usable_bases_on_bait','pct_usable_bases_on_target','fold_enrichment','zero_cvg_targets_pct',
             'fold_80_base_penalty','pct_target_bases_2x','pct_target_bases_10x','pct_target_bases_20x',
             'pct_target_bases_30x','hs_library_size','hs_penalty_10x','hs_penalty_20x','hs_penalty_30x',
             'at_dropout','gc_dropout','sample','library','read_group','pcttargetetX25','pcttargetetX50','pcttargetetX75',
             'pcttargetetX100')
  names <- paste(names.prefix,names,sep="")
  if(class(tmp)=="try-error"){
    values <- rep('NA',length(names));names(values) <- names
    return(values)
  }
  len <- grep("METRICS CLASS\tnet.sf.picard.analysis",tmp)+1 #rows to skip
  cols=strsplit(tmp[len],"\t")[[1]]
  vals=round(as.n(strsplit(tmp[len+1],"\t")[[1]][2:(length(names)-3)]),4)
  outFile2 <- gsub("hsMetrics.out","hsMetrics.perTarget.out",outFile)
  mean.cov <- try(system(sprintf("awk '{if(NR==1){for(i=1;i<=NF;i++){if($i==\"mean_coverage\"){COL=i;}}}else{print $COL;}}' %s",outFile2),intern=TRUE),silent=TRUE)
  covs <- round(sapply(c(25,50,75,100),function(x) sum(as.numeric(mean.cov) > x)/length(mean.cov)),4)
  vals <- c(vals,covs)
  names(vals) <- names
  ## multiply columns with pct by 100
  pct.cols <- grep("pct",names);vals[pct.cols] <- vals[pct.cols]*100
  cbind(as.character(vals),names)
  return(vals)
}


## /IACS1/GCC/QA/130306_SN208_0456_AC1REFACXX/MOON-LAL-AL622A706S3_5_TAGGAT/MOON-LAL-AL622A706S3_130306_SN208_0456_AC1REFACXX_s_5_TAGGAT.rg.sorted.recalibed.readCounts.tsv
getReadsOnTarget <- function(outFile){
    out2 <- try(read.table(outFile,header=TRUE),silent=TRUE)
}

picard.parse.IdxStats <- function(file){
   tmp=scan(file,what="character",sep="\n")
    tmp=tmp[-length(tmp)]
    ret <- try(do.call(rbind,lapply(1:length(tmp),function(i){
      strsplit(tmp[i],"\t|=| ")[[1]][c(1,4,7,10)]
      })))
   return(ret)
}

##ngs.getMask(runid); to setup seqLen
## 3194780 + 0 in total (QC-passed reads + QC-failed reads)
## 860 + 0 duplicates
## 1632467 + 0 mapped (51.10%:nan%)
## 3194780 + 0 paired in sequencing
## 1597390 + 0 read1
## 1597390 + 0 read2
## 1408008 + 0 properly paired (44.07%:nan%)
## 1472400 + 0 with itself and mate mapped
## 160067 + 0 singletons (5.01%:nan%)
## 74262 + 0 with mate mapped to a different chr
## 21821 + 0 with mate mapped to a different chr (mapQ>=5)
##out=lapply(flagstats,ngs.parseFlagstat);write.csv(cbind(bams,mat),file="/speedB/MiSeqData/Moon_QC_batch1/flagstat.csv",quote=FALSE)
samtools.parseFlagstat <- function(flagstatFile,genomeLen=as.numeric(getOption("ngs.genomeLen")),
                                   seqLen=as.numeric(getOption("ngs.seqLen")),names.prefix=""){
  ## GET total, mapped, properly paired, and itself
  flagstat <- try(readLines(flagstatFile))
  as.n=as.numeric
  names <- c("cov","totalreads", "mapped","paired","read1","read2",
             "properlypaired","bothmapped","singletons","matetodifferentchr")
  names <- paste(names.prefix,names,sep="")
  if(class(flagstat)=="try-error"){
    values <- rep('NA',length(names));names(values) <- names
    return(values)
  }
  totalReads <- as.n(strsplit(flagstat[grep("total",flagstat)[1]]," ")[[1]][1])
  cov <- round((totalReads*seqLen)/genomeLen ,2) # coverage
  mapped <- round(as.n(strsplit(flagstat[grep("mapped",flagstat)[1]],"\\+| ",perl=TRUE)[[1]][1])/as.n(totalReads) * 100,2)
  paired <- round(as.n(strsplit(flagstat[grep("paired",flagstat)[1]],"\\+| ",perl=TRUE)[[1]][1])/as.n(totalReads) * 100,2)
  read1 <- round(as.n(strsplit(flagstat[grep("read1",flagstat)[1]],"\\+| ",perl=TRUE)[[1]][1])/as.n(totalReads) * 100,2)
  read2 <- round(as.n(strsplit(flagstat[grep("read2",flagstat)[1]],"\\+| ",perl=TRUE)[[1]][1])/as.n(totalReads) * 100,2)
  pp <- round(as.n(strsplit(flagstat[grep("properly",flagstat)[1]],"\\+| ",perl=TRUE)[[1]][1])/as.n(totalReads) * 100,2)
  bothmapped <- round(as.n(strsplit(flagstat[grep("itself",flagstat)[1]],"\\+| ",perl=TRUE)[[1]][1])/as.n(totalReads) * 100,2)
  self <- round(as.numeric(strsplit(flagstat[grep("itself",flagstat)[1]], " ")[[1]][1])/as.numeric(totalReads) * 100, 2)
  singletons <- round(as.n(strsplit(flagstat[grep("singletons",flagstat)[1]],"\\+| ",perl=TRUE)[[1]][1])/as.n(totalReads) * 100,2)
  matetodifferentchr <- round(as.n(strsplit(flagstat[grep("with mate mapped to a different chr",
                                                          flagstat)[1]],"\\+| ",perl=TRUE)[[1]][1])/as.n(totalReads) * 100,2)
  out <- c(cov=cov,totalreads=totalReads, mapped=mapped,paired=paired,read1=read1,read2=read2,pp=pp,bothmapped=bothmapped,
           singletons=singletons,matetodifferentchr=matetodifferentchr)
  names(out) <- names
  return(out)
}
