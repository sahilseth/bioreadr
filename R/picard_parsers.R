## in-house


## outFile <- "/IACS1/GCC/QA/130306_SN208_0456_AC1REFACXX/MOON-LAL-AL603A666M8_1_TATCAG/MOON-LAL-AL603A666M8_130306_SN208_0456_AC1REFACXX_s_1_TATCAG.rg.sorted.estLib.out"

# "insert_size\tAll_Reads.fr_count\tAll_Reads.tandem_count"

#' @export
parse_est_libsize <- function(outFile){
  tmp=try(scan(outFile,what="character",sep="\n")) #read the output
  if(class(tmp)=="try-error"){return(list(libsize="NA"))} #at times when file is missing
  tab=try(as.matrix(do.call(rbind,strsplit(tmp[-(1:grep("duplication_group_count\t",tmp))],"\t"))))
  if(class(tab)=="try-error"){return(list(mrq="NA"))} #at times the table is empty
  tab=apply(tab,2,as.numeric);colnames(tab)=c("count","freq")
  return(tab)
}

## outFile <- "/IACS1/GCC/QA/130306_SN208_0456_AC1REFACXX/MOON-LAL-AL603A666M8_1_TATCAG/MOON-LAL-AL603A666M8_130306_SN208_0456_AC1REFACXX_s_1_TATCAG.rg.sorted.isize.out"

#  [9] "insert_size\tAll_Reads.fr_count"
#  [9] "insert_size\tAll_Reads.fr_count", "all_Reads.tandem_count"


#' @import magrittr
#' @export
parse_insert_size <- function(outFile){
  
  tmp <- try(scan(outFile, what="character",sep="\n", quiet = TRUE))

  out = c(med = NA, mu = NA, sd = NA, sd_trimmed = NA, tandem_per = NA)
  if(class(tmp)=="try-error") 
    return(out)
  
  # read the file, and get rows after insert_size
  start_from = grep("insert_size", tmp)
  tab = tmp[-c(1:start_from)] %>%
    strsplit("\t") %>%
    do.call(rbind, .) %>%
    apply(2, as.numeric) %>% as.data.frame()
  #colnames(tab2) = c("insert_size","All_Reads.fr_count")
  colnames(tab) = strsplit(tmp[start_from], "\t")[[1]]
  
  # get the summary:
  # median
  out['med'] = as.numeric(tab[which(sum(tab[,"All_Reads.fr_count"])*0.5 <= cumsum(tab[,"All_Reads.fr_count"]))[1], 1])
  out['mu'] <- dplyr::select(tab, freq = All_Reads.fr_count, value = insert_size) %>%
    fastmean() %>% round(2)
  out['sd'] <- dplyr::select(tab, freq = All_Reads.fr_count, value = insert_size) %>%
    fastRMSE() %>% round(2)

  # trimmed insert size SD, removing those above 500
  rows <- which(tab[,"insert_size"] <= 500 )
  out['sd_trimmed'] <- dplyr::select(tab[rows, ], freq = All_Reads.fr_count, value = insert_size) %>%
    fastRMSE() %>% round(2)

  # why??
  if(ncol(tab) > 2)
    out['tandem_per'] = round(sum(tab$All_Reads.tandem_count)/(sum(tab$All_Reads.tandem_count) + sum(tab$All_Reads.fr_count)), 3)

  return(out) 
}

## outFile <- "/IACS1/GCC/LevelII/130308_SN1222_0177_BD1UT7ACXX/MOON-LAL-AL639A739M9_130308_SN1222_0177_BD1UT7ACXX_s_6_AGTCAC.rg.sorted.marked.metrics"
## picard.parse.MarkDup(outFile)

#' @export
parse_mark_dup <- function(outFile, names.prefix=""){
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

#' @export
parse_hsmetrics <- function(outFile,names.prefix=""){
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


#file = "/rsrch2/iacs/tmp/151223_SN1120_0379_AC82C1ACXX/bams/KTakahashi-AZABMsubp-AfMDSLTHSC-12-T_C82C1ACXX-1-ACAAGCTA.bwa_reheaded.qual.arq.out"

#' @export
parse_read_qual <- function(file){
  df = params::read_sheet(file, ext = "tsv")
  colnames(df) = tolower(colnames(df))
  df
  
}

  
#' 
#' @export
#' 
parse_idxstats <- function(file){
   tmp=scan(file,what="character",sep="\n")
    tmp=tmp[-length(tmp)]
    ret <- try(do.call(rbind,lapply(1:length(tmp),function(i){
      strsplit(tmp[i],"\t|=| ")[[1]][c(1,4,7,10)]
      })))
   return(ret)
}

##ngs.getMask(runid); to setup read_length
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

#' Title
#'
#' @param x
#' @param genome_length
#' @param read_length
#' @param names.prefix
#'
#' @export
#'
parse_flagstat <- function(x,
                           genome_length = opts_flow$get("genome_length"),
                           read_length = opts_flow$get("read_length"),
                           verbose = opts_flow$get("verbose")){

  as.n=as.numeric

  if(is.null(genome_length))
    genome_length = 3^8

  check_args()

  read_length = as.n(read_length)

  ## GET total, mapped, properly paired, and itself
  f <- try(readLines(x))

  if(verbose > 1)
    print(f)


  totalReads <- as.n(strsplit(f[grep("total",f)[1]]," ")[[1]][1])
  cov <- round((totalReads*read_length)/genome_length, 2) # coverage
  #print(cov)
  mapped <- round(as.n(strsplit(f[grep("mapped",f)[1]],"\\+| ",perl=TRUE)[[1]][1])/as.n(totalReads) * 100,2)
  paired <- round(as.n(strsplit(f[grep("paired",f)[1]],"\\+| ",perl=TRUE)[[1]][1])/as.n(totalReads) * 100,2)
  read1 <- round(as.n(strsplit(f[grep("read1",f)[1]],"\\+| ",perl=TRUE)[[1]][1])/as.n(totalReads) * 100,2)
  read2 <- round(as.n(strsplit(f[grep("read2",f)[1]],"\\+| ",perl=TRUE)[[1]][1])/as.n(totalReads) * 100,2)
  pp <- round(as.n(strsplit(f[grep("properly",f)[1]],"\\+| ",perl=TRUE)[[1]][1])/as.n(totalReads) * 100,2)
  bothmapped <- round(as.n(strsplit(f[grep("itself",f)[1]],"\\+| ",perl=TRUE)[[1]][1])/as.n(totalReads) * 100,2)
  self <- round(as.numeric(strsplit(f[grep("itself",f)[1]], " ")[[1]][1])/as.numeric(totalReads) * 100, 2)
  singletons <- round(as.n(strsplit(f[grep("singletons",f)[1]],"\\+| ",perl=TRUE)[[1]][1])/as.n(totalReads) * 100,2)
  matetodifferentchr <- round(as.n(strsplit(f[grep("with mate mapped to a different chr",
                                                          f)[1]],"\\+| ",perl=TRUE)[[1]][1])/as.n(totalReads) * 100,2)
  out <- c(approx_cov = cov,
           total_reads = totalReads,
           mapped = mapped,
           paired = paired,
           read1 = read1,
           read2 = read2,
           properly_paired = pp,
           both_mapped = bothmapped,
           singletons = singletons,
           mate_on_diff_chr = matetodifferentchr)

  return(out)
}


fastmean <- function(dat) {
  with(dat, sum(freq*value)/sum(freq) )
}
fastRMSE <- function(dat) {
  mu <- fastmean(dat)
  with(dat, sqrt(sum(freq*(value-mu)^2)/(sum(freq)-1) ) )
}
