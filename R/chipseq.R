## in-house


## chipseq pipeline

## By Sahil Seth, Samir Amin
## Contributions from Kadir, Ayush

# require(flowr)
# require(params)
# require(ultraseq)


## use following as a example
if(FALSE){
  
  
  ## creating a sheet with all the samples in this folder
  library(drat)
  drat::addRepo("sahilseth")      
  install.packages("flowr")
  install.packages("ultraseq")
  
  fqdir = "/rsrch2/iacs/ngs_runs/150806_SN1120_0348_BC79KEACXX/fastqs/Y76I8n1Y76/Project_Pancreatic-MS132"
  fqsheet = create_fq_sheet(x=fqdir)
  ## for tesitng, getting the first samples:
  s = unique(fqsheet$samplename)[1]
  fqs = subset(fqsheet, samplename == s)$file
  
  ## we get two output
  out = chipseq(fqs = fqs, sample_name = s, out_path = fqdir)
  
  ## submit jobs to SHARK
  fobj2 = submit_flow(out$fobj, execute = TRUE)
  
}

#' A complete chipseq pipeline from fastq to MACS and Scripture
#' 
#' @param fqs a character vector of paths to fastqs for a sample
#' @param samplename name of this sample (a character vector, length one.)
chipseq <-function(fqs,
                   samplename = "sample_name",
                   out_path,
                   deffile = system.file("data/chipseq.def", package = "SaturnV"),
                   
                   bowtie_exe = "/risapps/rhel6/bowtie/1.1.2/bowtie",
                   java_exe = "/risapps/noarch/jdk/jdk1.8.0_45/bin/java",
                   igvtools_exe = "/rsrch2/genomic_med/krai/apps/igvtools/default/igvtools.jar",
                   scripture_exe = "/rsrch2/genomic_med/krai/apps/scripture/scripture-beta2.jar",
                   macs14_exe = "macs14", 
                   
                   bowtie_index = "/rsrch2/genomic_med/krai/annot/indexes/bowtie/hg19",
                   hg19table = "/rsrch2/genomic_med/krai/annot/hg19.table",
                   chr_prefix = "",
                   peak_fdr = "0.05"){
  
  ## check ALL arguments none of them should be NULL
  check_args()
  
  cmd_unzip <- sprintf("touch %s.running; gunzip -c %s > %s.fastq", samplename, paste(fqs, collapse = " "), samplename) 
  
  cmd_align <- sprintf("%s -n 1 -m 1 -S --best --strata -p 4 --chunkmbs 320 -t %s %s.fastq %s.sam",
                       bowtie_exe, bowtie_index, samplename, samplename)
  
  cmd_bam <- sprintf("module load samtools; samtools view -bS %s.sam > %s.bam",
                     samplename, samplename)
  
  ## primary output file
  bam = sprintf("%s.sorted.bam", samplename)
  
  cmd_sort <- sprintf("%s -m 3072M -T %s.sorted -o %s %s.bam",
                      samtools_exe, samplename, bam, samplename)
  
  cmd_index <- sprintf("%s index %s",
                       samtools_exe, bam)
  
  cmd_flagstat <- sprintf("%s flagstat %s > %s.flagstat",
                          samtools_exe, bam, bam)
  
  cmd_bam2bed <- sprintf("bedtools bamtobed -i %s > %s.sorted.bed",
                         bam, samplename)
  
  cmd_maketdf <- sprintf("%s -Xmx2g -jar %s count -w 25 -e 164 %s %s.tdf hg19",
                         java_exe, igvtools_exe, bam, samplename)
  
  cmd_macs14 <- sprintf("%s -t %s --nomodel -n %s",
                        macs14_exe, bam, samplename)
  
  
  ## scripture will have 25 commands, one for each chr
  scripture()
  
  
  
  ## creating a named list
  ## each element of the list can be of different length, to_flowmat can handle that.
  cmd_list = list(unzip = cmd_unzip, 
                  align = cmd_align,
                  bam = cmd_bam,
                  sort = cmd_sort, 
                  index = cmd_index,
                  flagstat = cmd_flagstat,
                  bam2bed = cmd_bam2bed,
                  maketdf = cmd_maketdf, 
                  macs = cmd_macs14)
                  #hg19tab = cmd_hg19table, 
                  #scripture = cmd_scripture, 
                  #scripture_bed = cmd_scripture_bed)
  
  ## give a named list and get a data.frame back !
  flowmat <- to_flowmat(x = cmd_list, samplename = samplename)
  
  #def = as.flowdef(deffile)
  ## create a flow object
  #fobj = to_flow(flowmat, def, flowname = "chipseq", flow_run_path = out_path)
  
  invisible(list(flowmat = flowmat))
}
