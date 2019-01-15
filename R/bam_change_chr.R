bam_change_chr <- function(bam,
                           outbam,
                           what = c("add", "remove"),
                           execute = FALSE,
                           samtool_exe = optsaturn$get("SaturnV_samtools")){
  what <- match.arg(what)
  
  message("> changing chr names for ", bam)
  
  if(what == "add"){
    if(missing(outbam))
      outbam <- gsub("bam$", "_chr.bam", bam)
    
    cmd <- paste(samtool_exe,  " view -H ", bam, " | sed -e 's/SN:1/SN:chr1/' | ",
                 "sed -e 's/SN:2/SN:chr2/' | sed -e 's/SN:3/SN:chr3/' | ",
                 "sed -e 's/SN:4/SN:chr4/' | sed -e 's/SN:5/SN:chr5/' | ",
                 "sed -e 's/SN:6/SN:chr6/' | sed -e 's/SN:7/SN:chr7/' | ",
                 "sed -e 's/SN:8/SN:chr8/' | sed -e 's/SN:9/SN:chr9/' | ",
                 "sed -e 's/SN:10/SN:chr10/' | sed -e 's/SN:11/SN:chr11/' | ",
                 "sed -e 's/SN:12/SN:chr12/' | sed -e 's/SN:13/SN:chr13/' | ",
                 "sed -e 's/SN:14/SN:chr14/' | sed -e 's/SN:15/SN:chr15/' | ",
                 "sed -e 's/SN:16/SN:chr16/' | sed -e 's/SN:17/SN:chr17/' | ",
                 "sed -e 's/SN:18/SN:chr18/' | sed -e 's/SN:19/SN:chr19/' | ",
                 "sed -e 's/SN:20/SN:chr20/' | sed -e 's/SN:21/SN:chr21/' | ",
                 "sed -e 's/SN:22/SN:chr22/' | sed -e 's/SN:X/SN:chrX/' | ",
                 "sed -e 's/SN:Y/SN:chrY/' | sed -e 's/SN:MT/SN:chrM/' | ",
                 "samtools reheader - ", bam, " > ", outbam, sep = "")
    
  }else if(what == "remove"){
    if(missing(outbam))
      outbam <- gsub("bam$", "nochr.bam", bam)
    
    cmd <- paste(samtool_exe,  " view -H ", bam, " | sed -e 's/SN:chr1/SN:1/' | ",
                 "sed -e 's/SN:chr2/SN:2/' | sed -e 's/SN:chr3/SN:3/' | ",
                 "sed -e 's/SN:chr4/SN:4/' | sed -e 's/SN:chr5/SN:5/' | ",
                 "sed -e 's/SN:chr6/SN:6/' | sed -e 's/SN:chr7/SN:7/' | ",
                 "sed -e 's/SN:chr8/SN:8/' | sed -e 's/SN:chr9/SN:9/' | ",
                 "sed -e 's/SN:chr10/SN:10/' | sed -e 's/SN:chr11/SN:11/' | ",
                 "sed -e 's/SN:chr12/SN:12/' | sed -e 's/SN:chr13/SN:13/' | ",
                 "sed -e 's/SN:chr14/SN:14/' | sed -e 's/SN:chr15/SN:15/' | ",
                 "sed -e 's/SN:chr16/SN:16/' | sed -e 's/SN:chr17/SN:17/' | ",
                 "sed -e 's/SN:chr18/SN:18/' | sed -e 's/SN:chr19/SN:19/' | ",
                 "sed -e 's/SN:chr20/SN:20/' | sed -e 's/SN:chr21/SN:21/' | ",
                 "sed -e 's/SN:chr22/SN:22/' | sed -e 's/SN:chrX/SN:X/' | ",
                 "sed -e 's/SN:chrY/SN:Y/' | sed -e 's/SN:chrM/SN:MT/' | ", 
                 "samtools reheader - ", bam, " > ", outbam, sep = "")
  }
  
  if(execute)
    system(cmd)
  else
    return(cmd)
  
  return(outbam)
}
