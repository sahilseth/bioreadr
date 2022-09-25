


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

