

get_chrs <- function(bam) {
    ## | grep -v MT |grep -v NC |grep -v GL
    cmd <- sprintf("samtools view -H %s | grep  '\\@SQ' | cut -f 2 | cut -d ':' -f 2", bam)
    chrs <- system(cmd,intern=TRUE)
    return (chrs)
}

get_chr_lens <-function(bam) {
    cmd <- sprintf("samtools view -H %s | grep  '\\@SQ' | cut -f 2,3 | cut -d ':' -f 3", bam)
    chrlens <- system(cmd,intern=TRUE)
    return (as.integer(chrlens))
}
