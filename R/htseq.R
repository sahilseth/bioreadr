## in-house

## Transcript-level count is not accurate as htseqcount discards the reads mapped to different isoforms
## cmds <- c(cmds, sprintf("%s -i 'transcript_id' %s %s %s > %s", htseqcount, option, bam, gtf,
##                         file.path(outdir, paste(oprefix, source, "ensembl", "htseqcount.transcript.txt", sep="."))))


opts_flow$set(
  htseq_opts = "-f bam -m union -a 10 -s no",
  gtf_mm10 = "/scratch/rists/hpcapps/reference/mouse/mm10/annotations/ensembl/mm10.GRCm38.75.gtf",
  gtf_hg19_ensembl = "/scratch/rists/hpcapps/reference/human/broad_hg19/annotations/ensemble/Homo_sapiens.GRCh37.75.gtf",
  gtf_hg19_tcga = "/scratch/rists/hpcapps/reference/human/tcga/TCGA.hg19.June2011.gtf"
  )

"

"

#' A wrapper around the RNASeq counting tool, htseq
#'
#' @param x a bam file
#' @param samplename samplename
#' @param outprefix output prefix to be used
#' @param htseq_opts options to be supplied to htseq, default retrieved using opts_flow$get("htseq_opts").
#' @param gtf a GTF file used
#'
#' @export
#'
htseq <- function(bam, samplename,
                  outprefix,
                  htseq_exe = opts_flow$get("htseq_exe"),
                  htseq_opts = opts_flow$get("htseq_opts"),
                  gtf = opts_flow$get("gtf")){

  check_args()

  prefix_gene = paste(outprefix, "htseqcount.gene.txt", sep=".")
  prefix_exon = paste(outprefix, "htseqcount.exon.txt", sep=".")
  
  cmd_gene <- sprintf("%s %s %s %s > %s",
                  htseq_exe, htseq_opts, bam, gtf, prefix_gene)
  cmd_exon <- sprintf("%s -i 'exon_id' %s %s %s > %s",
                      htseq_exe, htseq_opts, bam, gtf, prefix_exon)
  ## Transcript-level count is not accurate as htseqcount discards the reads mapped to different isoforms
  ## -i 'transcript_id'

  cmds = list(htseq_gene = cmd_gene,
              htseq_exon = cmd_exon)

  flowmat = to_flowmat(cmds, samplename = samplename)

  return(list(flowmat = flowmat, outfile = c(prefix_gene, prefix_exon)))
}
