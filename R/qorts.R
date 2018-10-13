

'
bam="../bams/RNA-334187-T_rg_reorder.bam"
  outpath="RNA-1614-T"

java -Xmx20G -jar ${quorts_exe} QC \
--singleEnded --maxReadLength 400 --generatePlots --numThreads 10 \
--chromSizes ${chrom_sizes} --noGzipOutput\
--addFunctions annotatedSpliceExonCounts,makeWiggles,makeAllBrowserTracks
${bam} \
${gtf} ${outpath}

'


qorts <- function(bam, samplename = opts_flow$get("samplename"),
                  outdir,
                  java_exe = opts_flow$get("java_exe"),
                  java_mem = "-Xmx20G",
                  java_tmp = opts_flow$get("java_tmp"),

                  qorts_jar="/rsrch2/iacs/apps/qorts/0.3.18/QoRTs.jar",
                  chrom_sizes="/scratch/rists/hpcapps/reference/human/ucsc_hg19/fastas/ucsc.hg19.sizes",
                  gtf="/scratch/rists/hpcapps/reference/human/hg19/annotations/ensemble/Homo_sapiens.GRCh37.75_chr.gtf",
                  qorts_opts = "--singleEnded --maxReadLength 400 --generatePlots --numThreads 10 --noGzipOutput --addFunctions annotatedSpliceExonCounts,makeWiggles,makeAllBrowserTracks"
                  ){

  # for ion torrent, it is impertive to put the max read length:
  # --maxReadLength 400
  # --isSingleEnd
  # and remove --stranded, since it is not stranded
  # its better for reads to be: --nameSorted, not sure WHY!

  check_args()

  qorts = sprintf("%s %s -Djava.io.tmpdir=%s -jar %s QC %s --chromSizes %s %s %s %s",
                    java_exe, java_mem, java_tmp, qorts_jar, qorts_opts, chrom_sizes, bam, gtf, outdir)

  ## --- INPUT is a NAMED list
  flowmat = to_flowmat(list(qorts = qorts), samplename)

  return(list(flowmat = flowmat, outfiles = outdir))
}
