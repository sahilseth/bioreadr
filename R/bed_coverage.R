

# http://www.gettinggeneticsdone.com/2014/03/visualize-coverage-exome-targeted-ngs-bedtools.html
# https://www.ampliseq.com/help/faq.action

# example exome cnv run
# funr ExomeLyzer::procExome \
# tumorBAM=/rsrch2/iacs/ngs_runs/1411_sarcomatoid/tmap/bams/WEX-1676-T_nochr.bam \
# normalBAM=/rsrch2/iacs/ngs_runs/1411_sarcomatoid/tmap/bams/WEX-1676-N_nochr.bam \
# bedName=/rsrch2/iacs/iacs_dep/sseth/projects2/jz_sarco/metadata/AmpliSeqExome.20140814/AmpliSeqExome.20131001.submitted_nochr.bed \
# odir=/rsrch2/iacs/ngs_runs/1411_sarcomatoid/exomecnv oprefix=WEX-1676-T__WEX-1676-N cores=2 MAPQ=10 flanking=250

# example sequenza run
# funr ExomePurity:::runSequenza sampleBam=/rsrch2/iacs/ngs_runs/1411_sarcomatoid/tmap/bams/WEX-973-T_nochr.bam \
# refBam=/rsrch2/iacs/ngs_runs/1411_sarcomatoid/tmap/bams/WEX-973-N_nochr.bam \
# bedFile=/rsrch2/iacs/iacs_dep/sseth/projects2/jz_sarco/metadata/AmpliSeqExome.20140814/AmpliSeqExome.20131001.submitted_nochr.bed \
# odir=/rsrch2/iacs/ngs_runs/1411_sarcomatoid/sequenza oprefix=WEX-973-T__WEX-973-N \
# fa=/scratch/rists/hpcapps/reference/human/broad_hg19/fastas/Homo_sapiens_assembly19.fasta

# bed=~/projects2/jz_sarco/metadata/CCP.20140814/CCP.20131001.designed.bed
# find CCP*bam | parallel 'bedtools coverage -hist -abam {} -b $bed | grep ^all > {}.hist.all.txt'
# example
# bedtools coverage -hist -abam CCP-973-T.bam -b $bed

coverage <- function(){
  
}