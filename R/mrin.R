#https://github.com/mskcc/RNAseqDB/blob/e431ca6610e4d03c19275d48e8054cfd61dcc68f/calc-expression.pl
# sub Run_mRIN {
#   my $bam = shift;
#   (-e $bam) or die "ERROR: $bam does not exist\n";
#   #`samtools view -b -F 1548 -q 30 $bam | $bedtools_dir/bamToBed -i stdin > sample.bed` if (!-e 'sample.bed');
#   if (!-e 'sample.bedGraph') {
#     `samtools view -b -F 1548 -q 30 $bam chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM | $bedtools_dir/bamToBed -i stdin > sample.bed` if (!-e 'sample.bed');
#     
#     `perl $mRIN_dir/tag2profile.pl -big -exact -of bedgraph -v sample.bed sample.bedGraph`;
#     `rm -f sample.bed` if (-e 'sample.bed' and -e 'sample.bedGraph');
#   }
#   `mkdir -p ks`;
#   #`perl ~/bin/mRIN/gen_transcript_cdf.pl -v ~/bin/mRIN/genes/refGene.rep.uniq.hg38.bed sample.bedGraph cdf/sample.cdf.bedGraph`;
#   #`perl ~/bin/mRIN/ks_test_uniform.pl    -v cdf/sample.cum.bedGraph ks/sample.ks.txt}`;
#   #`perl  /cbio/ski/schultz/home/wangq/bin/mRIN/gen_transcript_cdf.pl -v /cbio/ski/schultz/home/wangq/bin/mRIN/genes/refGene.rep.uniq.hg38.bed sample.bedGraph - | perl /cbio/ski/schultz/home/wangq/bin/mRIN/ks_test_uniform.pl -v - ks/sample.ks.txt`;
#   #`perl  /cbio/ski/schultz/home/wangq/bin/mRIN/gen_transcript_cdf.pl -v /cbio/ski/schultz/home/wangq/ref/transcript/gencode.v21/gencode.v21.annotation.bed sample.bedGraph - | perl /cbio/ski/schultz/home/wangq/bin/mRIN/ks_test_uniform.pl -v - ks/sample.ks.txt`;
#   `perl  $mRIN_dir/gen_transcript_cdf.pl -v $mRIN_dir/genes/refGene.rep.uniq.hg19.bed sample.bedGraph - | perl $mRIN_dir/ks_test_uniform.pl -v - ks/sample.ks.txt`;
#   `rm -f sample.bedGraph` if (-e 'sample.bedGraph' and -e 'ks/sample.ks.txt');
# }
