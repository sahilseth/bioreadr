# 
# 
# Pipeline Inputs
# 
# STAR, RSEM, and Kallisto indexes were all built with the same reference genome and annotation file.
# 
# Reference genome: HG38 (no alt analysis) with overlapping genes from the PAR locus removed (chrY:10,000-2,781,479 and chrY:56,887,902-57,217,415).
# ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_genbank/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh38/seqs_for_alignment_pipelines
# Annotation: Gencode V23 comprehensive annotation (CHR)
# http://www.gencodegenes.org/releases/23.html)
# Inputs were generated with the following commands
# 
# STAR
# sudo docker run -v $(pwd):/data quay.io/ucsc_cgl/star --runThreadN 32 --runMode genomeGenerate --genomeDir /data/genomeDir --genomeFastaFiles hg38.fa --sjdbGTFfile gencode.v23.annotation.gtf
# Kallisto
# sudo docker run -v $(pwd):/data quay.io/ucsc_cgl/kallisto index -i hg38.gencodeV23.transcripts.idx transcriptome_hg38_gencodev23.fasta
# RSEM
# sudo docker run -v $(pwd):/data --entrypoint=rsem-prepare-reference jvivian/rsem -p 4 --gtf gencode.v23.annotation.gtf hg38.fa hg38