#!/bin/bash
cd $1;

#files=$(ls */*/*.fastq.gz)
files=$(find . -name "*fastq.gz")

## ------- get summary per file
for file in $files; do
    echo "echo $file; gunzip -c $file | awk '{if(NR%4==2) tab[\$0]++} END {for (i in tab) {print i, tab[i]} }' > $(echo $file  | sed 's/fastq.gz/txt/g')"
done | parallel -j $cores

## ----- get total indexes
files=$(find . -name "*txt")
cat $files > total_index.txt

## --- get a sorted summary per table
gawk '{ind[$1]=ind[$1]+$2} END {sum=0; for(i in ind) sum=sum+ind[i]; for(i in ind) {print i, ind[i], ind[i]/sum*100}}' total_index.txt | sort -nrk2 - > total_summary.txt
