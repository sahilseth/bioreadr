#!/bin/bash -e


# usage: split_fq.sh number_of_chunks fastq
# use zcat OR cat to fetch the number of lines

#fq=/rsrch1/iacs/tmp/illumina_platinum/50x/NA12877/ERR194146_1.fastq.gz

## --- needs the following variables declared in name space
chunks=$1
fq=$2

if [ ! -f $fq ]; then
    echo "File not found!"
    exit 1
fi

echo `date` "; Working on $fq, to be split into $chunks chunks, getting total lines..."

## extenstion
if [ "${fq##*.}" == "gz" ]; then
  cat_cmd="zcat"
else
  cat_cmd="cat"
fi

# ------------------------ get number of lines ------------------------- #
tot_lines=$($cat_cmd $fq | wc -l )

# --------------------- num should be divisible by 4 ------------------- #
num=$(expr $tot_lines / $chunks)
num=$(expr $num % 4 + $num)

echo -e `date` "; The file would be split into $chunks chunks each with $num lines, splitting...."

# basename:
base="${fq%%.*}."


$cat_cmd $fq | split -l $num -d -a 3 - $base


#-d  adds digit prefixes instead of alphabetic
#-a 3  allows you to add more digits to output files.
#b  is the prefix to use
#In this example: output will be b000...b999


# ------------------------ rename files ------------------------------------ #
echo `date` "; Renaming files..."
end=$(expr $chunks - 1)
for i in $(seq 0 $end); do
  out_fls=$base$(printf "%03d" $i);
  mv $out_fls ${out_fls}.fastq
done


##~% FILE="example.tar.gz"
##~% echo "${FILE%%.*}"
##example
##~% echo "${FILE%.*}"
##example.tar
##~% echo "${FILE#*.}"
##tar.gz
##~% echo "${FILE##*.}"
##gz
