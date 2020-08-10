#!/bin/env bash

# this merged two GATK interval type files with a header like this:
# echo "usage:
#   gatk_cnv_cnts_concat.sh fl1.tsv fl2.tsv fl_merged.tsv
#   "

#' @HD     VN:1.6
#' @SQ     SN:1    LN:249250621
#' @SQ     SN:2    LN:243199373
#' @SQ     SN:3    LN:198022430
#' @SQ     SN:4    LN:191154276
#' @SQ     SN:5    LN:180915260
#' @SQ     SN:6    LN:171115067
#' @SQ     SN:7    LN:159138663
#' @RG     ID:GATKCopyNumber       SM:1004-T
#' CONTIG  START   END     LOG2_COPY_RATIO
#' 1       861072  861643  -0.125668

awk 'FNR==1 && NR!=1 {{ while (/^@|^CONTIG/) getline; }} 1 {{ print }}' $1 $2 1> $3