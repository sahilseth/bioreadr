% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bam_set.R
\name{paired_bam_set}
\alias{paired_bam_set}
\title{Create a single object, to be used by multiple modules}
\usage{
paired_bam_set(
  tumor_bam,
  normal_bam,
  outprefix,
  is_merged,
  split_by_chr,
  ref_fasta_path = opts_flow$get("ref_fasta_path")
)
}
\arguments{
\item{tumor_bam}{tumor bam file}

\item{normal_bam}{normal bam file}

\item{outprefix}{output prefix for paired sample analysis, example tumorname___normalname}

\item{is_merged}{are bam files merged - single bam file for all chromosomes [TRUE/FALSE]}

\item{split_by_chr}{should downstream analysis be split by chr.}

\item{ref_fasta_path}{path to reference fasta file}
}
\description{
Create a single object, to be used by multiple modules
}
