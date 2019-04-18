"Support for Copy Number Variations (CNVs) with GATK4
https://software.broadinstitute.org/gatk/documentation/article?id=11682
https://gatkforums.broadinstitute.org/dsde/discussion/11683/

https://gatkforums.broadinstitute.org/gatk/discussion/7387/description-and-examples-of-the-steps-in-the-acnv-case-workflow
"

# prepare intervals:

# 1. Collect raw counts data with PreprocessIntervals and CollectFragmentCounts
# Before collecting coverage counts that forms the basis of copy number variant detection, we define the resolution of the analysis with a genomic intervals list. The extent of genomic coverage and the size of genomic intervals in the intervals list factor towards resolution.
# 
# Preparing a genomic intervals list is necessary whether an analysis is on targeted exome data or whole genome data. In the case of exome data, we pad the target regions of the capture kit. In the case of whole genome data, we divide the reference genome into equally sized intervals or bins. In either case, we use PreprocessIntervals to prepare the intervals list.
# 
# For the tutorial exome data, we provide the capture kit target regions in 1-based intervals and set --bin-length to zero.

# gatk PreprocessIntervals \
# -L targets_C.interval_list \
# -R /gatk/ref/Homo_sapiens_assembly38.fasta \
# --bin-length 0 \
# --interval-merging-rule OVERLAPPING_ONLY \
# -O sandbox/targets_C.preprocessed.interval_list




# ![](https://us.v-cdn.net/5019796/uploads/editor/3z/gim58s5j2wmk.png)

# cnv-C: panel of normals for cnv really helps
# this is a critical part here
# used for 

# Step 1. Het Pulldown
# ** These instructions describe one method for Het Pulldown for matched samples. For more options, including tumor-only, please see: http://gatkforums.broadinstitute.org/gatk/discussion/7719/overview-of-getbayesianhetcoverage-for-heterozygous-snp-calling **
#   
#   Inputs
# control_bam -- BAM file for control sample (normal).
# case_bam -- BAM file for case sample (tumor).
# reference_sequence -- FASTA file for b37 reference.
# snp_file -- Picard interval list of common SNP sites at which to test for heterozygosity in the control sample .
# Outputs
# normal_het_pulldown -- TSV file with M entries containing ref/alt counts, ref/alt bases, etc., where M is the number of hets called in the control sample.
# tumor_het_pulldown -- TSV file with M entries containing ref/alt counts, ref/alt bases, etc. for sites in the case sample that were called as het in the control sample, where M is the number of hets called in the control sample.
# Format for both output files:
#   
#   CONTIG  POSITION        REF_COUNT       ALT_COUNT       REF_NUCLEOTIDE  ALT_NUCLEOTIDE  READ_DEPTH
# 1       809876  5       16      A       G       21
# 1       881627  23      12      G       A       35
# 1       882033  9       10      G       A       19
# 1       900505  26      24      G       C       50
# ....snip....
# Invocation
# java -jar <path_to_gatk_protected_jar> GetBayesianHetCoverage --reference <reference_sequence>
#   --snpIntervals <snp_file> --tumor <case_bam> --tumorHets <tumor_het_pulldown> --normal <control_bam>
#   --normalHets <normal_het_pulldown> --hetCallingStringency 30

# Step 2. Allelic CNV
# Inputs
# tumor_het_pulldown -- Generated in step 1.
# coverage_profile -- Tangent-normalized coverage TSV file obtained in the GATK CNV case workflow.
# called_segments -- Called-segments TSV file obtained in the GATK CNV case workflow.
# output_prefix -- Path and file prefix for creating the output files. For example, /home/lichtens/my_acnv_output/sample1
# Outputs
# acnv_segments -- TSV file with name ending with -sim-final.seg containing posterior summary statistics for log_2 copy ratio and minor-allele fraction in each segment. Using the above output_prefix, /home/lichtens/my_acnv_output/sample1-sim-final.seg
# acnv_cr_parameters -- TSV file with name ending with -sim-final.cr.param containing posterior summary statistics for global parameters of the copy-ratio model. Using the above output_prefix, /home/lichtens/my_acnv_output/sample1-sim-final.cr.param
# acnv_af_parameters -- TSV file with name ending with -sim-final.af.param containing posterior summary statistics for global parameters of the allele-fraction model. Using the above output_prefix, /home/lichtens/my_acnv_output/sample1-sim-final.af.param
# Other files containing intermediate results of the calculation are also generated.
# 
# Invocation
# java -Xmx8g -jar <path_to_gatk_protected_jar> AllelicCNV  --tumorHets <tumor_het_pulldown>
#   --tangentNormalized <coverage_profile> --segments <called_segments> --outputPrefix <output_prefix>
#   Step 3. Call CNLoH and Balanced Segments
# ** WARNING: This tool is experimental and exists primarily for internal Broad use. **
#   
#   Inputs
# tumor_het_pulldown -- Generated in step 1.
# acnv_segments -- Generated in step 2 (*-sim-final.seg).
# coverage_profile -- Tangent-normalized coverage TSV file obtained in the GATK CNV case workflow
# output_dir -- Directory for creating the output files. For example, /home/lichtens/my_acnv_cnlohcalls_output/
#   Outputs
# GATK-CNV-formatted seg file -- TSV file ending with -sim-final.cnv.seg. This file is formatted identically as the output of GATK CNV. Note that this implies that the allelic fraction values are not captured in this file.
# AllelicCapSeg-formatted seg file -- TSV file ending with -sim-final.acs.seg. This file is formatted identically as the output of Broad CGA AllelicCapSeg. Note that this file can be used as input to Broad-internal versions of ABSOLUTE.
# TITAN-compatible het file --TSV file ending with -sim-final.titan.het.tsv. This file can be used as the input to TITAN for the het read counts.
# TITAN-compatible copy-ratio file -- TSV file ending with -sim-final.titan.tn.tsv. This file can be used as the input to TITAN for the per-target copy-ratio estimates.
# Invocation
# java -Xmx8g -jar <path_to_gatk_protected_jar> CallCNLoHAndSplits  --tumorHets <tumor_het_pulldown>
#   --segments <acnv_segments> --tangentNormalized <coverage_profile> --outputDir <output_dir>
#   --rhoThreshold 0.2 --numIterations 10  --sparkMaster local[*]  

# def _titan_cn_file(cnr_file, work_dir, data):
#   """Convert CNVkit or GATK4 normalized input into TitanCNA ready format.
#     """
# out_file = os.path.join(work_dir, "%s.cn" % (utils.splitext_plus(os.path.basename(cnr_file))[0]))
# support_cols = {"cnvkit": ["chromosome", "start", "end", "log2"],
#   "gatk-cnv": ["CONTIG", "START", "END", "LOG2_COPY_RATIO"]}
# cols = support_cols[cnvkit.bin_approach(data)]
# if not utils.file_uptodate(out_file, cnr_file):
#   with file_transaction(data, out_file) as tx_out_file:
#   iterator = pd.read_table(cnr_file, sep="\t", iterator=True, header=0, comment="@")
# with open(tx_out_file, "w") as handle:
#   for chunk in iterator:
#   chunk = chunk[cols]
# chunk.columns = ["chrom", "start", "end", "logR"]
# if cnvkit.bin_approach(data) == "cnvkit":
#   chunk['start'] += 1
# chunk.to_csv(handle, mode="a", sep="\t", index=False)
# return out_file


# def _titan_het_file(vrn_files, work_dir, paired):
#   assert vrn_files, "Did not find compatible variant calling files for TitanCNA inputs"
# from bcbio.heterogeneity import bubbletree
# class OutWriter:
#   def __init__(self, out_handle):
#   self.writer = csv.writer(out_handle, dialect="excel-tab")
# def write_header(self):
#   self.writer.writerow(["Chr", "Position", "Ref", "RefCount", "Nref", "NrefCount", "NormQuality"])
# def write_row(self, rec, stats):
#   if rec.qual and float(rec.qual) > 0:
#   self.writer.writerow([rec.chrom, rec.pos, rec.ref, stats["tumor"]["depth"] - stats["tumor"]["alt"],
#                         rec.alts[0], stats["tumor"]["alt"], rec.qual])
# return bubbletree.prep_vrn_file(vrn_files[0]["vrn_file"], vrn_files[0]["variantcaller"],
#                                 work_dir, paired, OutWriter)
