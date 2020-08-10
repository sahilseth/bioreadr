
# https://gatk.broadinstitute.org/hc/en-us/articles/360040098592-CrosscheckFingerprints-Picard-
# This tool calculates the LOD score for identity check between "groups" of data 
# in the INPUT files as defined by the CROSSCHECK_BY argument. 
# A positive value indicates that the data seems to have come from the same 
# individual or, in other words the identity checks out. 
# The scale is logarithmic (base 10), so a LOD of 6 indicates that it is 1,000,000 
# more likely that the data matches the genotypes than not. 
# A negative value indicates that the data do not match. 
# A score that is near zero is inconclusive and can result from low coverage 
# or non-informative genotypes. Each group is assigned a sample identifier 
# (for SAM this is taken from the SM tag in the appropriate readgroup header line, 
# for VCF this is taken from the column label in the file-header. 
# After combining all the data from the same group together, 
# an all-against-all comparison is performed. 
# Results are categorized as one of EXPECTED_MATCH, EXPECTED_MISMATCH, 
# UNEXPECTED_MATCH, UNEXPECTED_MISMATCH, or AMBIGUOUS depending on the 
# LOD score and on whether the sample identifiers of the groups agree: 
# LOD scores that are less than LOD_THRESHOLD are considered mismatches, 
# and those greater than -LOD_THRESHOLD are matches (between is ambiguous). 
# If the sample identifiers are equal, the groups are expected to match. 
# They are expected to mismatch otherwise. 
# The identity check makes use of haplotype blocks defined in the 
# HAPLOTYPE_MAP file to enable it to have higher statistical power for 
# detecting identity or swap by aggregating data from several SNPs in the 
# haplotype block. This enables an identity check of samples with very 
# low coverage (e.g. ~1x mean coverage). ascat
# When provided a VCF, the identity check looks at the PL, GL and GT fields 
# (in that order) and uses the first one that it finds.
gatk_crosscheck <- function(bams,
                            gatk4_exe = opts_flow$get("gatk4_exe"),
                            hapmap_db = "~/ref/human/broad/fingerprint_maps/map_files/hg19_nochr.map"){

    gatk_crosscheck_opts = "--LOD_THRESHOLD=-5 --EXPECT_ALL_GROUPS_TO_MATCH=true --CROSSCHECK_BY SAMPLE "
    # singularity exec --bind /risapps --bind /rsrch3 /rsrch3/home/iacs/sseth/sifs/gatk_4.1.7.0.sif gatk CrosscheckFingerprints \
    # --INPUT=tmp/IPCT-S2006-MOON0051-Cap2023-8-ID01_190503-A00728-0034-BHJ3W3DMXX-1-ATCACG.bwa_recalibed.bam \
    # --INPUT=tmp/IPCT-S4008-MOON0051-Cap2043-4-ID01_190415-A00422-0045-AHK77HDSXX-4-ATCACG.bwa_recalibed.bam \
    # --HAPLOTYPE_MAP=~/ref/human/broad/fingerprint_maps/map_files/hg19_nochr.map --LOD_THRESHOLD=-5  --EXPECT_ALL_GROUPS_TO_MATCH=true --OUTPUT=tmp/sample.crosscheck_metrics
    bams = trk$bam
    bam_list = paste("--INPUT ", bams, sep = "", collapse = " ")
    hapmap_db = tools::file_path_as_absolute(hapmap_db)
    cmd = glue("{gatk4_exe} --java-options '-Xmx32g' CrosscheckFingerprints {bam_list} --OUTPUT gatk_crosscheck.metrics --HAPLOTYPE_MAP={hapmap_db} {gatk_crosscheck_opts}")
    system(cmd)

}