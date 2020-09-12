# Ssalar_population_comparison

#### These scripts include all analyses conducted in .....

File descriptions:
* hisat2_LC.sh/hisat2_SE_LH.sh: Produce BAM files of aligned reads for Lake Champlain individuals and for Lake Sebago and LaHave River pools
* gatk.sh: For Lake Champlain BAM files output from HISAT2, uses GATK v.3.8 and v4 to produce GVCF variant files
* genomicsdb.sh: Uses GVCF variants files produced by gatk.sh to call genotypes jointly (output is VCF files)
* variant_filtration.sh: Filters variants in VCF files produced by genomicsdb.sh and concatenates all loci into a single VCF file
