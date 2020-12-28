# Salmo_salar_population_comparisons

#### These scripts include all analyses conducted in [citation for paper comparing LC/SE/LH]

_Lake Champlain_

* gatk.sh: For Lake Champlain BAM files output from HISAT2, uses GATK v.3.8 and v4 to produce GVCF variant files
* genomicsdb.sh: Uses GVCF variants files produced by gatk.sh to call genotypes jointly (output is VCF files)
* variant_filtration.sh: Filters variants in VCF files produced by genomicsdb.sh and concatenates all loci into a single VCF file

_Lake Sebago and LaHave River_

* sra_download.sh: Downloads reads for Lake Sebago and LaHave River pools from NCBI SRA (requires sra_access_list.txt)
* trimmomatic.sh: Uses Trimmomatic to clean reads (script specific to Lake Sebago/LaHave River samples, but same settings were used for Lake Champlain individuals) and FastQC to check quality before and after cleaning
* popoolation2.sh: Creates mpileup and sync files

_All 3 populations_

* hisat2_LC.sh / hisat2_SE_LH.sh: Produce BAM files of aligned reads for Lake Champlain individuals / for Lake Sebago and LaHave River pools
