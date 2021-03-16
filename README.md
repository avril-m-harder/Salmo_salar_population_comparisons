# Salmo_salar_population_comparisons

#### These scripts include all analyses conducted in [citation for paper comparing LC/SE/LH]

_Lake Champlain_

* gatk.sh: For Lake Champlain BAM files output from HISAT2, uses GATK v.3.8 and v4 to produce GVCF variant files
* genomicsdb.sh: Uses GVCF variants files produced by gatk.sh to call genotypes jointly (output is VCF files)
* variant_filtration.sh: Filters variants in VCF files produced by genomicsdb.sh and concatenates all loci into a single VCF file
* LC_vcf_filtering.R: Filters VCF file produced by variant_filtration.sh based on GATK pass filters, depth, and number of samples genotyped requirements

_Lake Sebago and LaHave River_

* sra_download.sh: Downloads reads for Lake Sebago and LaHave River pools from NCBI SRA (requires sra_access_list.txt)
* trimmomatic.sh: Uses Trimmomatic to clean reads (script specific to Lake Sebago/LaHave River samples, but same settings were used for Lake Champlain individuals) and FastQC to check quality before and after cleaning
* popoolation2.sh: Creates mpileup and sync files
* SE_LH_sync_file_filtering.R: Filters sync file produced by Popoolation2 by depth requirements and checks for variation in read depth across pools

_All 3 populations_
* map_3_pop.R: Creates map of populations in Fig. 1
* hisat2_LC.sh / hisat2_SE_LH.sh: Produce BAM files of aligned reads for Lake Champlain individuals / for Lake Sebago and LaHave River pools
* Fig_03_components_ZBTB18.R: Creates components to assemble gene diagram for Figure 3.
* 3_population_analyses.R: Bulk of analyses for project from identification of final set of >23k SNPs through identification of outlier loci, heterozygosity and <i>F</i><sub>ST</sub> calculations, etc.
