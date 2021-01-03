# delvar_MMVD
Source code for manuscript: "The genetic consequences of dog breed formation - accumulation of deleterious genetic variation and fixation of mutations associated with myxomatous mitral valve disease in Cavalier King Charles spaniels"

## Fastq sequence data mapping
1. core_mapping_bwa_mem_GATK_3_7.pl <br/>
Perl script executing bwa, picard and GATK mapping procedure for sequencing data.

## Variant calling (SNIVs and INDELs)
1. VarCallScatter.pl <br/>
Perl script executing raw GATK variant calling on chromosome seqments.
2. Catvar.pl <br/>
Perl script concatinating scattered chromosome vcf-files to single chromosome vcf-files. 
