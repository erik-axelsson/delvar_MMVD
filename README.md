# delvar_MMVD
Source code for manuscript: "The genetic consequences of dog breed formation - accumulation of deleterious genetic variation and fixation of mutations associated with myxomatous mitral valve disease in Cavalier King Charles spaniels"

## Fastq sequence data mapping
1. core_mapping_bwa_mem_GATK_3_7.pl <br/>
Perl script executing bwa, picard and GATK mapping procedure for sequencing data.

## Variant calling (SNVs and INDELs)
1. VarCallScatter.pl <br/>
Perl script executing raw GATK variant calling on chromosome seqments.
2. Catvar.pl <br/>
Perl script executing concatination of scattered chromosome vcf-files to single chromosome vcf-files. 
3. Catvar_whole_genome.pl <br/>
Perl script executing concatination of per chromosome vcf-files to single whole genome vcf-file.
4. VarCalibrator_SNPs.sh <br/>
Bash script executing variant call recalibration for SNVs.
5. VarCalibrator_INDELs.sh <br/>
Bash script executing variant call recalibration for INDELs.
