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
6. ApplyReCalibrator_SNPs.sh <br/>
Bash script applying variant call recalibration scheme for SNVs.
7. ApplyReCalibrator_INDELs.sh <br/>
Bash script applying variant call recalibration scheme for INDELs.
8. snpEff.sh <br/>
Bash script executing snpEff program for variant effect prediction.

## GenomeStrip CNV calling
1. CNVDiscovery.per.chr.pl <br/>
Perl script executing GenomeStrip CNV discovery pipeline

## GenomeStrip deletions calling
1. SVPreprocess.per.chr.pl <br/>
Perl script executing GenomeStrip deletions preprocess pipeline
2. SVDiscovery.per.chr.pl <br/>
Perl script executing GenomeStrip deletions discovery pipeline
3. SVGenotyper.per.chr.pl <br/>
Perl script executing GenomeStrip deletions genotyping pipeline

## Liftover human genome (hg38) per base pair conservation scores (PhyloP) for canfam3 SNPs
1. wigFix_to_bed.pl <br/>
Perl script converting hg38 phylop values (downloaded from UCSC) from WigFix to BED format.
2. vcf_to_bed.pl <br/>
Perl script converting dog SNPs and INDELs from .vcf file to BED format. 
3. dog_to_human_liftover_batchmaker.pl <br/>
Perl script executing LiftOver of dog SNPs from CanFam3 to hg38.
4. hum_to_dog_liftover_batchmaker.pl <br/>
Perl script executing LiftOver of dog SNPs back from hg38 to canfam3 to test if they map back to original position.
5. sort_by_chromosome.pl <br/>
Perl script sorting LiftOver results based on hg38 chromosomes and print files with grep patterns.
6. grep_batch_maker.pl <br/>
Perl script executing Grep of phylop values for the dog SNPs.
7. dog_SNPs_phylop.pl <br/> 
Perl script printing file with canfam3 SNPs and phylop values.

## FST calculations

## VST calculations

## Relative amount of deleterious variation  
