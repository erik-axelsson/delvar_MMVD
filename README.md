# delvar_MMVD
Source code for manuscript: "The genetic consequences of dog breed formation - accumulation of deleterious genetic variation and fixation of mutations associated with myxomatous mitral valve disease in Cavalier King Charles spaniels"

## Fastq sequence data mapping
1. core_mapping_bwa_mem_GATK_3_7.pl <br/>
Perl script executing bwa, picard and GATK mapping procedure for sequencing data.

## Variant calling (SNVs and INDELs) and initial annotation
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
9. select_variants_batch_maker.pl <br/>
Perl script executing split of whole genome vcf file to per chromosome vcf file.

## Pipeline for Liftover of human genome (hg38) per base pair conservation scores (PhyloP) to canfam3 SNVs
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

## SNV and INDEL additional annotations (calculating FST and writing conservation scores to variant file)
1. vcf_breed_with_annotations_command_line_version.pl <br/>
Perl script calculating and reporting pariwise FST for all variants and breed pairs. Script also reports additional variant annotations such as evolutionary conservation, gene ontology, wolf-, cat- and andean fox alleles.  
2. vcf_breed_with_annotations_per_chr_batch_maker.pl <br/>
Perl script submitting "vcf_breed_with_annotations_command_line_version.pl" to slurm queue.

## Relative amount of deleterious variation  
1. genetic_var_distr_estimator_210112.pl <br/>
Perls script comparing relative abundance of delleterious alleles in breed paris using the R(A/B) statistics.
2. assign_p_value_to_purging_and_neg_sel.r <br/>
R code apending bonferronin corrected p-values to output table from 'genetic_var_distr_estimator_210112.pl'. 

## GenomeStrip CNV detection and annotation
1. CNVDiscovery.per.chr.pl <br/>
Perl script executing GenomeStrip CNV discovery pipeline
2. vcf_breed_CNV_with_annotations.pl <br/>
Perl script for annotating CNVs

## GenomeStrip deletion detection
1. SVPreprocess.per.chr.pl <br/>
Perl script executing GenomeStrip deletions preprocess pipeline
2. SVDiscovery.per.chr.pl <br/>
Perl script executing GenomeStrip deletions discovery pipeline
3. SVGenotyper.per.chr.pl <br/>
Perl script executing GenomeStrip deletions genotyping pipeline
4. vcf_breed_SV_FST_with_annotations.pl <br/>
Perl script for annotating deletions
