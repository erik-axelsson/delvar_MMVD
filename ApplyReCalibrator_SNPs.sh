#!/bin/bash                                                             
#SBATCH -A g2020004                                                            

#SBATCH -p core                                                                

#SBATCH -n 1                                                                   
#SBATCH -t 10:00:00
#SBATCH -J applyrecal
#SBATCH --mail-user=erik.axelsson@imbim.uu.se
#SBATCH -o /proj/snic2020-6-127/private/Variant_calling_2020_160DG/apply_recalibrator_SNPs.sh.log

source /home/erik/.bashrc.erik      

### Calibrate SNPs ###
java  -Xmx4g -jar  /sw/apps/bioinfo/GATK/3.7/GenomeAnalysisTK.jar -T ApplyRecalibration -R /proj/snic2020-6-127/private/Reference/canFam3.fa -input /proj/snic2020-6-127/nobackup/private/private/160DG_MEM_variantcalling/160DG.whole.genome.concat.raw.snps.indels.vcf -mode SNP --ts_filter_level 99.9 -recalFile /proj/snic2020-6-127/private/Variant_calling_2020_160DG/160DG.whole.genome.output.SNPs.recal -tranchesFile /proj/snic2020-6-127/private/Variant_calling_2020_160DG/160DG.whole.genome.output.SNPs.tranches -o /proj/snic2020-6-127/nobackup/private/private/160DG_MEM_variantcalling/160DG.whole.genome.concat.recalibrated.snps.raw.indels.vcf 



#java -jar GenomeAnalysisTK.jar \ 
#    -T ApplyRecalibration \ 
#    -R reference.fa \ 
#    -input raw_variants.vcf \ 
#    -mode SNP \ 
#    --ts_filter_level 99.0 \ 
#    -recalFile recalibrate_SNP.recal \ 
#    -tranchesFile recalibrate_SNP.tranches \ 
#    -o recalibrated_snps_raw_indels.vcf 

