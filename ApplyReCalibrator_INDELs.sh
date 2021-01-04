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
java  -Xmx4g -jar  /sw/apps/bioinfo/GATK/3.7/GenomeAnalysisTK.jar -T ApplyRecalibration -R /proj/snic2020-6-127/private/Reference/canFam3.fa -input /proj/snic2020-6-127/nobackup/private/private/160DG_MEM_variantcalling/160DG.whole.genome.concat.recalibrated.snps.raw.indels.vcf -mode INDEL --ts_filter_level 99.0 -recalFile /proj/snic2020-6-127/private/Variant_calling_2020_160DG/160DG.whole.genome.output.INDEL.recal -tranchesFile /proj/snic2020-6-127/private/Variant_calling_2020_160DG/160DG.whole.genome.output.INDEL.tranches -o /proj/snic2020-6-127/nobackup/private/private/160DG_MEM_variantcalling/160DG.99.9.recalibrated_variants.vcf 


#java -jar GenomeAnalysisTK.jar \ 
#    -T ApplyRecalibration \ 
#    -R reference.fa \ 
#    -input recalibrated_snps_raw_indels.vcf \ 
#    -mode INDEL \ 
#    --ts_filter_level 99.0 \ 
#    -recalFile recalibrate_INDEL.recal \ 
#    -tranchesFile recalibrate_INDEL.tranches \ 
#    -o recalibrated_variants.vcf 


