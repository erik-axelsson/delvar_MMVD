#!/bin/bash

#SBATCH -A g2020004                                                            

#SBATCH -p core                                                                

#SBATCH -n 1                                                                   
#SBATCH -t 3:00:00
#SBATCH -J snpEff
#SBATCH --mail-user=erik.axelsson@imbim.uu.se
#SBATCH -o /proj/snic2020-6-127/private/Variant_calling_2020_160DG/snpEff.sh.log

source /home/erik/.bashrc.erik      

### Annotate variants ###

java -Xmx4g -jar /sw/apps/bioinfo/snpEff/3.4/milou/snpEff.jar eff -v CanFam3.1.74 /proj/snic2020-6-127/nobackup/private/private/160DG_MEM_variantcalling/160DG.99.9.recalibrated_variants.vcf > /proj/snic2020-6-127/nobackup/private/private/160DG_MEM_variantcalling/160DG.99.9.recalibrated.variants.eff.vcf
