#! /usr/bin/perl -w

### Separate vcf file by chromosome   ###

## path to reference ##
$refa = "/proj/snic2020-6-127/private/Reference/canFam3.fa";
$refai = "/proj/snic2020-6-127/private/Reference/canFam3.fa.fai";

#manually determine which chromosomes to work with
@chromosomes=();
#@chromosomes=('6','8','15','17','18','23','24','25','26','30','32','33','34','35','36');
@chromosomes = (1..36);
push @chromosomes, 'X';
push @chromosomes, 'M';
#push @chromosomes, (Un);


foreach $chr (@chromosomes) {
    print "chromosome: $chr\n";
## run dir ##
    $run_dir = "/proj/snic2020-6-127/nobackup/private/private/160DG_MEM_variantcalling/";

## estimated run time ##

    $run_time="3:30:00";



    
    $vcf_file = "/proj/snic2020-6-127/nobackup/private/private/160DG_MEM_variantcalling/160DG.99.9.recalibrated.variants.eff.vcf";
    $chr_vcf_file = $run_dir."ALL_CHR/".$chr."_recalibrated_variants_eff.vcf";
    $output = $run_dir."ALL_CHR/".$chr."_vcf_split.sh";
    print "$vcf_file\n";
    print "$output\n";
    $log_file = $output.".log";
    if ($chr eq 'Un') {
	$chr = "/proj/snic2020-6-127/private/Reference/chrUN_contig_list.intervals";
    }else{
	$chr="chr".$chr;
    }
    open OUT, ">$output" or die;
    
    print OUT "\#!/bin/bash                                                             
\#SBATCH -A g2020004                                                             
\#SBATCH -p core                                                                 
\#SBATCH -n 1                                                                   
\#SBATCH -t $run_time
\#SBATCH -J select_variants
\#SBATCH --mail-user=erik.axelsson\@imbim.uu.se
\#SBATCH -o $log_file

source /home/erik/.bashrc.erik      
\n";
    
    print "working on chr: $chr\n";
    print OUT "\#\#\# Separate vcf file by chromosome \#\#\#\n";
    
### GATK select variants command for each chromosome file ### 
    
    print OUT "java -Xmx6g -jar /sw/apps/bioinfo/GATK/3.7/GenomeAnalysisTK.jar -T SelectVariants -R $refa --variant $vcf_file -L $chr -o $chr_vcf_file &&\n\n";
    
    print OUT "wait\n\n\n\n";
    
    system "chmod +x $output";
    system "sbatch $output";
}


sub numerically { 
    $a<=>$b;
}

#  java
#     -jar GenomeAnalysisTK.jar
#     -T HaplotypeCaller
#     -R reference/human_g1k_v37.fasta
#     -I sample1.bam [-I sample2.bam ...] \
#     -o output.raw.snps.indels.vcf
