#! /usr/bin/perl -w

### Liftover hg38 coordinates of dog SNPs back to canfam3 to if original canfam3 location is found  ###



#manually determine which chromosomes to work with
@chromosomes=();
#@chromosomes=('6','8','15','17','18','23','24','25','26','30','32','33','34','35','36');
@chromosomes = (1..37);
push @chromosomes, 'X';
push @chromosomes, 'Un';
#push @chromosomes, (38);


foreach $chr (@chromosomes) {
    print "chromosome: $chr\n";
## run dir ##
    $run_dir = "/proj/uppstore2017236/b2013119_nobackup/Liftover/";

## estimated run time ##

    $run_time="02:00:00";

    $runfile = $run_dir.$chr."_liftover.sh";
    print "$runfile\n";
    $log_file = $runfile.".log";
    $human_bed_file=$run_dir."Agrarian_Arctic_SNPs/canfam3_chr".$chr."_all_SNPs_hg38_coordinates.bed";
    $dog_bed_file=$run_dir."Agrarian_Arctic_SNPs/and_back_to_human/canfam3_chr".$chr."_all_SNPs_canfam3_coordinates.bed";
    $dog_unmapped_bed_file=$run_dir."Agrarian_Arctic_SNPs/and_back_to_human/canfam3_chr".$chr."_all_SNPs_canfam3_coordinates.unmapped.bed";
    $chain=$run_dir."hg38ToCanFam3.over.chain";

    open OUT, ">$runfile" or die;
    
    print OUT "\#!/bin/bash                                                             
\#SBATCH -A snic2017-7-392                                                             
\#SBATCH -p core                                                                 
\#SBATCH -n 1                                                                   
\#SBATCH -t $run_time
\#SBATCH -J lifover
\#SBATCH --mail-user=erik.axelsson\@imbim.uu.se
\#SBATCH -o $log_file

source /home/erik/.bashrc.erik      
\n";
    
    print "working on chr: $chr\n";
    print OUT "\#\#\# Liftover hg38 coordinates to canfam3 for dog SNPs\#\#\#\n";
    
### GATK select variants command for each chromosome file ### 

    print OUT "/proj/uppstore2017236/b2013119_nobackup/Liftover/liftOver $human_bed_file $chain $dog_bed_file $dog_unmapped_bed_file\n\n"; 
  
    print OUT "wait\n\n\n\n";
    
    system "chmod +x $runfile";
    system "sbatch $runfile";
}
