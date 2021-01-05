#! /usr/bin/perl -w

### Estimate average pairwise FST for variants separated by chromosome   ###



#manually determine which chromosomes to work with
@chromosomes=();
#@chromosomes=('6','8','15','17','18','23','24','25','26','30','32','33','34','35','36');
@chromosomes = (1..37);
push @chromosomes, 'X';
#push @chromosomes, 'Un';
push @chromosomes, (M);
#push @chromosomes, (38);


foreach $chr (@chromosomes) {
    print "chromosome: $chr\n";
## run dir ##
    $run_dir = "/proj/snic2020-6-127/private/Analyses/";

## estimated run time ##

    $run_time="2:00:00";

    $runfile = $run_dir.$chr."_vcf_breed.sh";
    print "$runfile\n";
    $log_file = $runfile.".log";
    open OUT, ">$runfile" or die;
    
    print OUT "\#!/bin/bash                                                             
\#SBATCH -A g2020004                                                             
\#SBATCH -p core                                                                 
\#SBATCH -n 1                                                                   
\#SBATCH -t $run_time
\#SBATCH -J estimate_pairwise_FST
\#SBATCH --mail-user=erik.axelsson\@imbim.uu.se
\#SBATCH -o $log_file

source /home/erik/.bashrc.erik      
\n";
    
    print "working on chr: $chr\n";
    print OUT "\#\#\# Estimate average pairwise fst for variants on selected chromosome \#\#\#\n";
    
### GATK select variants command for each chromosome file ### 
    
    print OUT "./vcf_breed_with_annotations_command_line_version.pl $chr\n";
    
    print OUT "wait\n\n\n\n";
    
    system "chmod +x $runfile";
    system "sbatch $runfile";
}


