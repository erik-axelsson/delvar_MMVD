#! /usr/bin/perl -w

### Grep phylop values for dog snps  ###



#manually determine which chromosomes to work with
@chromosomes=();
#@chromosomes=('6','8','15','17','18','23','24','25','26','30','32','33','34','35','36');
@chromosomes = (1..21);
push @chromosomes, 'X';
push @chromosomes, 'M';
#push @chromosomes, (22);


foreach $chr (@chromosomes) {
    print "chromosome: $chr\n";
## run dir ##
    $run_dir = "/proj/uppstore2017236/b2013119_nobackup/Liftover/Agrarian_Arctic_SNPs/";

## estimated run time ##

    $run_time="00:15:00";

    $runfile = $run_dir.$chr."_grep.sh";
    print "$runfile\n";
    $log_file = $runfile.".log";
    $phylop_file="/proj/uppstore2017236/b2013119/private/Analyses/PHYLOP/100WAY/chr".$chr."_hg38.phyloP100way.bed";
    $pattern_file=$run_dir."hg38_chr".$chr."_canfam3_coordinates_for_all_SNPs.grep.pattern.txt";
    $results=$run_dir."phylop_100way_for_hg38_chr".$chr."_canfam3_SNPs.txt";

    open OUT, ">$runfile" or die;
    
    print OUT "\#!/bin/bash                                                             
\#SBATCH -A snic2017-7-392                                                        
\#SBATCH -p core                                                                 
\#SBATCH -n 1                                                                   
\#SBATCH -t $run_time
\#SBATCH -J grep
\#SBATCH --mail-user=erik.axelsson\@imbim.uu.se
\#SBATCH -o $log_file

source /home/erik/.bashrc.erik      
\n";
    
    print "working on chr: $chr\n";
    print OUT "\#\#\# Grep phylop values for dog SNPs\#\#\#\n";
    
### GATK select variants command for each chromosome file ### 

    print OUT "grep -F -f $pattern_file $phylop_file >$results\n";
    
    print OUT "wait\n\n\n\n";
    
    system "chmod +x $runfile";
    system "sbatch $runfile";
}


