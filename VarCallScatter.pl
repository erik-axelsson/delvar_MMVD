#! /usr/bin/perl -w

### Variant calling using Haplotype Caller in GATK    ###
### Chromosomes are scattered to parallelise analyses ###

## Scatter chromosomes into windows ##
$window_size = 10000000;


## path to reference ##
$refa = "/proj/snic2020-6-127/private/Reference/canFam3.fa";
$refai = "/proj/snic2020-6-127/private/Reference/canFam3.fa.fai";



## Samples ##

$sample_file = "/proj/snic2020-6-127/nobackup/private/private/160DG_MEM_mapping/sequencing_summary_20131202.txt";


open IN, $sample_file or die;
while (<IN>) {
    if (/(.+?)\s/) {
	$sample = $1;
	$sample =~ s/\r//;
	$sample =~ s/\s//;
	$sample =~ s/\n//;
	$sample =~ s/\./\_/g;
	push @samples, $sample;
  }
}

foreach $sample (@samples) {
    print "$sample\n";
    $count++;
    print "$count\n";
}


## use specific samples ##
#@samples = ('B4');
#@samples = ('B4', 'B6', 'B7', 'B8');


## Genomic intervals for scattering ##
$interval_file = "/proj/snic2020-6-127/private/Reference/canFam3.fa.fai";
open IN, $interval_file;
while (<IN>){ 
    if (/(.+?)\t(.+?)\t/) {
	$chr = $1;
	$size = $2;
	$size =~ s/\s//g;
	$max_size = $size;
	$chr =~ s/chr//;
	for ($i=1; $i<$max_size; $i+=$window_size){
	    print "$chr\t$i\n";
	    push @{$CHR_INT{$chr}}, $i;
	}
	push @{$CHR_INT{$chr}}, $size;	
	$CHR_SIZE{$chr} = $size;
    }
} 
#die;
#foreach $chr (sort numerically keys %CHR_INT) {
#    foreach $int (@{$CHR_INT{$chr}}) {
#	print "$chr\t$int\n";
#    }
#}
#    print "$chr\t$CHR_SIZE{$chr}\n";
#}
#die;

#manually determine which chromosomes to work with
@chromosomes=();
@chromosomes = (1..37);
push @chromosomes, 'X';
push @chromosomes, 'M';
#push @chromosomes, (38);


# Paths to mapped bam files #
foreach $chr (@chromosomes) {
    print "chromosome: $chr\n";
## run dir ##
    $chr_dir = "/proj/snic2020-6-127/nobackup/private/private/160DG_MEM_variantcalling/".$chr."/";
    system "mkdir $chr_dir";
    $run_dir = "/proj/snic2020-6-127/nobackup/private/private/160DG_MEM_variantcalling/".$chr."/".$chr."_scatter/";
    system "mkdir $run_dir";

    $sample_bams = ();
    foreach $sample (@samples) {
	$sample_bam = "/proj/snic2020-6-127/nobackup/private/private/160DG_MEM_mapping/".$sample."/".$sample."_".$chr."_merged.sorted.dedup.indelrealigned.leftaligned.calibrated.bam";
#	print "$sample_bam\n";
	$sample_bams .= "-I ".$sample_bam." "; 
    }
    $count=0;
    $length = @{$CHR_INT{$chr}};
    print "$length\n";
#    die;

    foreach $int (@{$CHR_INT{$chr}}) {
	last if ($count==(@{$CHR_INT{$chr}}-1));
	$start = $int;
	$end = ${$CHR_INT{$chr}}[$count+1];
	if ($count>0) {
	    $start += 1;
	}
	$count++;
	$l_command = "chr".$chr.":".$start."-".$end;
	print "$chr\t$start\t$end\n";

	$varcall_file = $run_dir.$chr."_".$start."_".$end."_varcall.sh";
	$output = $run_dir.$chr."_".$start."_".$end.".raw.snps.indels.vcf";
	print "$varcall_file\n";
	print "$output\n";
	$log_file = $varcall_file.".log";
	open OUT, ">$varcall_file" or die;
	
	print OUT "\#!/bin/bash                                                             
\#SBATCH -A snic2017-7-392                                                             
\#SBATCH -p core                                                                 
\#SBATCH -n 1                                                                   
\#SBATCH -t 1-23:00:00
\#SBATCH -J varcallscatt
\#SBATCH --mail-user=erik.axelsson\@imbim.uu.se
\#SBATCH -o $log_file
\#SBATCH --mail-type=FAIL

source /home/erik/.bashrc.erik      
\n";
    
	print "working on chr: $chr\n";
	print OUT "\#\#\# Call variants using GATK Haplotype Caller \#\#\#\n";
	
### GATK haplotype caller command for each chromosome segment file ### 

	print OUT "time java -Xmx2g -d64 -jar /sw/apps/bioinfo/GATK/3.7/GenomeAnalysisTK.jar -T HaplotypeCaller -R $refa -L $l_command $sample_bams -o $output &&\n\n";

	print OUT "wait\n\n\n\n";

	system "chmod +x $varcall_file";
	system "sbatch $varcall_file";
    }
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
