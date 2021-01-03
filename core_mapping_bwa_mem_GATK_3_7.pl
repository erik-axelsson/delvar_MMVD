#!/usr/bin/perl -w 

### Mapping and quality control using bwa, piccard and GATK ###
### Allocate one core per sample ###


## Host specific ##

$system = "milou";

if ($system eq "milou") {
$node_ram = 128;
$cores = 16;
}else {
    if ($system eq "kalkyl") {
	$node_ram = 24;
	$cores = 8;
    }
}

## Process memory ##
$aln_mem = 3.5;
$sampe_mem = 5.4;
$max_cores_aln = int($node_ram/$aln_mem);
if ($max_cores_aln > $cores) {
    $max_cores_aln = $cores;
} 
$max_cores_sampe = int($node_ram/$sampe_mem);
if ($max_cores_sampe > $cores) {
    $max_cores_sampe = $cores;
} 

 
## path to reference ##
$refa = "/proj/b2013119/private/Reference/canFam3.fa";
$refai = "/proj/b2013119/private/Reference/canFam3.fa.fai";

## run dir ##
$run_dir = "/proj/b2013119/nobackup/private/160DG_MEM_mapping/";

## align all samples ##

#$sample_file = "/proj/b2013119/private/sequencing_summary_20131202.txt";
#$sample_file = "/proj/b2013119/private/Westies_unmapped.txt";
#$sample_file = "/proj/b2013119/private/Mapping_analysis/failed_mapping_20140123.txt"
$sample_file = "/proj/b2013119/private/Mapping_analysis/160DG_bwa_mem_failed_170201.txt";

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

#foreach $sample (@samples) {
#    print "$sample\n";
#    $count++;
#    print "$count\n";
#}


## align specific samples ##
#@samples = ('3_33_97');
#@samples = ('B4', 'B6', 'B7', 'B8');


## chromosomes included ##
# automatically load chromosomes
#@chromosomes=();
#open IN, "/proj/b2013119/private/Reference/canFam3.fa.fai" or die;
#while (<IN>) {
#    if (/chr(.+?)\t/) {
#	push @chromosomes, $1;
#    }
#}

#manually determine which chromosomes to work with
@chromosomes=();
#push @chromosomes, 'Un';
push @chromosomes, 'X';
push @chromosomes, 'M';
push @chromosomes, (1..38);


# Path to all fastq files #
@files = </proj/b2013119/private/Sequencing_reports/*report.txt>;
foreach $infile (@files) {
#  print $infile . "\n";
  #$infile = "130522_D00118_0084_BC23GRACXX_report.txt";
  open IN, $infile or die;
  if ($infile =~ m/Sequencing_reports\/(.+)\_report\.txt/) {
    $dir = $1;
  }
 line:while (<IN>) {
#    print;
    if ((/^3\_/) or (/^[A-Z]/)){
      chomp;
      $id = $_;
      $id =~ s/\s//g;
      $id =~ s/\r//g;
#      $id =~ s/_/./g;
#      print "$id\n";
      push @{$ORG{$id}}, $dir; 
    }
  }  
}



### align sample by sample ###
foreach $sample (@samples) {    
    unless (($sample eq '3_33_97') or ($sample eq '3_33_97')) {
    foreach $id (sort keys %ORG) {
	if ($sample eq $id) {
	    @sample_bams = @sampe_list = @sort_bam_list = ();
# Make batch script for aligning sample
	    $sample_dir = "/proj/b2013119/nobackup/private/160DG_MEM_mapping/".$id;
	    system "mkdir $sample_dir";
	    $align_file = $run_dir.$id."/".$id."_align_GATK_3_7_GATK_only.sh";
	    print "$align_file\n";
	    $log_file = $align_file.".log";
	    $error_file=$align_file.".errors";
	    open OUT, ">$align_file" or die;

	    print OUT "\#!/bin/bash                                                             
\#SBATCH -A b2013119                                                             
\#SBATCH -p core                                                                 
\#SBATCH -n 1                                                                   
\#SBATCH -t 2-05:00:00
\#SBATCH -J align
\#SBATCH --mail-type=FAIL
\#SBATCH --mail-user=erik.axelsson\@imbim.uu.se
\#SBATCH -o $log_file

set -e
set -o pipefail

source /home/erik/.bashrc.erik      
\n";

	    print "working on sample: $id\n";
	    print OUT "\#\#\# Align reads \#\#\#\n";
	    
# Mapp each fastq file, pair paired-end reads (bwa mem) and convert to bam and sort and index bam file (Samtools) 
#	    print OUT "\#\#\#Map reads\#\#\#\n";
	    @sam = ();
	    @targets = ();
	    @realigned = ();
	    @left_aligned = ();
	    foreach $dir (@{$ORG{$id}}) {
		print "$dir\t";

		$c_dir = "/proj/b2013119/private/".$dir;
		$c_dir .= "/Sample_".$id;
		@fastaq = glob("$c_dir/*");
		%FILE_PAIRS = ();
		foreach $file (@fastaq) {
		    print "$file\n";
		    if ($file =~ /.+?\_L(\d+)\_R.+?\.fastq\.gz/){
			$lane = $1;
			print "$lane\n";
			push @{$FILE_PAIRS{$lane}}, $file;
		    }
		}
		foreach $lane (sort keys %FILE_PAIRS) {
#		    $sai_pair = "";
		    $fastq_pair = "";
		    foreach $file (@{$FILE_PAIRS{$lane}}) {
			if ($file =~ /.+\/(\d{6})\_.+\/((.+?)\_R\d+\_.+)\.fastq\.gz/){
			    $date = $1;
			    $full_id = $2;
			    $bam_embryo = $date."_".$3;
			    $fastq = $file;
			    $sam = $run_dir.$id."/".$date."_".$full_id.".sam";
			    push @sam, $sam;
			    $read_group = "\'\@RG\\tID:".$date.$lane.$id."\\tPL:illumina\\tLB:".$id."\\tSM:".$id."\'";		
#			    $sam_pair .= $sam." ";
			    $fastq_pair .= $fastq." ";
			    print "printing $fastq\n";
			}
		    }
		    $lane_ba = $run_dir.$id."\/".$bam_embryo;
		    $lane_bam = $run_dir.$id."\/".$bam_embryo.".bam";

		    ####### ARGUMENT PRINTED #################
		    print OUT "#time bwa mem -M -R $read_group $refa $fastq_pair 2> $error_file | samtools view -bS - | samtools sort - $lane_ba &&\n\n";

		    ####### ARGUMENT PRINTED #################
		    print OUT "#time samtools index $lane_bam &&\n\n";
		    push @sample_bams, $lane_bam;
#		    $sampe_command = "time samtools view -bS -o $lane_bam $sam"; 
#		    $sort_bam_command = "time samtools sort $lane_bam";
#		    push @sampe_list, $sampe_command;
#		    push @sort_bam_list, $sort_bam_command;
		}
	    }
	    
   
## convert sam to bam  
#	    print OUT "\#\#\# Convert sam to bam \#\#\#\n";
#	    foreach $sampe (@sampe_list){
#		print OUT "$sampe &&\n\n";
#	    }


## Sort bam
#	    print OUT "\#\#\# Sort bam \#\#\#\n";
#	    foreach $bam (@sort_bam_list){
#		print OUT "$bam &&\n\n";
#	    }


## clean up .sai files
#	    print OUT "\#\#\# Clean-up .sai files \#\#\#\n";
#	    foreach $sai (@sai) {
#		print OUT "rm $sai &&\n\n";
#	    }
	    
	    
# Piccard MergeSamFiles to combine all bam files
	    print OUT "\#\#\# Merge lanes \#\#\#\n";
	    $merged_bam = $run_dir.$id."\/".$id."_merged.sorted.bam";
	    $merge_sam_files = "time java -Xmx6g -d64 -Djava.io.tmpdir=/scratch -jar /pica/sw/apps/bioinfo/picard/1.92/kalkyl/MergeSamFiles.jar ";
	    foreach $bam (@sample_bams){
		$merge_sam_files .= "INPUT=".$bam." "; 
	    }
	    $merge_sam_files .= "OUTPUT=".$merged_bam." SORT_ORDER=coordinate ASSUME_SORTED=false VALIDATION_STRINGENCY=LENIENT";

	    ####### ARGUMENT PRINTED #################
	    print OUT "#$merge_sam_files &&\n\n";



# Index merged bam file
	    ####### ARGUMENT PRINTED #################
	    print OUT "#time samtools index $merged_bam &&\n\n\n\n";



## clean up lane bams
	    print OUT "\#\#\# Clean-up lane bams \#\#\#\n";
	    foreach $lane_bam (@sample_bams){
		####### ARGUMENT PRINTED #################
		print OUT "#rm $lane_bam &&\n\n";
	    }
	    ####### ARGUMENT PRINTED #################
	    print OUT "#wait\n\n\n\n";


# Piccard MarkDuplicates (before splitting into chromosomes to mark interchromosomal duplications)
	    print OUT "\#\#\# Mark duplicates\#\#\# \n";
	    $dup_rem = $run_dir.$id."\/".$id."_merged.sorted.dedup.bam";
	    $metrix = $run_dir.$id."\/".$id."_dedup_metrix.txt";
	    ####### ARGUMENT PRINTED #################
	    print OUT "#time java -Xmx6g -d64 -Djava.io.tmpdir=/scratch -jar /pica/sw/apps/bioinfo/picard/1.92/kalkyl/MarkDuplicates.jar INPUT=$merged_bam OUTPUT=$dup_rem VALIDATION_STRINGENCY=LENIENT METRICS_FILE=$metrix && \n\n";
	    
	    ####### ARGUMENT PRINTED #################
	    print OUT "#wait\n\n\n\n";
	    
# Index merged deduped bam file
	    print OUT "\#\#\# Index deduped files \#\#\#\n";
	    ####### ARGUMENT PRINTED #################
	    print OUT "#time samtools index $dup_rem &&\n\n";
	    ####### ARGUMENT PRINTED #################
	    print OUT "#wait\n\n\n\n";

## Clean up merged bam
	    print OUT "\#\#\# Clean-up merged bam \#\#\#\n";
	    ####### ARGUMENT PRINTED #################
	    print OUT "#rm $merged_bam && \n\n";
	    ####### ARGUMENT PRINTED #################
	    print OUT "#wait\n\n\n\n";

######## Split workflow on chromosomes from hereon ############

	    print OUT "\#\#\# Rerun GATK part of pipeline using newest version (3.7) as old (2.7.2) exited with unknown error. \#\#\#\n";


# GATK RealignerTargetCreator
	    print OUT "\#\#\# Identify targets for indel realigning \#\#\#\n";
	    foreach $chr (@chromosomes) {
	    $ind_real_output = $run_dir.$id."\/".$chr."_forIndelRealigner.intervals";
	    push @targets, $ind_real_output; 
	    ####### ARGUMENT PRINTED #################
	    print OUT "time java -Xmx2g -d64 -jar /pica/sw/apps/bioinfo/GATK/3.7/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $refa -L chr".$chr." -I $dup_rem -o $ind_real_output --known /pica/v8/b2013119/private/Annotations/Indels/Axelsson_2013_indels_canfam3.vcf &&\n\n";
	    }
	    ####### ARGUMENT PRINTED #################
	    print OUT "wait\n\n\n\n";


# GATK IndelRealigner
	    print OUT "\#\#\# Indel realign \#\#\#\n";
	    foreach $chr (@chromosomes) {
	    $ind_real_output = $run_dir.$id."\/".$chr."_forIndelRealigner.intervals";
	    $realigned = $run_dir.$id."\/".$id."_".$chr."_merged.sorted.dedup.indelrealigned.bam";
	    push @realigned, $realigned;
	    ####### ARGUMENT PRINTED #################
	    print OUT "time java -Xmx6g -d64 -jar /pica/sw/apps/bioinfo/GATK/3.7/GenomeAnalysisTK.jar -T IndelRealigner -R $refa -L chr".$chr." -I $dup_rem -targetIntervals $ind_real_output -o $realigned -known /pica/v8/b2013119/private/Annotations/Indels/Axelsson_2013_indels_canfam3.vcf &&\n\n";
	    }
	    ####### ARGUMENT PRINTED #################
	    print OUT "wait\n\n\n\n";



# Clean-up indel targets
	    print OUT "\#\#\# Clean-up Indeltargetintervals \#\#\#\n";
	    foreach $target (@targets) {
		####### ARGUMENT PRINTED #################
		print OUT "rm $target && \n\n";
	    }
	    ####### ARGUMENT PRINTED #################
	    print OUT "wait\n\n\n\n";


## Clean-up marked duplicates bam
#	    print OUT "\#\#\# Clean-up deduped bam \#\#\#\n";
#	    print OUT "rm $dup_rem && \n\n";
#	    print OUT "wait\n\n\n\n";


#GATK leftalign
	    print OUT "\#\#\# Left align \#\#\#\n";
	    foreach $chr (@chromosomes) {
	    $left_aligned = $run_dir.$id."\/".$id."_".$chr."_merged.sorted.dedup.indelrealigned.leftaligned.bam";
	    $realigned = $run_dir.$id."\/".$id."_".$chr."_merged.sorted.dedup.indelrealigned.bam";
	    push @left_aligned, $left_aligned;
	    ####### ARGUMENT PRINTED #################
	    print OUT "time java -Xmx6g -d64 -jar /sw/apps/bioinfo/GATK/3.7/GenomeAnalysisTK.jar -T LeftAlignIndels -R $refa -L chr".$chr." -I $realigned -o $left_aligned &&\n\n";
	    }
	    ####### ARGUMENT PRINTED #################
	    print OUT "wait\n\n\n\n";


# Clean-up indel realigned bams
	    print OUT "\#\#\# Clean-up indel realigned bams \#\#\#\n";
	    foreach $realigned (@realigned) {
		####### ARGUMENT PRINTED #################
		print OUT "rm $realigned && \n\n";
	    }
	    ####### ARGUMENT PRINTED #################
	    print OUT "wait\n\n\n\n";


#GATK BaseRecalibrator on all chromosomes simultaneously for maximum accuracy
	    print OUT "\#\#\# Baserecalibrator \#\#\#\n";
	    $recal_input ="";
	    foreach $chr (@chromosomes) {
	    $left_aligned = $run_dir.$id."\/".$id."_".$chr."_merged.sorted.dedup.indelrealigned.leftaligned.bam";
	    $recal_input .= "-I ".$left_aligned." ";
	    }

	    $recal_table = $run_dir.$id."\/".$id."_recal_data.table";	    
	    ####### ARGUMENT PRINTED #################
	    print OUT "time java -Xmx2g -d64 -jar /pica/sw/apps/bioinfo/GATK/3.7/GenomeAnalysisTK.jar -T BaseRecalibrator -R $refa $recal_input -knownSites /pica/v8/b2013119/private/Annotations/Axelsson_2013_SNPs_canfam3.vcf -knownSites /pica/v8/b2013119/private/Annotations/Indels/Axelsson_2013_indels_canfam3.vcf -o $recal_table &&\n\n";
	    ####### ARGUMENT PRINTED #################
	    print OUT "wait\n\n\n\n";



#GATK PrintReads for each chromosome separately
# -Xmx1g fastest on milou
	    print OUT "\#\#\# Print reads \#\#\#\n";
	    foreach $chr (@chromosomes) {
	    $left_aligned = $run_dir.$id."\/".$id."_".$chr."_merged.sorted.dedup.indelrealigned.leftaligned.bam";
	    $base_cal_bam = $run_dir.$id."\/".$id."_".$chr."_merged.sorted.dedup.indelrealigned.leftaligned.calibrated.bam"; 
	    ####### ARGUMENT PRINTED #################
	    print OUT "time java -Xmx1g -d64 -jar /pica/sw/apps/bioinfo/GATK/3.7/GenomeAnalysisTK.jar -T PrintReads -R $refa -L chr".$chr." -I $left_aligned -BQSR $recal_table -o $base_cal_bam &&\n\n";
	    }
	    ####### ARGUMENT PRINTED #################
	    print OUT "wait\n\n\n\n";



# Clean-up left aligned bams
	    print OUT "\#\#\# Clean-up left aligned bams \#\#\#\n";
	    foreach $left_aligned (@left_aligned) {
		####### ARGUMENT PRINTED #################
		print OUT "rm $left_aligned && \n\n";
	    }
	    ####### ARGUMENT PRINTED #################
	    print OUT "wait\n\n\n\n";

#GATK reduce reads
# skip this step as this only saves 50% file size but uses lots of CPU
#	    $reduced_reads = $run_dir.$id."\/".$id."_merged.sorted.dedup.indelrealigned.leftaligned.calibrated.reduced.bam";
#	    print OUT "java -Xmx2g -d64 -jar /pica/sw/apps/bioinfo/GATK/3.7/GenomeAnalysisTK.jar -T ReduceReads -R $refa -I $base_cal_bam -o $reduced_reads\n\n"; 

	    close OUT;
	    system "chmod +x $align_file";
	    system "sbatch $align_file";

	}
    }
}	
}


