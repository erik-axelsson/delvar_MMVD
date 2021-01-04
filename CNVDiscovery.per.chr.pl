#! /usr/bin/perl -w

##manually determine which chromosomes to work with
@chromosomes=();
#push @chromosomes, 'Un';
push @chromosomes, 'X';
#push @chromosomes, '38';
#push @chromosomes, (2..36);

## set size parameters for discovery ###
## parameters now set as recommended for 1000genomes project (6-8x cov.) ## 
$tw_size = 5000; #set the range of deletion sizes that will be searched for 
$tw_overlap = 2500;
$max_rgl=2500;
$b_precision=200;
$mr_length=2500;


## target all samples ##
$sample_file = "/proj/b2013119/private/sequencing_summary_20131202.txt";
#$sample_file = "/proj/b2013119/private/Westies_unmapped.txt";
#$sample_file = "/proj/b2013119/private/Mapping_analysis/failed_mapping_20140123.txt"
#$sample_file = "/proj/b2013119/private/Mapping_analysis/160DG_bwa_mem_failed_170201.txt";

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
#@samples = ('3_33_97', 'B4');
#@samples = ('B4', 'B6', 'B7', 'B8');

foreach $chr (@chromosomes) {
    
    $run_dir = "/proj/b2013119/nobackup/private/160DG_GenomeSTRIP/CNVDiscovery/".$chr."/";
    system "mkdir $run_dir";
    $bam_list = $run_dir."bam.list";
    system "rm $bam_list";
    open BAM, ">>$bam_list" or die;  
    $log_dir = $run_dir."Logdir/";
    $SVP_file = $run_dir.$chr."_CNVDiscovery.sh";
    system "rm $SVP_file";
    $log_file = $SVP_file.".log";
#$output_vcf = $run_dir.".SVDiscovery.vcf";
    $interval_list=$run_dir."interval.list";
    system "rm $interval_list";
    
    open LIST, ">>$interval_list" or die;

#foreach $chr (@chromosomes) {
    print LIST "chr$chr\n";
    $meta_data_dir = "-md /proj/b2013119/nobackup/private/160DG_GenomeSTRIP/SVPreprocess/$chr/Metadata/ ";
    foreach $sample (@samples) {
	$input_bam = "/proj/b2013119/nobackup/private/160DG_MEM_mapping/".$sample."/".$sample."_".$chr."_merged.sorted.dedup.indelrealigned.leftaligned.calibrated.bam"; 
	print BAM "$input_bam\n";
    }
#}
    close BAM;
    

    open OUT, ">>$SVP_file" or die;
    print OUT "#!/bin/bash
#SBATCH -A b2013119                                                             
#SBATCH -p core                                                                 
#SBATCH -n 1                                                                   
#SBATCH -t 5-01:00:00
#SBATCH -J GENOMEStrip
#SBATCH --mail-user=erik.axelsson\@imbim.uu.se
#SBATCH -o $log_file

#source /home/erik/.bashrc.erik      

module load bwa/0.5.9 
module swap java slurm-drmaa/latest
module load slurm-drmaa
#module load R

outdir=/proj/b2013119/private/GenomeSTRIP/Reference_metadata_bundle 
readLength=100 
reference=canFam3.fa
export SV_DIR=/sw/apps/bioinfo/GenomeSTRiP/2.00.1710/milou/svtoolkit
#export GATK_DIR=/sw/apps/bioinfo
echo \$\{outdir\}
echo \$\{readLength\}
echo \$\{reference\}


#These executables must be on your path.
which java > /dev/null || exit 1
which bwa > /dev/null || exit 1

#The directory containing libbwa.so must be on your LD_LIBRARY_PATH
export LD_LIBRARY_PATH=\$\{SV_DIR\}/bwa:\$\{LD_LIBRARY_PATH\}

#The directory containing libdrmaa.so must be on your LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/sw/apps/build/slurm-drmaa/1.0.7/lib:\$\{LD_LIBRARY_PATH\}


#classpath=\"\$\{SV_DIR\}/lib/SVToolkit.jar:\$\{GATK_DIR\}/GATK/3.7/GenomeAnalysisTK.jar:\$\{GATK_DIR\}/GATK-Queue/3.7/milou/Queue.jar\"

classpath=\"\$\{SV_DIR\}/lib/SVToolkit.jar:\$\{SV_DIR\}/lib/gatk/GenomeAnalysisTK.jar:\$\{SV_DIR\}/lib/gatk/Queue.jar:\$\{SV_DIR\}/lib/gatk/drmaa-common-1.0.jar\"

 
localReference=\$\{outdir\}/`echo \$\{reference\} | awk -F / '{ print \$NF }'`
echo \$\{localReference\}


java -Xmx6g -cp \$\{classpath\} org.broadinstitute.gatk.queue.QCommandLine -S \$\{SV_DIR\}/qscript/discovery/cnv/CNVDiscoveryPipeline.q -S \$\{SV_DIR\}/qscript/SVQScript.q -cp \$\{classpath\} -gatk \$\{SV_DIR\}/lib/gatk/GenomeAnalysisTK.jar -configFile \$\{SV_DIR\}/conf/genstrip_parameters.txt -R /proj/b2013119/private/GenomeSTRIP/Reference_metadata_bundle/canFam3.fa -intervalList $interval_list -I $bam_list $meta_data_dir -jobLogDir $log_dir -ploidyMapFile /proj/b2013119/private/GenomeSTRIP/Reference_metadata_bundle/reference.ploidymap.txt -genomeMaskFile /proj/b2013119/private/GenomeSTRIP/Reference_metadata_bundle/canFam3.mask.fa -genderMapFile /proj/b2013119/private/GenomeSTRIP/Reference_metadata_bundle/160DG_gender_map_file.txt -rmd /proj/b2013119/private/GenomeSTRIP/Reference_metadata_bundle -tilingWindowSize $tw_size -tilingWindowOverlap $tw_overlap -maximumReferenceGapLength $max_rgl -boundaryPrecision $b_precision -minimumRefinedLength $mr_length -l DEBUG -runDirectory $run_dir -jobRunner Drmaa -gatkJobRunner Drmaa -jobNative \"-A b2013119  -p core -n 1 -t 5:00:00 -J GenomeSTRiP\" -run
\n";
	

    close OUT or die;
    system "chmod +x $SVP_file";
    system "sbatch $SVP_file";
}

#unused arguments
#-O $output_vcf
