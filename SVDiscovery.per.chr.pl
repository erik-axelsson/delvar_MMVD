#! /usr/bin/perl -w

##manually determine which chromosomes to work with
@chromosomes=();
#push @chromosomes, 'Un';
push @chromosomes, 'X';
#push @chromosomes, '37';
#push @chromosomes, (1..36);

## set size parameters for discovery ###
## usually 100bp to 100 kb and 100kb to 10 Mb ## 
$min_size = 100; #set the range of deletion sizes that will be searched for 
$max_size = 100000;
$name = $min_size."_".$max_size;
$window_padding = 3000; #overlap between intervals set above
$window_size = 5000000; #run discovery not on entire chromosome but in windows


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
#    $run_dir = "/proj/b2013119/private/GenomeSTRIP/SVDiscovery/$chr/";
    $run_dir = "/proj/b2013119/nobackup/private/160DG_GenomeSTRIP/SVDiscovery/$chr/";
    mkdir $run_dir;
#    $run_dir = "/proj/b2013119/private/GenomeSTRIP/SVDiscovery/$chr/$name/";
    $run_dir = "/proj/b2013119/nobackup/private/160DG_GenomeSTRIP/SVDiscovery/$chr/$name/";
    mkdir $run_dir;
#    $run_run_dir = "/proj/b2013119/private/GenomeSTRIP/SVDiscovery/$chr/$name/Rundir/";
    $run_run_dir = "/proj/b2013119/nobackup/private/160DG_GenomeSTRIP/SVDiscovery/$chr/$name/Rundir/";
    mkdir $run_run_dir;
#    $pre_processing_dir = "/proj/b2013119/private/GenomeSTRIP/SVPreprocess/$chr/";
    $pre_processing_dir = "/proj/b2013119/nobackup/private/160DG_GenomeSTRIP/SVPreprocess/$chr/";
    $bam_list = $run_dir.$chr."_bam.list";
    system "rm $bam_list";
    open BAM, ">>$bam_list" or die;  
    foreach $sample (@samples) {
	$input_bam = "/proj/b2013119/nobackup/private/160DG_MEM_mapping/".$sample."/".$sample."_".$chr."_merged.sorted.dedup.indelrealigned.leftaligned.calibrated.bam"; 
	print BAM "$input_bam\n";
    }
    close BAM;
    $meta_data_dir = $pre_processing_dir."Metadata/";
    $log_dir = $run_dir."Logdir/";
    $SVP_file = $run_dir.$chr."_SVDiscovery.sh";
    system "rm $SVP_file";
    $log_file = $SVP_file.".log";
    $output_vcf = $run_dir.$chr."_".$name.".SVDiscovery.vcf";
    open OUT, ">>$SVP_file" or die;
    print OUT "#!/bin/bash
#SBATCH -A b2013119                                                             
#SBATCH -p core                                                                 
#SBATCH -n 1                                                                   
#SBATCH -t 1-01:00:00
#SBATCH -J GENOMEStrip
#SBATCH --mail-user=erik.axelsson\@imbim.uu.se
#SBATCH -o $log_file

#source /home/erik/.bashrc.erik      

module load bwa/0.5.9 
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

#classpath=\"\$\{SV_DIR\}/lib/SVToolkit.jar:\$\{GATK_DIR\}/GATK/3.7/GenomeAnalysisTK.jar:\$\{GATK_DIR\}/GATK-Queue/3.7/milou/Queue.jar\"

classpath=\"\$\{SV_DIR\}/lib/SVToolkit.jar:\$\{SV_DIR\}/lib/gatk/GenomeAnalysisTK.jar:\$\{SV_DIR\}/lib/gatk/Queue.jar\"

 
localReference=\$\{outdir\}/`echo \$\{reference\} | awk -F / '{ print \$NF }'`
echo \$\{localReference\}


java -Xmx6g -cp \$\{classpath\} org.broadinstitute.gatk.queue.QCommandLine -S \$\{SV_DIR\}/qscript/SVDiscovery.q -S \$\{SV_DIR\}/qscript/SVQScript.q -cp \$\{classpath\} -gatk \$\{SV_DIR\}/lib/gatk/GenomeAnalysisTK.jar -configFile \$\{SV_DIR\}/conf/genstrip_parameters.txt -R /proj/b2013119/private/GenomeSTRIP/Reference_metadata_bundle/canFam3.fa -L chr$chr -I $bam_list -md $meta_data_dir -bamFilesAreDisjoint true -jobLogDir $log_dir -ploidyMapFile /proj/b2013119/private/GenomeSTRIP/Reference_metadata_bundle/reference.ploidymap.txt -genomeMaskFile /proj/b2013119/private/GenomeSTRIP/Reference_metadata_bundle/canFam3.mask.fa -genderMapFile /proj/b2013119/private/GenomeSTRIP/Reference_metadata_bundle/160DG_gender_map_file.txt -rmd /proj/b2013119/private/GenomeSTRIP/Reference_metadata_bundle -O $output_vcf -minimumSize $min_size -maximumSize $max_size -windowPadding $window_padding -windowSize $window_size -l DEBUG -runDirectory $run_run_dir -run
\n";
	
    close OUT or die;
    system "chmod +x $SVP_file";
    system "sbatch $SVP_file";
}


