#!/usr/bin/perl -w

### Concatinate multiple .vcf files into one single for  each chromsome ###
### Remember to change from core to node job depending on size of vcfs. ####

#@chromosomes = (5,6,8,10,12,13,16,17,19,21,22);
#@chromosomes = (1..37);
#push @chromosomes, (M);
#push @chromosomes, (X);
push @chromosomes, (1);


## path to reference ##
$refa = "/proj/snic2020-6-127/private/Reference/canFam3.fa";
$refai = "/proj/snic2020-6-127/private/Reference/canFam3.fa.fai";


foreach $chr (@chromosomes) {
    $catvar_command = (); 
    @new_files = ();
   $run_dir = "/proj/snic2020-6-127/nobackup/private/private/160DG_MEM_variantcalling/".$chr."/".$chr."_scatter/";
    $bash_file = $run_dir.$chr."_catvar.sh";
    $output_file = $run_dir.$chr.".concat.raw.snps.indels.vcf";
    $log_file = $bash_file.".log";
    system "rm $output_file";
   @files = <$run_dir*.vcf>;
    $begin = 1;
    $count = 0;
    ### sort files numerically ###
    while ($count<@files) {
#	print "Hej\n";
	$index= 0;
	foreach $file (@files) {
	    if ($file =~ /^.+?\_(\d+)\_/) {
		$start = $1;
		print "testing $start\n";
		if ($start==$begin) {
		    $file = $files[$index];
		    push @new_files, $file;
		    print "$start matched $begin shift file: $file\n";
		    $begin+=10000000;
		    if ($count==0) {
			$begin+=1;
		    }
		    print "$begin\n";
		    $count++;
#		    die;
		}		
	    }
	    $index++;
	}
    }
   foreach $file (@new_files) {
       $catvar_command .= "-V ".$file." ";
       print "$file\n";
   }
    open OUT, ">$bash_file" or die;
    print OUT "#!/bin/bash
#SBATCH -A snic2017-7-392
#SBATCH -p core                                                                
#SBATCH -n 1                                                                   
#SBATCH -t 0:15:00
#SBATCH -J varcat
#SBATCH --mail-user=erik.axelsson\@imbim.uu.se
#SBATCH -o $log_file

source /home/erik/.bashrc.erik      

### Concat variants ###
java -cp  /sw/apps/bioinfo/GATK/3.7/GenomeAnalysisTK.jar org.broadinstitute.gatk.tools.CatVariants -R $refa $catvar_command -out $output_file -assumeSorted\n";

    system "chmod +x $bash_file";
    system "sbatch $bash_file";

}

sub numerically {
    $a<=>$b;
}
