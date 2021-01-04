#! /usr/bin/perl -w


################################################################################
### This script annotates nonredundent structural variants called using      ###### the SVGenotyper pipeline within GenomeStrip.                             ###
### Annotations include: Breedwise FST, gene, transcript, effect, phastcons  ###
### GO-terms and allele frequencies.                                         ###
################################################################################
 

$SV_TYPE="DEL"; #"DEL" for SVGenotyper pipeline, "CNV" for CNVDiscovery pipe line.
$min_fst=0.0; #set fst threshold for reporting structural variants 
$promoter=1000; #define promoter region as region within X distance from transcription start 


#include these chromosomes in analyses
@chromosomes=();
#push @chromosomes, 'Un';
#push @chromosomes, 'X';
#push @chromosomes, '38';
push @chromosomes, (1..38);




#annotation files
$gene_file = "/proj/uppstore2017236/b2013119/private/Analyses/ANNOTATIONS/ensGene_canfam3.txt";
$phastcon_file = "/proj/uppstore2017236/b2013119/private/Analyses/ANNOTATIONS/canfam3_phastConsElements100way.bed";
$GO_file = "/proj/uppstore2017236/b2013119/private/Analyses/ANNOTATIONS/Canis_GO_annotation.txt";


### load gene annotations ###
open IN, $gene_file;
while (<IN>) {
  @gene = split, /\t/;
  $chr=$gene[2];
  $chr=~s/chr//;
  $strand=$gene[3];
  $transcr_start=$gene[4];
  $transcr_stop=$gene[5];
  $transcript=$gene[1];
  $translation_start=$gene[6];
  $translation_stop=$gene[7];
  @exon_start=split(/,/, $gene[9]);
  @exon_end=split(/,/, $gene[10]);
  $GENE_ID{$chr}{$transcr_start} = $gene[12];
  $TRANSCRIPT_TYPE{$chr}{$transcr_start} = $gene[13];
  $TRANSCRIPT_ID{$chr}{$transcr_start} = $gene[1];
  $STRAND{$chr}{$transcr_start} = $strand;
  $TRANSCRIPTION_START{$chr}{$transcr_start} = $transcr_stop;
  $TRANSCRIPT{$chr}{$transcr_start}=$transcript;
  $TRANSLATION_START{$chr}{$transcr_start}=$translation_start;
  $TRANSLATION_STOP{$chr}{$transcr_start}=$translation_stop;
  $EXON_START{$chr}{$transcr_start}=[@exon_start];
  $EXON_END{$chr}{$transcr_start}=[@exon_end];
}
close IN;
print "Gene annotations loaded\n"; 

### Load GO annotations ###
open IN, $GO_file;
while (<IN>) {
#  print;
  if (/(.*)\t(.*)\t(.*)\t(.*)\t(.*)\t(.*)\t(.*)\t(.*)\t(.*)\t(.*)\n/) {
    $gene_id=$1;
    $gene_id=~s/\s//g;
    $gene_name=$2;
    $description=$3;
    $GO_term=$8;
    $GO_domain=$9;
    $GENE_NAME{$gene_id}=$gene_name;
    $GENE_DESCRIPTION{$gene_id}=$description;
    if ($GO_domain eq "biological_process") {
      $GO{$gene_name}.=$GO_term.", ";
      $GO{$gene_id}.=$GO_term.", ";
#      print "$gene_id\t$transcript_id\t$GO_term\t$GO_domain\n";
#      print "$GO{$gene_id}\n";
    }
  }
}      
close IN;
print "GO annotations loaded\n";


### load phast con annotations ###
open IN, $phastcon_file;
while (<IN>) {
  @phast = split;
  $chr=$phast[0];
  $chr=~ s/chr//;
  $phast_start = $phast[1];
  $phast_end = $phast[2];
  $lod = $phast[3];
  $lod =~ s/lod\=//;
  $PHAST_START{$chr}{$phast_start}=1;
  $PHAST_END{$chr}{$phast_end}=1;
  $PHAST_LOD{$chr}{$phast_start}=$lod;
  $PHAST_LOD{$chr}{$phast_end}=$lod;
}
close IN;
print "Phastcons loaded\n";





foreach $chr (sort @chromosomes) {
    $infile = "/proj/uppstore2017236/b2013119_nobackup/private/160DG_GenomeSTRIP/SVGenotyper/".$chr."/100_100000/".$chr.".dels.genotypes.annotated.vcf";
    push @infiles, $infile;
}


open OUT, ">>/proj/uppstore2017236/b2013119/private/GenomeSTRIP/Analyses/ANNOTATED_SV/".$SV_TYPE."/FST_".$SV_TYPE."_190424.txt";

print OUT "chr\tstart\tend\tlength\tid\tNA\tNA\tgene_name\tgene\ttranscript\teffect\tdescription\tphastcon\tphastcon_lod\tbeagle_fst\tckcs_fst\tgs_fst\tgr_fst\tlr_fst\tpoodle_fst\trotw_fst\twhwt_fst\ttotal_FST\tGO-annotation\tbeagle_A\tbeagle_B\tbeagle_sz\tckcs_A\tckcs_B\tckcs_sz\tgs_A\tgs_B\tgs_sz\tgr_A\tgr_B\tgr_sz\tlr_A\tlr_B\tlr_sz\tpoodle_A\tpoodle_B\tpoodle_sz\trotw_A\trotw_B\trotw_sz\twhwt_A\twhwt_B\twhwt_sz\n";       


### Breed o sample identifier ###
#$breed = "Beagle";
#$breed = "Poodle";
#$breed = "LabradorRetriever";
#$breed = "GoldenRetriever";
#$breed = "GermanShephard";
#$breed = "WestHighlandWhiteTerrier";
#$breed = "Rottweiler";
#$breed = "CavalierKingCharlesSpaniel";

# List of ids with samples to move
 
$sample_list = "/proj/uppstore2017236/b2013119/private/160_DG_rawdata/Sequencing_reports/sequencing_summary_20131202.txt";
open IN, $sample_list or die;
while (<IN>) {
    if (/(.+?)\t.+?\t(.+?)\n/) {
        $id = $1;
        $breed_id = $2;
#        print "$id\t$breed_id\n";
        $BREED{$id}=$breed_id;
#       if ($breed eq $breed_id) {
	push @samples, $id;
#       }
    }
}
#die;

$i=0;
while ($samples[$i]) {
#    print "$samples[$i]\n";
    $i++;
}





%SAMPLE_ORDER=();
@BREED_ORDER=();
%DOGS=();

############################# Read .vcf files ################################

foreach $infile (@infiles) {
    print "$infile\n";
    open IN, $infile or die;
    $count=0;
  LINE: while (<IN>) {
      $line = $_;
#   last if ($count==1000);
      next LINE if ($line=~/^\##/);
      ### Extract sample order and sample-breed connection ###
      if ($line=~/^\#CH/) {
#        print"$line\n";
	  @header = split('\t', $line);
	  $i=$count=0;
#        print "@header\n";
	  while ($header[$i]) {
	      $sample=$header[$i];
	      $sample =~ s/\n//;
	      $sample =~ s/\s//;
#		  print "$sample\n";
	      if (exists $BREED{$sample}) {
		  $count++;
#		    print "$count\t$sample\n";
		  push @{$SAMPLE_ORDER{$BREED{$sample}}}, $i; #Position of all samples belonging to each breed in vcf file 
		  push @BREED_ORDER, $BREED{$sample}; #Array with breed position in vcf file
		  $DOGS{$i}=$sample; #dog order and sample id
		  $BREED_COUNT{$BREED{$sample}}++;
	      }
	  $i++;
	  }
      }
      
#        print "@BREED_ORDER\n";
      
      
      ### Each SNP line by line ###
      if ($line!~/^#/) { 
	  %FREQ=%A=%B=%ALLELE_FREQS=%SITE_VISE_IND_COUNT=();
	  $alt_count=$sample_tracker=$allele_count=$tot_freq_test=$effect=$high_fst_flag=0;
	  @genotype=@effect=@ids=@rt=();
	  $eff="UNDEF";
	  $text=$pre_text=$prot_id=$the_effect=$effe=$the_rest="";
	  foreach $breed (sort keys %BREED_COUNT) {
#           $REF_ALLELE_FREQ{$breed}=0;
	      $CHROMOSOMES{$breed}=0;
	  }
	  $sv_count++;
	  @header = split('\t', $line);
	  $chr=$header[0];
	  $chr=~ s/chr//;
	  next LINE unless ($header[6] eq 'PASS');
#       print "$line\n";
#	  if ($header[6] eq PASS) { # Passed filtering
#	      $passed++;
#	  }
	  $max = @header-1;
	  @info = split /;/, $header[7];

### only include nonredundent calls
	  next LINE unless ($info[20] eq 'GSDUPLICATES=NA');
#	  print "$info[20]\n";
	  $pos = $header[1];
	  
	  
	  ### Go through all samples ###
	  for ($i=9; $i<=$max; $i++) {
	      $sv_flag=0;
#	    if($BREED_ORDER[$i-9] eq Poodle) {
#		print "$header[1]\t$i\t$BREED_ORDER[$i-9]\t$header[$i]\n";
#	    }

#Example of an individual call at deletion
#GT:CNF:FT:GL:GP:GQ:GSPC 0/0:2.9824:PASS:-0.01,-1.82,-24.54:-0.00,-5.01,-27.73:50:0
#Genotype : fractional copy number : Per-sample genotyoe filter : Genotyope likelihood with no freq prior : Genotyoe likelihood : Genotyoe quality : Number of suportin gread pairs for variant allele 

	      if ($header[$i]=~/(.+?)\:(.+?)\:(.+?)\:(.+?)\:(.+?)\:(.+?)\:(.+?)/) {
		  $GT=$1;
		  $CN=$2;  #fractional copynumber
		  $FI=$3;  #per sample genotype filter
		  $GLN=$4; #Genotype likelihood with no freq prior
		  $GL=$5;  #Genotype likelihood
                  $GQ=$6;  #Genotype quality
		  $DP_ALT=$7; #Number of read paris supporting alternative alleles
		  
		  next if ($FI eq 'LQ'); #individual genotype filer
#		  $GL=$2;
#		  $GQ=$3;
#		  $FT=$4;
#		  $RC=$5;
#		  $DR=$6;
#		  $DV=$7;
#		  $RR=$8;
#		  $RV=$9;
#		  print "$header[2]\t$i\t$BREED_ORDER[$i-9]\t$GT\n";
		  
		  @genotype=split /\//, $GT;
#		  print "@genotype\n";
		  
		  # Count allele freqs and number of structural variants per individual 
		  foreach $genotype (@genotype) {
		      unless ($genotype eq ".") {
			  $ALLELE_FREQS{$BREED_ORDER[$i-9]}{$genotype}++;
			  $CHROMOSOMES{$BREED_ORDER[$i-9]}++; # Number of chromosomes sampled per breed
			  if ($genotype == 1) {	# individual is carrier of structural variant
			      $alt_count++;
			      $sv_flag=1;
			  }
		      }
		  }
		  if ($sv_flag == 1) {#individual is carrier
		      $SITE_VISE_IND_COUNT{$DOGS{$i}}++;
		  }
	      }
	  }

	  ### End each sample loop ###

	  # Skip singleton CNVs (CNVs observed in only one individual across dog pop.
	  next LINE unless ($alt_count>1);

	  

	  ### Store breed specific allele frequencies for SNP ###
	  
	  foreach $breed (sort keys %CHROMOSOMES) {
#	      $text.="$breed\t";
#	      print "$header[2]\t$breed\t$ALLELE_FREQS{$breed}{'0'}\t$CHROMOSOMES{$breed}\n";
	      for ($j=0; $j<2; $j++){
		  if (defined ($ALLELE_FREQS{$breed}{$j}) && ($CHROMOSOMES{$breed})) {         
		      $freq=$ALLELE_FREQS{$breed}{$j}/$CHROMOSOMES{$breed};
		      if ($j==0) {
			  $tot_freq_test+=$freq;
			  
			  ## Count number of segregating deldetions per breed
			  if (($freq<1) && ($freq>0)) {
			      #unless ref state is fixed deletion is segregating
			      $SEGREGATING{$breed}++;
			  }
		      }
		  }else{
		      $freq=0;
		  }
#                   print OUT "GT\t$j\tAF\t$freq\t";
		  $text.="$freq\t";
	      }
#               print OUT "SZ\t$CHROMOSOMES{$breed}\t";
	      $text.="$CHROMOSOMES{$breed}\t";
	  }


	 ###count number of private deletions
	  $fixed=0;
	  foreach $breed (sort keys %CHROMOSOMES) {
	      if (defined ($ALLELE_FREQS{$breed}{'0'}) && ($CHROMOSOMES{$breed})) {         
		  $freq=$ALLELE_FREQS{$breed}{'0'}/$CHROMOSOMES{$breed};
		  if ($freq eq 1) {
		      $fixed++;
		  }else{
		      $private_breed=$breed;
		  }
	      } else {
		  # if breed lacks reference individuals
		  if (defined ($ALLELE_FREQS{$breed}{'1'}) && ($CHROMOSOMES{$breed})) {         
		      $freq=1-$ALLELE_FREQS{$breed}{'1'}/$CHROMOSOMES{$breed};
		      if ($freq eq 1) {
			  $fixed++;
		      }else{
			  $private_breed=$breed;
		      }
		  }
	      }
	  }

	  if ($fixed==7) {
	      $PRIVATE{$private_breed}++;
	  }

          #skip nonvariable sites
	  next LINE if (($tot_freq_test==0) || ($tot_freq_test==8));
	  
	  $end_SV=$info[2];
	  $end_SV=~s/END=//;
	  $sv_start=$header[1];
	  $sv_stop=$end_SV;
	  $length_SV=$end_SV-$header[1];

	  # populate results array
	  $results="chr\tstart\tend\tlength\tid\tNA\tNA\tgene_name\tgene\ttranscript\teffect\tdescription\tphastcon\tphastcon_lod\tbeagle_fst\tckcs_fst\tgs_fst\tgr_fst\tlr_fst\tpoodle_fst\trotw_fst\twhwt_fst\ttotal_FST\tGO-annotation";
#\tbeagle_A\tbeagle_B\tbeagle_sz\tckcs_A\tckcs_B\tckcs_sz\tgs_A\tgs_B\tgs_sz\tgr_A\tgr_B\tgr_sz\tlr_A\tlr_B\tlr_sz\tpoodle_A\tpoodle_B\tpoodle_sz\trotw_A\trotw_B\trotw_sz\twhwt_A\twhwt_B\twhwt_sz\n";
	  @results=split /\t/, $results;
	  $it=0;
	  foreach $item (sort @results) {
	      $rt[$it]="";
	      $it++;
#	      print "$item\n";
	  }
	  $rt[0]="$header[0]";
	  $rt[1]="$header[1]";
	  $rt[2]="$end_SV";
	  $rt[3]="$length_SV";
	  $rt[4]="$header[2]";
#	  $pre_text= "$header[0]\t$header[1]\t$end_SV\t$length_SV\t$header[2]\t";
	  






######################## Annotate with regards to genes ###########################

	  $within_gene=$conserved=$annotation_flag=0;
	  %GENE_TO_GO=();

	  #match against every gene on chromosome
	  foreach $transcr_start (sort numerically keys %{$TRANSCRIPTION_START{$chr}}) {
	      $annotated=0;


#print OUT "chr\tstart\tend\tlength\tid\tNA\tNA\tNA\tgene\ttranscript\teffect\tNA\tphastcon\tphastcon_lod\tbeagle_fst\tckcs_fst\tgs_fst\tgr_fst\tlr_fst\tpoodle_fst\trotw_fst\twhwt_fst\ttotal_FST\tGO-annotation\tbeagle_A\tbeagle_B\tbeagle_sz\tckcs_A\tckcs_B\tckcs_sz\tgs_A\tgs_B\tgs_sz\tgr_A\tgr_B\tgr_sz\tlr_A\tlr_B\tlr_sz\tpoodle_A\tpoodle_B\tpoodle_sz\trotw_A\trotw_B\trotw_sz\twhwt_A\twhwt_B\twhwt_sz\n";       



	      ## promoter region ##
	      if ($STRAND{$chr}{$transcr_start} eq "+") {
		  if ((($sv_start>($transcr_start-$promoter)) && ($sv_start<$transcr_start)) || ( ($sv_stop>($transcr_start-$promoter)) && ($sv_stop<$transcr_start))) {
		      $annotation_flag=3;
		      $within_gene=1;
		      $promoter_count++;
		      $gene_id=$GENE_ID{$chr}{$transcr_start};
		      if ($gene_id=~/(.+?)\.+\d+?/) {
			  $gene_id=$1;
		      }
		      $rt[7]="$GENE_NAME{$gene_id}";
		      $rt[11]="$GENE_DESCRIPTION{$gene_id}";
		      $rt[8]="$GENE_ID{$chr}{$transcr_start}";
		      $rt[9]="$TRANSCRIPT_ID{$chr}{$transcr_start}";
		      $rt[10]="promoter";
		      $GENE_TO_GO{$GENE_ID{$chr}{$transcr_start}}++;
		  }
	      }
	      if ($STRAND{$chr}{$transcr_start} eq "-") {
		  if ((($sv_start>$TRANSCRIPTION_START{$chr}{$transcr_start}) && ($sv_start<($TRANSCRIPTION_START{$chr}{$transcr_start}+$promoter))) || (($sv_stop>$TRANSCRIPTION_START{$chr}{$transcr_start}) && ($sv_stop<($TRANSCRIPTION_START{$chr}{$transcr_start}+$promoter)))) {
#		      ${$RESULTS{$chr}{$sv_start}}[3]= "promoter\t";
		      $annotation_flag=3;
		      $promoter_count++;
		      $within_gene=1;
#		      $pre_text.= "\t\t\t$GENE_ID{$chr}{$transcr_start}\t$TRANSCRIPT_ID{$chr}{$transcr_start}\t";
#		      $pre_text.= "promoter\t";
		      $gene_id=$GENE_ID{$chr}{$transcr_start};
		      if ($gene_id=~/(.+?)\.+\d+?/) {
			  $gene_id=$1;
		      }
		      $rt[7]="$GENE_NAME{$gene_id}";
		      $rt[11]="$GENE_DESCRIPTION{$gene_id}";
    		      $rt[8]="$GENE_ID{$chr}{$transcr_start}";
		      $rt[9]="$TRANSCRIPT_ID{$chr}{$transcr_start}";
		      $rt[10]="promoter";
#		      push (@{$RESULTS_TWO{$chr}{$sv_start}}, "$GENE_ID{$chr}{$transcr_start}\t$TRANSCRIPT_ID{$chr}{$transcr_start}\t");
#		      push (@{$RESULTS_TWO{$chr}{$sv_start}}, "promoter\t");
		      $GENE_TO_GO{$GENE_ID{$chr}{$transcr_start}}++;
		  }
	      }
	      
	      ## sv spans entire gene: sv_start before gene_start and sv_stop after gene_stop ##
	      if ($sv_start<=$transcr_start && $sv_stop>=$TRANSCRIPTION_START{$chr}{$transcr_start}) {
		  $entire_gene++;
		  $coding++;
		  $within_gene=1;
	           #check if transcript type is coding (cmpl) or noncoding (none)
		  if (($annotation_flag<5) && ($TRANSCRIPT_TYPE{$chr}{$transcr_start}=~/cmpl/)) {
#		      ${$RESULTS{$chr}{$sv_start}}[3]="coding\t";
		      $annotation_flag=5;
		      $gene_id=$GENE_ID{$chr}{$transcr_start};
		      if ($gene_id=~/(.+?)\.+\d+?/) {
			  $gene_id=$1;
		      }
		      $rt[7]="$GENE_NAME{$gene_id}";
		      $rt[11]="$GENE_DESCRIPTION{$gene_id}";
		      $rt[8]="$GENE_ID{$chr}{$transcr_start}";
		      $rt[9]="$TRANSCRIPT_ID{$chr}{$transcr_start}";
		      $rt[10]="coding";
		      $annotated=1;
		  }

		  if (($annotation_flag<4) && ($TRANSCRIPT_TYPE{$chr}{$transcr_start}=~/none/)) {
#		      ${$RESULTS{$chr}{$sv_start}}[3]="coding\t";
		      $annotation_flag=4;
		      $gene_id=$GENE_ID{$chr}{$transcr_start};
		      if ($gene_id=~/(.+?)\.+\d+?/) {
			  $gene_id=$1;
		      }
		      $rt[7]="$GENE_NAME{$gene_id}";
		      $rt[11]="$GENE_DESCRIPTION{$gene_id}";
		      $rt[8]="$GENE_ID{$chr}{$transcr_start}";
		      $rt[9]="$TRANSCRIPT_ID{$chr}{$transcr_start}";
		      $rt[10]="noncoding transcript";
		      $annotated=1;
		  }

		  $GENE_TO_GO{$GENE_ID{$chr}{$transcr_start}}++;
	      }
	      
	      
	      ## sv spans part of transcribed gene: sv_start, sv_stop or both inside transcribed gene ##
	      if ((($sv_start>$transcr_start) && ($sv_start<$TRANSCRIPTION_START{$chr}{$transcr_start})) || (($sv_stop>$transcr_start) && ($sv_stop<$TRANSCRIPTION_START{$chr}{$transcr_start}))) {
		  $genic++;
		  $within_gene=1;
#	  push (@{$RESULTS{$chr}{$sv_start}}, "$chr,$sv_start,$sv_stop");
#	  print "inside gene\t";
		  $gene_id=$GENE_ID{$chr}{$transcr_start};
		  if ($gene_id=~/(.+?)\.+\d+?/) {
		      $gene_id=$1;
		  }
		  $rt[7]="$GENE_NAME{$gene_id}";
		  $rt[11]="$GENE_DESCRIPTION{$gene_id}";
		  $rt[8]="$GENE_ID{$chr}{$transcr_start}";
		  $rt[9]="$TRANSCRIPT_ID{$chr}{$transcr_start}";
#		  $rt[10]="transcribed";
#		  $pre_text.= "\t\t\t$GENE_ID{$chr}{$transcr_start}\t$TRANSCRIPT_ID{$chr}{$transcr_start}\t";
#		  $pre_text.= "transcribed\t";		 
#		  push (@{$RESULTS_TWO{$chr}{$sv_start}}, "$GENE_ID{$chr}{$transcr_start}\t$TRANSCRIPT_ID{$chr}{$transcr_start}\t");
#		  push (@{$RESULTS_TWO{$chr}{$sv_start}}, "inside gene\t");
		  $i=0;
 		  $GENE_TO_GO{$GENE_ID{$chr}{$transcr_start}}++;
		  


		  
		  ## sv spans entire coding part of gene. genic:transcribed:entire coding gene ##
		  if (($sv_start<=$TRANSLATION_START{$chr}{$transcr_start}) && ($sv_stop>=$TRANSLATION_STOP{$chr}{$transcr_start})) {
		      $coding++;
#		      print "\tcoding";
		      #check if transcript type is coding (cmpl) or noncoding (none)
		      if (($annotation_flag<5) && ($TRANSCRIPT_TYPE{$chr}{$transcr_start}=~/cmpl/)) { #if previous annotation is of lower weight
			  $annotation_flag=5; #update annotatioin flag with new higher weight
			  $rt[10]="coding"; #update effect with new high weight effect
			  $annotated=1; #flag variant as annotated 
		      }

		      if (($annotation_flag<4) && ($TRANSCRIPT_TYPE{$chr}{$transcr_start}=~/none/)) {
			  $annotation_flag=4;
			  $rt[10]="noncoding transcript";
			  $annotated=1;
		      }

		      $GENE_TO_GO{$GENE_ID{$chr}{$transcr_start}}++;
		  }
		  


		  ### match against single exons ###

		  foreach $exon_start (@{$EXON_START{$chr}{$transcr_start}}){	    	   

#		      if ($TRANSCRIPT_ID{$chr}{$transcr_start}=~ /ENSCAFT00000000999/) {
		    
#			  print STDERR "		      if (($sv_start>$exon_start && $sv_start<${$EXON_END{$chr}{$transcr_start}}[$i]) || ($sv_stop>$exon_start && $sv_stop<${$EXON_END{$chr}{$transcr_start}}[$i]) || ($sv_start<=$exon_start && $sv_stop>=${$EXON_END{$chr}{$transcr_start}}[$i])) {\n";
#		      }


			  ## sv spans entire or part of transcribed exon. gene:transcribed:part of transcribed exon ##
		      if ((($sv_start>$exon_start) && ($sv_start<${$EXON_END{$chr}{$transcr_start}}[$i])) || (($sv_stop>$exon_start) && ($sv_stop<${$EXON_END{$chr}{$transcr_start}}[$i])) || (($sv_start<=$exon_start) && ($sv_stop>=${$EXON_END{$chr}{$transcr_start}}[$i]))) {
			  $transcribed++;
			  
			  
			  ## genic:transcribed:part of coding exon ##
			  if (($sv_start>$TRANSLATION_START{$chr}{$transcr_start} && $sv_start<$TRANSLATION_STOP{$chr}{$transcr_start}) || ($sv_stop>$TRANSLATION_START{$chr}{$transcr_start} && $sv_stop<$TRANSLATION_STOP{$chr}{$transcr_start})) {
			      $coding++;
			      $GENE_TO_GO{$GENE_ID{$chr}{$transcr_start}}++;
			      #check if transcript type is coding (cmpl) or noncoding (none)
			      if (($annotation_flag<5) && ($TRANSCRIPT_TYPE{$chr}{$transcr_start}=~/cmpl/)) {
				  $annotated=1;
				  $annotation_flag=5;
				  $rt[10]="coding";
			      }
			      if (($annotation_flag<4) && ($TRANSCRIPT_TYPE{$chr}{$transcr_start}=~/none/)) {
				  $annotation_flag=4;
			  $rt[10]="noncoding transcript";
				  $annotated=1;
			      }
			      
			  }


		      if ($TRANSCRIPT_ID{$chr}{$transcr_start}=~ /ENSCAFT00000048400/) {
			  print "$TRANSCRIPT_ID{$chr}{$transcr_start}\t$STRAND{$chr}{$transcr_start}\n";
		     

			  if (($sv_start<=$transcr_start) && ($sv_stop=>$TRANSLATION_START{$chr}{$transcr_start})) {
			      print 			  "genic:transcribed:entire 5UTR\n";
print STDERR "			  if ($sv_start<=$transcr_start && $sv_stop=>$TRANSLATION_START{$chr}{$transcr_start}) {
\n";
			  }


			  ## genic:transcribed:entire 5UTR ##
			  if (($sv_start<=$transcr_start) && ($sv_stop=>$TRANSLATION_START{$chr}{$transcr_start})) {
			      $GENE_TO_GO{$GENE_ID{$chr}{$transcr_start}}++;
			      if ($STRAND{$chr}{$transcr_start} eq "+") {
				  if ($annotation_flag<4) {
				      $annotation_flag=4;
				      $annotated=1;
				      $rt[10]="5UTR";
				  }
				  $FIVE_UTR++;
			      }
			      if ($STRAND{$chr}{$transcr_start} eq "-") {
				  if ($annotation_flag<4) {
#				      ${$RESULTS{$chr}{$sv_start}}[3]="3UTR\t";
				      $annotation_flag=4;
				      $annotated=1;
				      $rt[10]="3UTR";
				  }	
				  $THREE_UTR++;
			      }
			  }


			  if (($sv_start>$transcr_start && $sv_start<$TRANSLATION_START{$chr}{$transcr_start}) || ($sv_stop>$transcr_start && $sv_stop<$TRANSLATION_START{$chr}{$transcr_start})) {
print "			  ## genic:transcribed:part of 5UTR ##\n";
			  }
			  
			  ## genic:transcribed:part of 5UTR ##
			  if (($sv_start>$transcr_start && $sv_start<$TRANSLATION_START{$chr}{$transcr_start}) || ($sv_stop>$transcr_start && $sv_stop<$TRANSLATION_START{$chr}{$transcr_start})) {
			      $GENE_TO_GO{$GENE_ID{$chr}{$transcr_start}}++;
			      if ($STRAND{$chr}{$transcr_start} eq "+") {
				  if ($annotation_flag<4) {
#				      ${$RESULTS{$chr}{$sv_start}}[3]="5UTR\t";
				      $annotation_flag=4;
				      $annotated=1;
				      $rt[10]="5UTR";
				  }
				  $FIVE_UTR++;
			      }
			      if ($STRAND{$chr}{$transcr_start} eq "-") {
				  if ($annotation_flag<4) {
#				      ${$RESULTS{$chr}{$sv_start}}[3]="3UTR\t";
				      $annotation_flag=4;
				      $annotated=1;
				      $rt[10]="3UTR";
				  }	
				  $THREE_UTR++;
			      }
			  }


			  if ($sv_start<$TRANSLATION_STOP{$chr}{$transcr_start} && $sv_stop>$TRANSCRIPTION_START{$chr}{$transcr_start}) {
print "			  ## genic:transcribed:entire 3UTR ##\n";
			  }

			  
			  ## genic:transcribed:entire 3UTR ##
			  if (($sv_start<$TRANSLATION_STOP{$chr}{$transcr_start}) && ($sv_stop>$TRANSCRIPTION_START{$chr}{$transcr_start})) {
			      $GENE_TO_GO{$GENE_ID{$chr}{$transcr_start}}++;
			      if ($STRAND{$chr}{$transcr_start} eq "+") {
				  if ($annotation_flag<4) {
#				      ${$RESULTS{$chr}{$sv_start}}[3]="3UTR\t";
				      $annotation_flag=4;
				      $rt[10]="3UTR";
				      $annotated=1;
				  }
				  $THREE_UTR++;
			      }
			      if ($STRAND{$chr}{$transcr_start} eq "-") {
				  if ($annotation_flag<4) {
#				      ${$RESULTS{$chr}{$sv_start}}[3]="5UTR\t";
				      $annotation_flag=4;
				      $rt[10]="5UTR";
				      $annotated=1;
				  }
				  $FIVE_UTR++;
			      }
			  }
			  


			  if (($sv_start>$TRANSLATION_STOP{$chr}{$transcr_start} && $sv_start<$TRANSCRIPTION_START{$chr}{$transcr_start}) || ($sv_stop>$TRANSLATION_STOP{$chr}{$transcr_start} && $sv_stop<$TRANSCRIPTION_START{$chr}{$transcr_start})) {
print "			  ## genic:transcribed:part of 3UTR ##\n";
print STDERR "			  if (($sv_start>$TRANSLATION_STOP{$chr}{$transcr_start} && $sv_start<$TRANSCRIPTION_START{$chr}{$transcr_start}) || ($sv_stop>$TRANSLATION_STOP{$chr}{$transcr_start} && $sv_stop<$TRANSCRIPTION_START{$chr}{$transcr_start})) {
\n";
			  }

			  ## genic:transcribed:part of 3UTR ##
			  if (($sv_start>$TRANSLATION_STOP{$chr}{$transcr_start} && $sv_start<$TRANSCRIPTION_START{$chr}{$transcr_start}) || ($sv_stop>$TRANSLATION_STOP{$chr}{$transcr_start} && $sv_stop<$TRANSCRIPTION_START{$chr}{$transcr_start})) {
			      $GENE_TO_GO{$GENE_ID{$chr}{$transcr_start}}++;
			      if ($STRAND{$chr}{$transcr_start} eq "+") {
				  if ($annotation_flag<4) {
#				      ${$RESULTS{$chr}{$sv_start}}[3]="3UTR\t";
				      $annotation_flag=4;
				      $annotated=1;
				      $rt[10]="3UTR";
				  }
				  $THREE_UTR++;
			      }
			      if ($STRAND{$chr}{$transcr_start} eq "-") {
				  if ($annotation_flag<4) {
#				      ${$RESULTS{$chr}{$sv_start}}[3]="5UTR\t";
				      $annotation_flag=4;
				      $annotated=1;
				      $rt[10]="5UTR";
				  }
				  $FIVE_UTR++;
			      }
			  }
		      }
		      }
		      $i++; #keeps track of exonnumber
		  }
	  
		  #if within transcribed part of gene but still not annotated: intron 
		  if ($annotated == 0) {
		      if ($annotation_flag<2) {
			  $annotation_flag=2;
			  $rt[10]="intron";
			  $annotated=1;
		      }
		  }
		  
		  
	      }
	      
	      
	  }
	  #	print "$chr\t$sv_start\t$sv_stop\t$beagle_fst\t$ckcs_fst\t$gs_fst\t$gr_fst\t$lr_fst\t$poodle_fst\t$rotw_fst\t$whwt_fst\t$beagle_sz\t$ckcs_sz\t$gs_sz\t$gr_sz\t$lr_sz\t$poodle_sz\t$rotw_sz\t$whwt_sz\n";
	  if ($within_gene == 0) {
	      if ($annotation_flag<1) {
		  $annotation_flag=1;
		  $rt[10]="intergenic";
	      }
	      $inter_genic++;
	  }
	  
	  
	  
	  


	  
#################### Check for overlap with PHASTCONS ######################

	  $top_lod=$conserved=0;
	  for ($j=$sv_start; $j<=$sv_stop;$j++) {
	      if ((exists $PHAST_START{$chr}{$j}) || (exists $PHAST_END{$chr}{$j})) {
		  $conserved=1;
		  if ($PHAST_LOD{$chr}{$j}>$top_lod) {
		      $top_lod=$PHAST_LOD{$chr}{$j};
		  }
		  #	  $coordinate="$PHAST_LOD{$chr}{$j}\t".$j.",";
	      }
	  }
	  if ($conserved==1) {
	      $rt[12]="conserved";
	      $rt[13]="$top_lod";
#	      $pre_text.= "conserved\t$top_lod\t";
	      $conserved_sv++;
	  }else{ 
#	      $pre_text.= "\t\t\t";
	  }
     




######################### GO-annotations ################################

      foreach $gene_id (sort keys %GENE_TO_GO) {
#	  $gene_id =$GENE_ID{$chr}{$transcr_start};
	  if ($gene_id=~/(.+?)\..+/) {
	    $gene_id=$1;
	  }
	  if (exists $GO{$gene_id}) {
	      $rt[23]="$GO{$gene_id}";
#	    push (@{$RESULTS_TWO{$chr}{$sv_start}}, "$gene_id\t$GO{$gene_id}");
	  }else{
	      $rt[23]="NO GO-ANNOTATION";
#	    push (@{$RESULTS_TWO{$chr}{$sv_start}}, "$gene_id\tNO GO-ANNOTATION");
	  }
	}





################################# FST ###########################################
	  
	  ### Calculate 1) total FST across all pops and 2) average pairwise FST for each breed 
###
	  $allele_count=2; #Assume that SV always have only 2 alleles
	  $breed_A=$n1=$p1=$Tr=$Tn=$Tn1=$Tn=$Tn1=$Tp=$Tp_avg=$Tfst=$Tc_to_2=$Tc=$Ts=0;
	  $array_pos=14;
#           print "P\n";        
	  ## Start setting total FST variables ##
	  foreach $breed (sort keys %CHROMOSOMES) {
	      $Tr++; #total number of breeds analysed at locus;
	      $Tn+=$CHROMOSOMES{$breed}; # total average sample size across all breeds 
	      if (defined $ALLELE_FREQS{$breed}{0}) {
		  $Tp+=$ALLELE_FREQS{$breed}{0};  #total average sample freq. of reference allele 
	      }   
	  }   
#           $Tp_avg=$Tp_avg/$Tr; #total average allele freq. of reference allele
	  $Tn=$Tn/$Tr; # total average sample size across all breeds 
	  $Tp_avg=$Tp/($Tr*$Tn);      #total average sample freq. of reference allele 
#           print "samples size: $Tr\tAvg sample size: $Tn\tavg. allele freq:$Tp_avg\t";
	  #draw main breed
	  
	  foreach $breed (sort keys %CHROMOSOMES) {
	      
	      ## More  total FST variables ##
	      $Tn1 = $CHROMOSOMES{$breed}; #sample size of breed A
	      if (defined $ALLELE_FREQS{$breed}{0}) {
		  $Tp1 = $ALLELE_FREQS{$breed}{0}/$CHROMOSOMES{$breed}; # allele freq of ref allele in breed
	      }else{
		  $Tp1=0;
	      }
#               print "$Tp1\t";
	      $Ts+= ($Tn1*($Tp1-$Tp_avg)**2)/(($Tr-1)*$Tn); #the sample variance of reference allele freqs over populaitons.
#               print STDERR "$Ts+= ($Tn1*($Tp1-$Tp_avg)**2)/(($Tr-1)*$Tn)\n";
	      $Tdev+= ($Tn1-$Tn)**2; #standard deviation of sample sizes
#               print "sample varianse of allele freqs:$Ts\t";
	      ## End total FST variables
	      
	      
	      
	      $breed_A=$breed; 
	      $AVG_FST{$breed_A}=$AVG_FST_EW{$breed_A}=$ALT_AVG_FST{$breed_A}=$SEC_ALT_AVG_FST{$breed_A}=0;
	      $breed_comparisons=0;
#               print OUT "$breed_A\t";
	      ### The pariwise comparisons ###
	      foreach $breed (sort keys %CHROMOSOMES) {
		  next if ($breed eq $breed_A);
		  $breed_B=$breed;
		  $pavg_weigh=0;
		  $Hs_breed_A=$Hs_breed_B=$Ht_both=$Weir_top=$Weir_base=$pavg_weigh=0;
		  $FST{$breed_A}{$breed_B}=$FST_EW{$breed_A}{$breed_B}=0;
		  #FST estimated for each allele
		  for ($k=0; $k<$allele_count; $k++){
#                       print "$breed_A\t$breed\tallele\t$k\t";
		      $r = 2; #number of samples/populations
		      $n1 = $CHROMOSOMES{$breed_A}; #sample size of breed A
		      if (defined $ALLELE_FREQS{$breed_A}{$k}) {                      
			  $p1 = $ALLELE_FREQS{$breed_A}{$k}/$CHROMOSOMES{$breed_A}; # allele freq of allele $k in breed A
		      }else{ 
			  $p1 = 0;
		      }
#                       print "$p1\t";
		      $n2=$p2=$n_avg=$p_avg=$dev=$c=$c_to_2=0;
		      $n2 = $CHROMOSOMES{$breed_B}; #sample size of breed B
		      if (defined $ALLELE_FREQS{$breed_B}{$k}) {                      
			  $p2 = $ALLELE_FREQS{$breed_B}{$k}/$CHROMOSOMES{$breed_B}; # allele freq of allele $k in breed B
		      }else{ 
			  $p2 = 0;
		      }
#                       print "$p2\t";
		      if (($n1>0) && ($n2>0))  {
			  ### Alternative strategy for more accurate weighting of multiallelic FST estimates
			  $Hs_breed_A+=$p1*$p1;
			  $Hs_breed_B+=$p2*$p2;
			  $Ht_both+=(($p1+$p2)/2)*(($p1+$p2)/2);
			  ### Ends here ###
			  
			  
			  
			  if ($p1!=$p2){
			      $n_avg = ($n1+$n2)/2; #average sample size                      
			      $p_avg = ($p1*$n1+$p2*$n2)/($r*$n_avg); #average  allele freq. in both populations
			      $pavg_weigh+=$p_avg; #sum of average allele freqs for all alleles at site in the two pos compared
			      $s = (($n1*($p1-$p_avg)**2)/(($r-1)*$n_avg))+(($n2*($p2-$p_avg)**2)/(($r-1)*$n_avg)); #the sample variance of reference allele freqs over populaitons.
			      $dev = sqrt((($n1-$n_avg)**2+($n2-$n_avg)**2)/$r);
			      $c = $dev/$n_avg;
			      $c_to_2 = $c**2; #the squared coefficient of variation of sample sizes
			      
			      ##### the formula for calculationg FST #####
			      $FST{$breed_A}{$breed_B} += ($s-(1/(2*$n_avg-1))*(($p_avg*(1-$p_avg))-((($r-1)/$r)*$s)))   /    ((1-(2*$n_avg*$c_to_2)/((2*$n_avg-1)*$r))*($p_avg*(1-$p_avg))+(1+(2*$n_avg*($r-1)*$c_to_2)/((2*$n_avg-1)*$r))*($s/$r));
			      
			      $FST_EW{$breed_A}{$breed_B} += (($s-(1/(2*$n_avg-1))*(($p_avg*(1-$p_avg))-((($r-1)/$r)*$s)))   /    ((1-(2*$n_avg*$c_to_2)/((2*$n_avg-1)*$r))*($p_avg*(1-$p_avg))+(1+(2*$n_avg*($r-1)*$c_to_2)/((2*$n_avg-1)*$r))*($s/$r)))*$p_avg;
			      
			      ### Second Aaternative strategy for more accurate weighting of multiallelic FST estimates
#                               $Weir_top+=($s-(1/(2*$n_avg-1))*(($p_avg*(1-$p_avg))-((($r-1)/$r)*$s)));
#                               $Weir_base+=((1-(2*$n_avg*$c_to_2)/((2*$n_avg-1)*$r))*($p_avg*(1-$p_avg))+(1+(2*$n_avg*($r-1)*$c_to_2)/((2*$n_avg-1)*$r))*($s/$r));
			      
			      
			  }else{
			      $FST{$breed_A}{$breed_B} += 0;
			  }
		      }
#                       print "$FST{$breed_A}{$breed_B}\n";
		  }
		  #average FST across all alleles
		  $FST{$breed_A}{$breed_B}=$FST{$breed_A}{$breed_B}/$allele_count;
		  if ($pavg_weigh>0) {
		      $FST_EW{$breed_A}{$breed_B}=$FST_EW{$breed_A}{$breed_B}/$pavg_weigh;
		  }else{
		     $FST_EW{$breed_A}{$breed_B}=0;
		  }
#                   print "\navg_fst: $FST{$breed_A}{$breed_B}\n";
		  ### Alternative strategy for more accurate weighting of multiallelic FST estimates
		  $Hs_breed_A=1-$Hs_breed_A;
		  $Hs_breed_B=1-$Hs_breed_B;
		  $Ht_both=1-$Ht_both;
		  unless ($Ht_both==0) {
		      $Alt_fst=($Ht_both-($Hs_breed_A+$Hs_breed_B)/2)/$Ht_both;
		      $ALT_AVG_FST{$breed_A}+=$Alt_fst;
		  }else{ 
		      $ALT_AVG_FST{$breed_A}+=0;
		  }
		  ### Ends here ###
		  
		  ### Second alternative strategy for more accurate weighting of multiallelic FST estimates
#                   $Second_alt_fst=$Weir_top/$Weir_base;
#                   $SEC_ALT_AVG_FST{$breed_A}+=$Second_alt_fst;
		  
		  ### Average pairwise FST for breed_A vs all breeds and alleles ###
		  $AVG_FST{$breed_A}+=$FST{$breed_A}{$breed_B};
		  $AVG_FST_EW{$breed_A}+=$FST_EW{$breed_A}{$breed_B};
		  $breed_comparisons++;
#                   print "$breed_B\t$FST{$breed_A}{$breed_B}\t";
	      }
	      $AVG_FST{$breed_A}=$AVG_FST{$breed_A}/$breed_comparisons;               
	      $AVG_FST_EW{$breed_A}=$AVG_FST_EW{$breed_A}/$breed_comparisons;

              ### only report SVs for which at least one breed has extreme fst ###
	      if ($AVG_FST_EW{$breed_A}>$min_fst) {
		  $high_fst_flag=1;
	      }

	      $ALL_AVG_FST_EW{$chr}{$pos}{$breed_A}=$AVG_FST_EW{$breed_A};
	      $ALT_AVG_FST{$breed_A}=$ALT_AVG_FST{$breed_A}/$breed_comparisons;               
#               $SEC_ALT_AVG_FST{$breed_A}=$SEC_ALT_AVG_FST{$breed_A}/$breed_comparisons;       
	      $rt[$array_pos]="$AVG_FST_EW{$breed_A}";
	      $array_pos++;
#	      $pre_text.= "$AVG_FST_EW{$breed_A}\t";
	  }
	  ## Cont. of total FST variables
	  $Tdev= sqrt($Tdev/($Tr-1)); #standard deviation of sample sizes
	  $Tc = $Tdev/$Tn;
	  $Tc_to_2 = $Tc**2; #the squared coefficient of variation of sample sizes
#           print STDERR "$pos\t$Tfst = ($Ts-(1/(2*$Tn-1))*(($Tp_avg*(1-$Tp_avg))-((($Tr-1)/$Tr)*$Ts)))   /    ((1-(2*$Tn*$Tc_to_2)/((2*$Tn-1)*$r))*($Tp_avg*(1-$Tp_avg))+(1+(2*$Tn*($r-1)*$Tc_to_2)/((2*$Tn-1)*$Tr))*($Ts/$r));\n";
	  unless (($Ts==0) || ($Tn<5)) {
	      $Tfst = ($Ts-(1/(2*$Tn-1))*(($Tp_avg*(1-$Tp_avg))-((($Tr-1)/$Tr)*$Ts)))   /    ((1-(2*$Tn*$Tc_to_2)/((2*$Tn-1)*$r))*($Tp_avg*(1-$Tp_avg))+(1+(2*$Tn*($r-1)*$Tc_to_2)/((2*$Tn-1)*$Tr))*($Ts/$r));
	  }else{
	      $Tfst=0;
	  }
	  $rt[$array_pos]="$Tfst";
	  $array_pos++;
#	  $pre_text.="$Tfst";
	  $results_length=@rt;
	  if ($high_fst_flag==1) {
	      for ($i=0; $i<$results_length; $i++) {
		  print OUT "$rt[$i]\t";
	      }
	      print OUT "$text\n";
#	  	  print OUT "$pre_text\t\t$text\n";
	  }

#       print "total FST:$Tfst\n";
#               print "$pos\t$length_1\t$length_2\t$A\t$B\t$chromosomes\tminor allele frequency=$freq\t$allele_freq_vcf\t$tal\n";




###Count segregating structural variants with potential functional effects (defined as either affecting coding seqence or phastcons
	      
	  foreach $breed (sort keys %CHROMOSOMES) {
	      for ($j=0; $j<2; $j++){
		  if (defined ($ALLELE_FREQS{$breed}{$j}) && ($CHROMOSOMES{$breed})) {         
		      $freq=$ALLELE_FREQS{$breed}{$j}/$CHROMOSOMES{$breed};
		      if ($j==0) {
			  
			  ## Count number of segregating deletions per breed
			  if (($freq<1) && ($freq>0)) {
			      #unless ref state is fixed deletion is segregating
			      if (($rt[10] eq 'coding')||($rt[12] eq 'conserved')) {
				  $FUNCTIONAL_SEGREGATING{$breed}++;
			      }
			  }
		      }else{
			  $freq=0;
		      }
		  }
	      }	  
	  }
	  if (($rt[10] eq 'coding')||($rt[12] eq 'conserved')) {
	      foreach $id (sort keys %SITE_VISE_IND_COUNT) {
		  $TOTAL_IND_FUNC{$id}++;
	      }
	  }
	  
      }
  }
}
#    die if ($pos==69688517);
#           }
#       }
    
    
close OUT;
                                
open OUT, ">>/proj/uppstore2017236/b2013119/private/GenomeSTRIP/Analyses/ANNOTATED_SV/".$SV_TYPE."/Segregating_and_private_per_breed_".$SV_TYPE."_190424.txt";
 
foreach $breed (sort keys %SEGREGATING) {
    if ($breed=~/[a-z]/) {
	print OUT "$breed\t$SEGREGATING{$breed}\t$PRIVATE{$breed}\n";
    }
}
close OUT;

open OUT, ">>/proj/uppstore2017236/b2013119/private/GenomeSTRIP/Analyses/ANNOTATED_SV/".$SV_TYPE."/Functional_segregating_per_breed_".$SV_TYPE."_190424.txt";

foreach $breed (sort keys %FUNCTIONAL_SEGREGATING) {
    if ($breed=~/[a-z]/) {
	print OUT "$breed\t$FUNCTIONAL_SEGREGATING{$breed}\n";
    }
}

close OUT;    

open OUT, ">>/proj/uppstore2017236/b2013119/private/GenomeSTRIP/Analyses/ANNOTATED_SV/".$SV_TYPE."/Functional_segregating_per_individual".$SV_TYPE."_190424.txt";

foreach $ind (sort keys %TOTAL_IND_FUNC) {
    print OUT "$ind\t$BREED{$ind}\t$TOTAL_IND_FUNC{$ind}\n";
}


close OUT;    





sub numerically {
    $a<=>$b;
}






### Key to interpreting GenomeStrip vcf files ###


##INFO=<ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around END for imprecise variants">
##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS for imprecise variants">
##INFO=<ID=END,Number=1,Type=Integer,Description="End coordinate of this variant">
##INFO=<ID=GSCOHERENCE,Number=1,Type=Float,Description="Value of coherence statistic">
##INFO=<ID=GSCOHFN,Number=1,Type=Float,Description="Coherence statistic per pair">
##INFO=<ID=GSCOHPVALUE,Number=1,Type=Float,Description="Coherence metric (not a true p-value)">
##INFO=<ID=GSCOORDS,Number=4,Type=Integer,Description="Original cluster coordinates">
##INFO=<ID=GSCORA6,Number=1,Type=Float,Description="Correlation with array intensity from Affy6 arrays">
##INFO=<ID=GSCORI1M,Number=1,Type=Float,Description="Correlation with array intensity from Illumina 1M arrays">
##INFO=<ID=GSCORNG,Number=1,Type=Float,Description="Correlation with array intensity from NimbleGen arrays">
##INFO=<ID=GSDEPTHCALLS,Number=.,Type=String,Description="Samples with discrepant read pairs or low read depth">
##INFO=<ID=GSDEPTHCALLTHRESHOLD,Number=1,Type=Float,Description="Read depth threshold (median read depth of samples with discrepant read pairs)">
##INFO=<ID=GSDEPTHNOBSSAMPLES,Number=1,Type=Integer,Description="Number of samples with discrepant read pairs in depth test">
##INFO=<ID=GSDEPTHNTOTALSAMPLES,Number=1,Type=Integer,Description="Total samples in depth test">
##INFO=<ID=GSDEPTHOBSSAMPLES,Number=.,Type=String,Description="Samples with discrepant read pairs in depth test">
##INFO=<ID=GSDEPTHPVALUE,Number=1,Type=Float,Description="Depth p-value using chi-squared test">
##INFO=<ID=GSDEPTHPVALUECOUNTS,Number=4,Type=Integer,Description="Depth test read counts (carrier inside event, carrier outside event, non-carrier inside, non-carrier outside)">
##INFO=<ID=GSDEPTHRANKSUMPVALUE,Number=1,Type=Float,Description="Depth p-value using rank-sum test">
##INFO=<ID=GSDEPTHRATIO,Number=1,Type=Float,Description="Read depth ratio test">
##INFO=<ID=GSDMAX,Number=1,Type=Integer,Description="Maximum value considered for DOpt">
##INFO=<ID=GSDMIN,Number=1,Type=Integer,Description="Minimum value considered for DOpt">
##INFO=<ID=GSDOPT,Number=1,Type=Integer,Description="Most likely event length">
##INFO=<ID=GSDSPAN,Number=1,Type=Integer,Description="Inner span length of read pair cluster">
##INFO=<ID=GSELENGTH,Number=1,Type=Integer,Description="Effective length">
##INFO=<ID=GSEXPMEAN,Number=1,Type=Float,Description="Expected read depth sample mean">
##INFO=<ID=GSGMMWEIGHTS,Number=.,Type=Float,Description="Genotyping depth model cluster weights">
##INFO=<ID=GSM1,Number=1,Type=Float,Description="Genotyping depth model parameter M1">
##INFO=<ID=GSM2,Number=2,Type=Float,Description="Genotyping depth model parameters M2[0],M2[1]">
##INFO=<ID=GSMEMBNPAIRS,Number=1,Type=Integer,Description="Number of pairs used in membership test">
##INFO=<ID=GSMEMBNSAMPLES,Number=1,Type=Integer,Description="Number of samples used in membership test">
##INFO=<ID=GSMEMBOBSSAMPLES,Number=.,Type=String,Description="Samples participating in membership test">
##INFO=<ID=GSMEMBPVALUE,Number=1,Type=Float,Description="Membership p-value">
##INFO=<ID=GSMEMBSTATISTIC,Number=1,Type=Float,Description="Value of membership statistic">
##INFO=<ID=GSNDEPTHCALLS,Number=1,Type=Integer,Description="Number of samples with discrepant read pairs or low read depth">
##INFO=<ID=GSNHET,Number=1,Type=Integer,Description="Number of heterozygous snp genotype calls inside the event">
##INFO=<ID=GSNHOM,Number=1,Type=Integer,Description="Number of homozygous snp genotype calls inside the event">
##INFO=<ID=GSNNOCALL,Number=1,Type=Integer,Description="Number of snp genotype non-calls inside the event">
##INFO=<ID=GSNPAIRS,Number=1,Type=Integer,Description="Number of discrepant read pairs">
##INFO=<ID=GSNSAMPLES,Number=1,Type=Integer,Description="Number of samples with discrepant read pairs">
##INFO=<ID=GSNSNPS,Number=1,Type=Integer,Description="Number of snps inside the event">
##INFO=<ID=GSOUTLEFT,Number=1,Type=Integer,Description="Number of outlier read pairs on left">
##INFO=<ID=GSOUTLIERS,Number=1,Type=Integer,Description="Number of outlier read pairs">
##INFO=<ID=GSOUTRIGHT,Number=1,Type=Integer,Description="Number of outlier read pairs on right">
##INFO=<ID=GSREADGROUPS,Number=.,Type=String,Description="Read groups contributing discrepant read pairs">
##INFO=<ID=GSREADNAMES,Number=.,Type=String,Description="Discrepant read pair identifiers">
##INFO=<ID=GSRPORIENTATION,Number=1,Type=String,Description="Read pair orientation">
##INFO=<ID=GSSAMPLES,Number=.,Type=String,Description="Samples contributing discrepant read pairs">
##INFO=<ID=GSSNPHET,Number=1,Type=Float,Description="Fraction of het snp genotype calls inside the event">
##INFO=<ID=HOMLEN,Number=.,Type=Integer,Description="Length of base pair identical micro-homology at event breakpoints">
##INFO=<ID=HOMSEQ,Number=.,Type=String,Description="Sequence of base pair identical micro-homology at event breakpoints">
##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">
##INFO=<ID=NOVEL,Number=0,Type=Flag,Description="Indicates a novel structural variation">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">



##ALT=<ID=DEL,Description="Deletion">
##FILTER=<ID=COHERENCE,Description="GSCOHPVALUE == \NA\ || GSCOHPVALUE <= 0.01">
##FILTER=<ID=COVERAGE,Description="GSDEPTHCALLTHRESHOLD == \NA\ || GSDEPTHCALLTHRESHOLD >= 1.0">
##FILTER=<ID=DEPTH,Description="GSDEPTHRATIO == \NA\ || GSDEPTHRATIO > 0.8 || (GSDEPTHRATIO > 0.63 && (GSMEMBPVALUE == \NA\ || GSMEMBPVALUE >= 0.01))">
##FILTER=<ID=DEPTHPVAL,Description="GSDEPTHPVALUE == \NA\ || GSDEPTHPVALUE >= 0.01">
##FILTER=<ID=LQ,Description="Low Quality Genotyping filter">
##FILTER=<ID=PAIRSPERSAMPLE,Description="GSNPAIRS <= 1.1 * GSNSAMPLES">
##FORMAT=<ID=CN,Number=1,Type=Integer,Description="Copy number genotype for imprecise events">
##FORMAT=<ID=CNF,Number=1,Type=Float,Description="Estimate of fractional copy number">
##FORMAT=<ID=CNL,Number=.,Type=Float,Description="Copy number likelihoods with no frequency prior">
##FORMAT=<ID=CNP,Number=.,Type=Float,Description="Copy number likelihoods">
##FORMAT=<ID=CNQ,Number=1,Type=Float,Description="Copy number genotype quality for imprecise events">
##FORMAT=<ID=FT,Number=1,Type=String,Description="Per-sample genotype filter">
##FORMAT=<ID=GL,Number=G,Type=Float,Description="Genotype likelihoods with no frequency prior">
##FORMAT=<ID=GP,Number=G,Type=Float,Description="Genotype likelihoods">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=GSPC,Number=1,Type=Integer,Description="Number of supporting read pairs for structural variant alleles">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">


#Example of an individual call at deletion
#GT:CNF:FT:GL:GP:GQ:GSPC 0/0:2.9824:PASS:-0.01,-1.82,-24.54:-0.00,-5.01,-27.73:50:0
#Genotype : fractional copy number : Per-sample genotyoe filter : Genotyope likelihood with no freq prior : Genotyoe likelihood : Genotyoe quality : Number of suportin gread pairs for variant allele 

