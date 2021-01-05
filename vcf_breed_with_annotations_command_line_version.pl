#! /usr/bin/perl



######################
### Set parameters ###
######################

$chromosome=$ARGV[0];
$snp_flag=0; #set to 1 if only SNPs are to be analysed. set to "0" if all variants are to be included 
$GT_Q=3; #minimum genotype quality to include call in individual genome 
$DP_C=1; #minimum sequence depth for position in individual 
$effect_sort=0; #0 means no filter. Set to 1 if only SNPs with particular effect are to be analysed.
$select_eff="."; #select type of SNP to analyse ("." for all SNPs).
$min_sz=0; #minimum number of chromosomes sampled per breed to report variable site in output file
$min_fst=0; #set fst threshold for reporting site.
$promoter=1000; #distance in bp from transcription start site used to define promoter region.
 
#$infile = "/proj/b2013119/nobackup/Maja/TangoA_CFA13region_recalibrated_variants.eff.vcf";

#$infile = "/proj/b2013119/private/Variant_calling/10000_Chr1_recalibrated_variants.eff.vcf";

$variant_infile = "/proj/snic2020-6-127/nobackup/private/private/160DG_MEM_variantcalling/ALL_CHR/".$chromosome."_recalibrated_variants_eff.vcf";

print "dog variants: $variant_infile\n";



############################################
### Pre load positions of variable sites ###
############################################

### Read .vcf file ###
open IN, $variant_infile or die;
while (<IN>) {
    $line = $_;
    if ($line!~/^#/) { 
	@header = split('\t', $line);
	$chr=$header[0];
	$chr=~ s/chr//;
	$pos=$header[1];
	$PRE_CHECK_VAR{$chr}{$pos}=1;
    }
}

print "Variable sites pre-loadinbg finished\n"; 

########################
### LOAD ANNOTATIONS ###
########################

$gene_file = "/proj/snic2020-6-127/private/Analyses/ANNOTATIONS/ensGene_canfam3.txt";
$hundred_way_phylop_file = "/proj/snic2020-6-127/nobackup/private/Liftover/PHYLOP/canfam3_SNPs_with_100WAY_phylop.txt";
$fortysix_way_phylop_file = "/proj/snic2020-6-127/nobackup/private/Liftover/PHYLOP/canfam3_SNPs_with_46WAY_phylop.txt";
$GO_file = "/proj/snic2020-6-127/private/Analyses/ANNOTATIONS/Canis_GO_annotation.txt";
$phastcon_file="/proj/snic2020-6-127/private/Analyses/ANNOTATIONS/canfam3_phastConsElements100way.bed";

$cat_file="/proj/snic2020-6-127/nobackup/private/CatToDog/canFam3.felCat8.net.axt";


### Load aligned cat genome ###
open IN, $cat_file or die;
$line=1;
reader:while (<IN>) {
    chomp;
    unless (/^#/) {
	if (/\d+/){
            @info=split, /\t/;
            $chr=$info[1];
            $chr=~s/chr//;
	    next reader unless ($chr==$chromosome);
            $start=$info[2];
            $end=$info[3];
            $len=$end-$start;
#            print "$len\n";
            $line=2;
	}
        if ((/^[A-Z,a-z]+/) && ($line==3)){
            $cat=$_;
            $cat_length=length($cat);
 #           print "cat:$cat\n";
            $line=1;
            $j=0; #keeps track of dog base                                                                 
            for ($i=0; $i<$cat_length; $i++) { #keeps track of cat base                                    
		unless ((substr $dog, $i ,1) eq "-"){#only store cat base if aligned to dog (skip gaps)    
#               print "$i\n";
		    if (exists $PRE_CHECK_VAR{$chr}{$start+$j}) {        
			$CAT{$chr}{$start+$j}=substr $cat, $i ,1;
		    }
                    $j++;
		}
            }
	}
        if ((/^[A-Z,a-z]+/) && ($line==2)){
            $dog=$_;
            $line++;
#            print "dog:$dog\n";
	}
    }
#    foreach $chr (sort keys %CAT) {
#	foreach $pos (sort keys %{$CAT{$chr}}) {
#           print "$chr\t$pos\t$CAT{$chr}{$pos}\n";
#	}
#    }
}


print "Cat genome loaded\n";


### load wolf genotype ###
$wolf_vcf="/proj/snic2020-6-127/nobackup/private/private/Public_canid_variant_calling/chr".$chromosome."/".$chromosome.".concat.raw.snps.indels.vcf";
print "wolf variants: $wolf_vcf\n";
open IN, $wolf_vcf;
while (<IN>) {
#    print;
    if (/^chr/) {
#      print;
	@_=split, /\t/;
	$canfam3_chr=$canfam3_pos="";
	$canfam3_chr=$_[0];
	$canfam3_chr=~s/chr//;
	if ($canfam3_chr==$chromosome) {
	    $canfam3_pos=$_[1];
#	    print "$canfam3_chr\t$canfam3_pos\t$_[7]\n";
	    @info=split /;/, $_[7];
#	    print "@info\n";
	    if ($info[1]=~/AF\=(.+)/) {
		$wolf_alt_freq=$1;
#		print "$wolf_alt_freq\n";
#	print "$canfam3_chr\t$canfam3_pos\t$wolf_alt_freq\n";
		$WOLF_ALT_FREQ{$canfam3_chr}{$canfam3_pos}=$wolf_alt_freq;
	    }  
	}
    }
}
print "wolf data loaded\n";
close IN;


### load andean fox genotype ###
$af_vcf="/proj/snic2020-6-127/nobackup/private/private/Public_canid_variant_calling/AF_genotype_given_alleles/chr".$chromosome."/".$chromosome.".concat.raw.snps.indels.vcf";
print "wolf variants: $wolf_vcf\n";
open IN, $af_vcf;
while (<IN>) {
#    print;
    if (/^chr/) {
#      print;
	@_=split, /\t/;
	$canfam3_chr=$canfam3_pos="";
	$canfam3_chr=$_[0];
	$canfam3_chr=~s/chr//;
	if ($canfam3_chr==$chromosome) {
	    $canfam3_pos=$_[1];
	    $andean_ref=$_[3];
	    $andean_alt=$_[4];
#	    print "$canfam3_chr\t$canfam3_pos\t$_[7]\n";
	    @info=split /;/, $_[7];
#	    print "@info\n";
	    if ($info[1]=~/AF\=(.+)/) {
		$af_alt_freq=$1;
#		print "$af_alt_freq\n";
#	print "$canfam3_chr\t$canfam3_pos\t$wolf_alt_freq\n";
		if ($af_alt_freq==1) {
		    $ANDEAN{$canfam3_chr}{$canfam3_pos}=$andean_alt;
		}
		if ($af_alt_freq==0) {
		    $ANDEAN{$canfam3_chr}{$canfam3_pos}=$andean_ref;
		}
		if ($af_alt_freq==0.5) {
		    $ANDEAN{$canfam3_chr}{$canfam3_pos}=$andean_ref."/".$andean_alt;
		}

	    }else{
		$ANDEAN{$canfam3_chr}{$canfam3_pos}="NA";
	    }
	}
    }
}
print "andean fox data loaded\n";
close IN;


#foreach $chr(sort numerically keys %WOLF_ALT_FREQ) {
#    foreach $pos (sort numerically keys %{$WOLF_ALT_FREQ{$chr}}) {
#	print "$chr\t$pos\t$WOLF_ALT_FREQ{$chr}{$pos}\n";
#    }
#}


### load 100 way phylop scores ###
open IN, $hundred_way_phylop_file;
while (<IN>) {
    @phylop=split,/\t/;
    $chr=$phylop[0];
    $chr=~s/\s//;
    if ($chr==$chromosome) {
	$pos=$phylop[1];
	$pos=~s/\s//;
	$phylo=$phylop[4];
	$phylo=~s/\s//;
	$PHYLOP{$chr}{$pos}=$phylo;
	$PHYLOP_HUM_CHR{$chr}{$pos}=$phylop[2];
	$PHYLOP_HUM_POS{$chr}{$pos}=$phylop[3];
    }
}
print "100 way phylop scores loaded\n";
close IN;


### load 46 way phylop scores ###
open IN, $fortysix_way_phylop_file;
while (<IN>) {
  @phylop=split,/\t/;
  $chr=$phylop[0];
  $chr=~s/\s//;
  if ($chr==$chromosome) {
      $pos=$phylop[1];
      $pos=~s/\s//;
      $phylo=$phylop[4];
      $phylo=~s/\s//;
      $PHYLOP20{$chr}{$pos}=$phylo;
#  $PHYLOP20_HUM_CHR{$chr}{$pos}=$phylop[2];
#  $PHYLOP20_HUM_POS{$chr}{$pos}=$phylop[3];
  }
}
print "46 way phylop scores loaded\n";
close IN;




### load gene annotations ###
open IN, $gene_file;
while (<IN>) {
  @gene = split, /\t/;
  $chr=$gene[2];
  $chr=~s/chr//;
  if ($chr==$chromosome) {
      $strand=$gene[3];
      $transcr_start=$gene[4];
      $transcr_stop=$gene[5];
      $transcript=$gene[1];
      $translation_start=$gene[6];
      $translation_stop=$gene[7];
      @exon_start=split(/,/, $gene[9]);
      @exon_end=split(/,/, $gene[10]);
      $GENE_ID{$chr}{$transcr_start} = $gene[12];
      $TRANSCRIPT_ID{$chr}{$transcr_start} = $gene[1];
      $STRAND{$chr}{$transcr_start} = $strand;
      $TRANSCRIPTION_START{$chr}{$transcr_start} = $transcr_stop;
      $TRANSCRIPT{$chr}{$transcr_start}=$transcript;
      $TRANSLATION_START{$chr}{$transcr_start}=$translation_start;
      $TRANSLATION_STOP{$chr}{$transcr_start}=$translation_stop;
      $EXON_START{$chr}{$transcr_start}=[@exon_start];
      $EXON_END{$chr}{$transcr_start}=[@exon_end];
  }
}
close IN;
print "Gene annotations loaded\n";


### Load GO annotations ###
open IN, $GO_file;
while (<IN>) {
#  print;
  if (/(.*)\t(.*)\t(.*)\t(.*)\t(.*)\t(.*)\t(.*)\t(.*)\t(.*)\t(.*)\n/) {
#      print "1:$1\t2:$2\t3:$3\t4:$4\t5:$5\t6:$6\t7:$7\t8:$8\t9:$9\t10:$10\n";
    $gene_id=$1;
    $gene_id=~s/\s//g;
    $gene_name=$2;
    $description=$3;
    $GO_term=$8;
    $GO_domain=$9;
#    print "$GO_term\t";
    $GO_term=~s/\s/_/g;
#    print "$GO_term\t";
#    $GO_domaine=~s/\s//g;
#    print "$GO_domain\n";
    $GENE_NAME{$gene_id}=$gene_name;
    $GENE_DESCRIPTION{$gene_id}=$description;
    if ($GO_domain eq "biological_process") {
      $GO{$gene_name}.=$GO_term.",";
      $GO{$gene_id}.=$GO_term.",";
#      print "$gene_id\t$transcript_id\t$GO_term\t$GO_domain\t$GO{$gene_name}\n";
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
  if ($chr==$chromosome) {
      $phast_start = $phast[1];
      $phast_end = $phast[2];
      $lod = $phast[3];
      $lod =~ s/lod\=//;
      $PHAST_START{$chr}{$phast_start}=1;
      $PHAST_END{$chr}{$phast_end}=1;
      $PHAST_LOD{$chr}{$phast_start}=$lod;
      $PHAST_LOD{$chr}{$phast_end}=$lod;
  }
}
close IN;
print "Phastcon scores loaded\n";





################################################
#### START lOADING VARIANTS FROM VCF FILE ######
################################################


### Breed o sample identifier ###
#$breed = "Beagle"
#$breed = "Poodle";
#$breed = "LabradorRetriever";
#$breed = "GoldenRetriever";
#$breed = "GermanShephard";
#$breed = "WestHighlandWhiteTerrier";
#$breed = "Rottweiler";
#$breed = "CavalierKingCharlesSpaniel";

# List of ids with samples to move
 
$sample_list = "/proj/snic2020-6-127/nobackup/private/private/160DG_MEM_mapping/sequencing_summary_20131202.txt";
open IN, $sample_list or die;
while (<IN>) {
    if (/(.+?)\t.+?\t(.+?)\n/) {
	$id = $1;
	$breed_id = $2;
	print "$id\t$breed_id\n";
	$BREED{$id}=$breed_id;
#	if ($breed eq $breed_id) {
	    push @samples, $id;
#	}
    }
}
#die;

$i=0;
while ($samples[$i]) {
    print "$samples[$i]\n";
    $i++;
}

%SAMPLE_ORDER=();
@BREED_ORDER=();
open OUT, ">>/proj/snic2020-6-127/nobackup/private/private/160DG_MEM_variantcalling/FST_chr".$chromosome."_all_variants_sz".$min_sz."_annotated.txt";
open SOUT, ">>/proj/snic2020-6-127/nobackup/private/private/160DG_MEM_variantcalling/all_pariwise_FST_chr".$chromosome."_all_variants_sz".$min_sz."_annotated.txt";

print OUT "chr\tpos\thum_chr\thum_pos\tref\talt\tcat\tandean_fox\twolf_allele_freq\tgene\tsnpeff_effect\tall_effects\t100way\t46way\tbeagle_fst\tckcs_fst\tgs_fst\tgr_fst\tlr_fst\tpoodle_fst\trotw_fst\twhwt_fst\ttotal_FST\tGO-annotation\tbeagle_A\tbeagle_B\tbeagle_sz\tckcs_A\tckcs_B\tckcs_sz\tgs_A\tgs_B\tgs_sz\tgr_A\tgr_B\tgr_sz\tlr_A\tlr_B\tlr_sz\tpoodle_A\tpoodle_B\tpoodle_sz\trotw_A\trotw_B\trotw_sz\twhwt_A\twhwt_B\twhwt_sz\n"; 	


### Read .vcf file ###
open IN, $variant_infile or die;
$count=$all_pair_FST_flag=0;
LINE: while (<IN>) {
    $line = $_;
#   last if ($count==1000);

    ### Extract header information ###
    if ($line=~/^\#CH/) {
	print"$line\n";
	@header = split('\t', $line);
	$i=$count=0;
	print "@header\n";
	while ($header[$i]) {
#	    print "$i\n";
	    $sample = $header[$i];
	    $sample =~ s/\n//;
	    if (exists $BREED{$sample}) {
		$count++;
		print "$count\t$sample\n";
		push @{$SAMPLE_ORDER{$BREED{$sample}}}, $i; #Position of all samples belonging to each breed in vcf file 
		push @BREED_ORDER, $BREED{$sample}; #Array with breed position in vcf file
		$BREED_COUNT{$BREED{$sample}}++;
	    }
	    $i++;
    }

	print "@BREED_ORDER\n";
	#reads header line
#	die;
    }

    ### Each SNP line by line ###
    if ($line!~/^#/) { 
	%FREQ=%A=%B=%ALLELE_FREQS=%FST=%FST_EW=();
	$sample_tracker=$allele_count=$tot_freq_test=$effect=$high=$pos=$chr=0;
	@genotype=@effect=@ids=();
	$eff="UNDEF";
	$text=$prot_id=$the_effect=$effe=$the_rest=$all_fst_comparisons="";
	foreach $breed (sort keys %BREED_COUNT) {
#	    $REF_ALLELE_FREQ{$breed}=0;
	    $CHROMOSOMES{$breed}=0;
	}
	$count++;
	@header = split('\t', $line);
	$chr=$header[0];
	$chr=~ s/chr//;
	$pos=$header[1];
#	print "$line\n";
	if ($header[6] eq PASS) { # Passed filtering
	    $length_1=length($header[3]);
	    $length_2=length($header[4]);
	    @alternative_alleles=split /,/, $header[4];
	    $allele_count=@alternative_alleles+1;
#	    print "@alternative_alleles\t$allele_count\n";
	    ### If $snp_flag is activated only biallelic SNPs are analysed, otherwise all variants including INDELS are analysed ###
	    if ($snp_flag==1) {
		next LINE unless (($length_1==1) && ($length_2 == 1));
	    }
	    $max = @header-1;
#	print "$max\n";
#		print "$header[7]\n";
	    @info = split /;/, $header[7];
#		    $allele_freq_vcf=$info[1];
	    if ($info[1]=~/AF\=(\d\.\d+e*\-*\d*)/) {
		$allele_freq_vcf=$1; #Allele frequency reported by GATK
	    }
#	    print "$header[1]\n";
	    ### Optional variant effect filter ###
#	    if (($info[16]=~/EFF\=(.+)/) ||($info[17]=~/EFF\=(.+)/) ||($info[12]=~/EFF\=(.+)/)) { #SNP effect
	    if ($header[7]=~/EFF\=(.+)\;*/) { #SNP effect
		$effect=$1;
	    }
#	    print "$effect\n";
	    if ($effect_sort==1) {
		next LINE unless ($effect=~/$select_eff/);
	    }
	    @effect=split /,/, $effect;

### PASTED IN FROM SNP_ANNOTATE.PL ###
#extract most relevant snpeff effect of snp
	    $effect_count=@effect;
      #      print "$effect_count\n";
	    foreach $effe (@effect) {
#	print "$effe\n";
#	print "$effect_count\n";
		#	die if ($effe =~/STOP/);
	    }
	    if ($effect_count==1) {
		foreach $effe (@effect) { 
		    @e_split = split/\|/, $effe;
		    $gene_name=$e_split[5];
#	  print "$effe\n";
		    if ($effe=~ /(.+?)\((.+?)\)/) {
			$the_effect=$1;
			$the_rest=$2;
		    }
		}
	    }
	    
	    
	    
	    if ($effect_count>1) {
		if ($effect[0]=~ /(.+?)\((.+?)\)/) {
		    $the_effect=$1;
		    $the_rest=$2;
		    @e_split = split/\|/, $effect[0];
		    $gene_name=$e_split[5];
		}
		foreach $effe (@effect) { 
		    if ($effe=~ /(.+?)\((HIGH.+?)\)/) {
#	    print "HIGH!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n$effe\n";
			#	    die;
			$the_effect=$1;
			$the_rest=$2;
			$high=1;
			@e_split = split/\|/, $effe;
			$gene_name=$e_split[5];
			
		    }
		}
		unless ($high==1) {
		    foreach $effe (@effect) { 
			if ($effe=~ /(.+?)\((MODERATE.+?)\)/) {
			    $the_effect=$1;
			    $the_rest=$2;
			    $high=1;
			    @e_split = split/\|/, $effe;
			    $gene_name=$e_split[5];
			    
			}
		    }
		}
	    }
    
	    
#print "$chr\t$pos\t$ref\t$alt\t$gene_name\t$the_effect\n";
#die unless ($gene eq "NA");
#die if ($high_eff =~ /STOP/);
	    unless ($gene_name=~ /./) {
		$gene_name="NA";
	    }
#wolf alternative allele freq
	    if (exists $WOLF_ALT_FREQ{$chr}{$pos}) {
		$wolf=$WOLF_ALT_FREQ{$chr}{$pos};
	    }else{
		$wolf="0";
	    }

#andean fox base
	    if (exists $ANDEAN{$chr}{$pos}) {
		$andean=$ANDEAN{$chr}{$pos};
	    }else{
		$andean="NA";
	    }

#cat base
	    if (exists $CAT{$chr}{$pos}) {
		$cat=$CAT{$chr}{$pos};
	    }else{
		$cat="NA";
	    }
	    # Get human hg38 position of dog variant
	    if (exists $PHYLOP_HUM_CHR{$chr}{$pos}) {
		$hum_chr=$PHYLOP_HUM_CHR{$chr}{$pos};
		$hum_pos=$PHYLOP_HUM_POS{$chr}{$pos};

	    } else {
		$hum_chr="NA";
		$hum_pos="NA";
### PRINT DIRECTLY TO OUTPUT FILE !!!
	    }
#      push (@{$RESULTS{$chr}{$sv_start}}, $sv_start);
#      push (@{$RESULTS{$chr}{$sv_start}}, $sv_stop);
#      push (@{$RESULTS{$chr}{$sv_start}}, "0");
#      print "$chr\t$pos\t$the_effect\t";

	    
	    
	    
	    ## GO-annotations ##
#      foreach $gene_id (sort keys %GENE_TO_GO) {
##	  $gene_id =$GENE_ID{$chr}{$transcr_start};
#	  if ($gene_id=~/(.+?)\..+/) {
#	    $gene_id=$1;
#	  }
	    if (exists $GO{$gene_name}) {
		$go_annotation=$GO{$gene_name};
#	    push (@{$RESULTS_TWO{$chr}{$pos}}, "$gene_id\t$GO{$gene_id}");
	    }else{
		$go_annotation="NA";
#	    push (@{$RESULTS_TWO{$chr}{$pos}}, "$gene_id\tNO GO-ANNOTATION");
	    }
#	}
	    
	    
	    
	    ### Check for overlap with 100 way Phylop ###
#print STDERR "      if (exists $PHYLOP{$chr}{$pos}) {\n";
#      print "$chr\t$pos\n";
	    if (exists $PHYLOP{$chr}{$pos}) {
#	print "PHYLOP scoure\n";
		$phylo=$PHYLOP{$chr}{$pos};
#	push (@{$RESULTS{$chr}{$pos}}, "$phylo\t");
		$phylop_sv++;
	    }else{
		$phylo="NA";
#	push (@{$RESULTS{$chr}{$pos}}, "NA\t");
#	print "NO svore\n";
	    }
	    
	    
	    ### Check for overlap with 20 way Phylop ###
#print STDERR "      if (exists $PHYLOP{$chr}{$pos}) {\n";
#      print "$chr\t$pos\n";
	    if (exists $PHYLOP20{$chr}{$pos}) {
#	print "PHYLOP scoure\n";
		$twenty_phylo=$PHYLOP20{$chr}{$pos};
#	push (@{$RESULTS{$chr}{$pos}}, "$twenty_phylo");
		$twenty_phylop_sv++;
	    }else{ 
		$twenty_phylo="NA";
#	push (@{$RESULTS{$chr}{$pos}}, "NA");
#	print "NO svore\n";
	    }
#     }
	    

	    
	    
	    
	    
	    ### sample by sample in each vcf file line ###
	    
	    for ($i=9; $i<=$max; $i++) {
#		    if($BREED_ORDER[$i-9] eq Poodle) {
#	    print "$pos\t$i\t$BREED_ORDER[$i-9]\t$header[$i]\n";
#		    }
		
		### Each sample at a time ###
		if ($header[$i]=~/(.+?)\:(.+?)\:(\d+?)\:(.+?)\:(.+)/) {
		    $GT=$1; #genotype
		    $AD=$2; #allele depth
		    $DP=$3; #total depth
		    $GQ=$4; #genotype quality
		    $PL=$5;
		    $DEPTH{$DP}++; ## Individual depth at SNP positions. Use as proxy for seq depth in individual genomes.
#			print "Hej\t$pos\tsample: $i\t$DP\n";
		    $pos = $header[1];
		    @genotype=split /\//, $GT;
#			$if ($DP=~/\d+/) { # Genotypes with missing depth data are not considered
		    if ($GQ>=$GT_Q) { ### set genotyping quality threshold
			if ($DP>=$DP_C){ ## Coverage threshold for genotypes. If 1 genotypes cannot be called.
			  
			    
			    ## Genotype for each individual ##
			    if ($DP>1) {
				foreach $genotype (@genotype) {
#				    print STDERR "$BREED_ORDER[$i-9]\t$genotype\t$ALLELE_FREQS{$BREED_ORDER[$i-9]}{$genotype}++;\n";
#				    print STDERR "$BREED_ORDER[$i-9]\tchromosomes\t$CHROMOSOMES{$BREED_ORDER[$i-9]}++;\n";
				    $ALLELE_FREQS{$BREED_ORDER[$i-9]}{$genotype}++;
				    $CHROMOSOMES{$BREED_ORDER[$i-9]}++; # Number of chromosomes sampled per breed
				}
			    }else{
#				    print STDERR "$BREED_ORDER[$i-9]\t$genotype[0]\t$ALLELE_FREQS{$BREED_ORDER[$i-9]}{$genotype[0]}++;\n";
#				    print STDERR "$BREED_ORDER[$i-9]\tchromosomes\t$CHROMOSOMES{$BREED_ORDER[$i-9]}++;\n";
				$ALLELE_FREQS{$BREED_ORDER[$i-9]}{$genotype[0]}++;
				$CHROMOSOMES{$BREED_ORDER[$i-9]}++; # Number of chromosomes sampled per breed
			    }
#				    if ($GT =~ /0\/0/) {
#					$A{$BREED_ORDER[$i-9]}+=2;
#				    }
#				    if ($GT =~ /0\/1/) {
#					$A{$BREED_ORDER[$i-9]}+=1;
#					$B{$BREED_ORDER[$i-9]}+=1;
#				    }
#				    if ($GT =~ /1\/0/) {
#					$A{$BREED_ORDER[$i-9]}+=1;
#					$B{$BREED_ORDER[$i-9]}+=1;
#				    }
#				    if ($GT =~ /1\/1/) {
#					$B{$BREED_ORDER[$i-9]}+=2;
#				    }
#		print "$GT\t$AD\t$DP\t$GQ\t$PL\n";
			}
#			    }
		    }
		} else {
		    $DEPTH{0}++;
		}
		
		### End each sample loop ###
		
		
		$sample_tracker = $i;
	    }
	    ### End all samples loop ###




	    ### Store breed specific allele frequencies for SNP ###
	    
	    foreach $breed (sort keys %CHROMOSOMES) {
		
### Folded allele frequency spectrum at biallelic sites in breed ### 
		if ($allele_count==2) {           
		    foreach $genotype (sort keys %{$ALLELE_FREQS{$breed}}) {               
			$freq=$ALLELE_FREQS{$breed}{$genotype}/$CHROMOSOMES{$breed};
			last; #store allele freq of first allele
		    }
		    if ($freq>0.5) {
			$freq=1-$freq;
		    }
		    $FREQ{$breed}=$freq;
		    
		    ### Frequencies are binned as we are at most sampling 40 chromosomes per breed ###
		    if ($FREQ{$breed}<=0.05) {
			$ALLELE_FREQS_BIN{$breed}{0.05}++;
		    }
		    if ($FREQ{$breed}>0.05 && $FREQ{$breed}<=0.1) {
			$ALLELE_FREQS_BIN{$breed}{0.1}++;
		    }
		    if ($FREQ{$breed}>0.1 && $FREQ{$breed}<=0.15) {
			$ALLELE_FREQS_BIN{$breed}{0.15}++;
		    }
		    if ($FREQ{$breed}>0.15 && $FREQ{$breed}<=0.2) {
			$ALLELE_FREQS_BIN{$breed}{0.2}++;
		    }
		    if ($FREQ{$breed}>0.2 && $FREQ{$breed}<=0.25) {
			$ALLELE_FREQS_BIN{$breed}{0.25}++;
		    }
		    if ($FREQ{$breed}>0.25 && $FREQ{$breed}<=0.3) {
			$ALLELE_FREQS_BIN{$breed}{0.3}++;
		    }
		    if ($FREQ{$breed}>0.3 && $FREQ{$breed}<=0.35) {
			$ALLELE_FREQS_BIN{$breed}{0.35}++;
		    }
		    if ($FREQ{$breed}>0.35 && $FREQ{$breed}<=0.4) {
			$ALLELE_FREQS_BIN{$breed}{0.4}++;
		    }
		    if ($FREQ{$breed}>0.4 && $FREQ{$breed}<=0.45) {
			$ALLELE_FREQS_BIN{$breed}{0.45}++;
		    }
		    if ($FREQ{$breed}>0.45 && $FREQ{$breed}<=0.5) {
			$ALLELE_FREQS_BIN{$breed}{0.5}++;
		    }
		}
		
		
		### Frequencies of all alleles in breed ###
#		print OUT "$breed\t";
#		$text.="$breed\t";
		    for ($j=0; $j<$allele_count; $j++){
			if (defined ($ALLELE_FREQS{$breed}{$j}) && ($CHROMOSOMES{$breed})) {	     
			    $freq=$ALLELE_FREQS{$breed}{$j}/$CHROMOSOMES{$breed};
			    if ($j==0) {
				$tot_freq_test+=$freq;
			    }
			}else{
			    $freq=0;
			    if ($j==0) {
				$tot_freq_test+=$freq;
			    }
			}
#		    print OUT "GT\t$j\tAF\t$freq\t";
			$text.="$freq\t";
		}
#		print OUT "SZ\t$CHROMOSOMES{$breed}\t";
		$text.="$CHROMOSOMES{$breed}\t";
		next LINE if ($CHROMOSOMES{$breed}<$min_sz);
	    }

### skip sites that are not variable in data set ###
	    next LINE if ($tot_freq_test==0);
	    

#### Print site info and annotations ####  
	    print OUT "$chr\t$pos\t$hum_chr\t$hum_pos\t$header[3]\t$header[4]\t$cat\t$andean\t$wolf\t$gene_name\t$the_effect\t$effect\t$phylo\t$twenty_phylo\t";

	    
	    ### Calculate 1) total FST across all pops and 2) average pairwise FST for each breed ###
#	    print OUT "Avg FST\t";
	    $breed_A=$n1=$p1=$Tr=$Tn=$Tn1=$Tn=$Tn1=$Tp=$Tp_avg=$Tfst=$Tc_to_2=$Tc=$Ts=0;
#	    print "P\n";	
	    ## Start setting total FST variables ##
	    foreach $breed (sort keys %CHROMOSOMES) {
		$Tr++; #total number of breeds analysed at locus;
		$Tn+=$CHROMOSOMES{$breed}; # total average sample size across all breeds 
		if (defined $ALLELE_FREQS{$breed}{0}) {
		$Tp+=$ALLELE_FREQS{$breed}{0};	#total average sample freq. of reference allele	
		}   
	    }	
#	    $Tp_avg=$Tp_avg/$Tr; #total average allele freq. of reference allele
	    $Tn=$Tn/$Tr; # total average sample size across all breeds 
	    $Tp_avg=$Tp/($Tr*$Tn);	#total average sample freq. of reference allele	
#	    print "samples size: $Tr\tAvg sample size: $Tn\tavg. allele freq:$Tp_avg\t";
		#draw main breed

	    foreach $breed (sort keys %CHROMOSOMES) {

		## More  total FST variables ##
		$Tn1 = $CHROMOSOMES{$breed}; #sample size of breed A
		if (defined $ALLELE_FREQS{$breed}{0}) {
		    $Tp1 = $ALLELE_FREQS{$breed}{0}/$CHROMOSOMES{$breed}; # allele freq of ref allele in breed
		}else{
		    $Tp1=0;
		}
#		print "$Tp1\t";
		$Ts+= ($Tn1*($Tp1-$Tp_avg)**2)/(($Tr-1)*$Tn); #the sample variance of reference allele freqs over populaitons.
#		print STDERR "$Ts+= ($Tn1*($Tp1-$Tp_avg)**2)/(($Tr-1)*$Tn)\n";
		$Tdev+= ($Tn1-$Tn)**2; #standard deviation of sample sizes
#		print "sample varianse of allele freqs:$Ts\t";
		## End total FST variables


		
		$breed_A=$breed; 
		$AVG_FST{$breed_A}=$AVG_FST_EW{$breed_A}=$ALT_AVG_FST{$breed_A}=$SEC_ALT_AVG_FST{$breed_A}=0;
		$breed_comparisons=0;
#		print OUT "$breed_A\t";

		### The pariwise comparisons ###
		foreach $breed (sort keys %CHROMOSOMES) {
		    next if ($breed eq $breed_A);
		    $breed_B=$breed;
#		    next if ((exists $FST_EW{$breed_B}) && (defined $FST_EW{$breed_B}{$breed_A}));
#		    print STDERR "		    next if ((exists $FST_EW{$breed_B}) && (defined $FST_EW{$breed_B}{$breed_A}));\n";
		    $pavg_weigh=0;
		    $Hs_breed_A=$Hs_breed_B=$Ht_both=$Weir_top=$Weir_base=$pavg_weigh=0;
		    $FST{$breed_A}{$breed_B}=$FST_EW{$breed_A}{$breed_B}=0;
		    #FST estimated for each allele
		    for ($k=0; $k<$allele_count; $k++){
#			print "$breed_A\t$breed\tallele\t$k\t";
			$r = 2; #number of samples/populations
			$n1 = $CHROMOSOMES{$breed_A}; #sample size of breed A
			if (defined $ALLELE_FREQS{$breed_A}{$k}) {			
			    $p1 = $ALLELE_FREQS{$breed_A}{$k}/$CHROMOSOMES{$breed_A}; # allele freq of allele $k in breed A
			}else{ 
			    $p1 = 0;
			}
#			print "$p1\t";
			$n2=$p2=$n_avg=$p_avg=$dev=$c=$c_to_2=0;
			$n2 = $CHROMOSOMES{$breed_B}; #sample size of breed B
			if (defined $ALLELE_FREQS{$breed_B}{$k}) {			
			    $p2 = $ALLELE_FREQS{$breed_B}{$k}/$CHROMOSOMES{$breed_B}; # allele freq of allele $k in breed B
			}else{ 
			    $p2 = 0;
			}
#			print "$p2\t";
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
#				$Weir_top+=($s-(1/(2*$n_avg-1))*(($p_avg*(1-$p_avg))-((($r-1)/$r)*$s)));
#				$Weir_base+=((1-(2*$n_avg*$c_to_2)/((2*$n_avg-1)*$r))*($p_avg*(1-$p_avg))+(1+(2*$n_avg*($r-1)*$c_to_2)/((2*$n_avg-1)*$r))*($s/$r));
				
				
			    }else{
				$FST{$breed_A}{$breed_B} += 0;
			    }
			}
#			print "$FST{$breed_A}{$breed_B}\n";
		    }
		    #average FST across all alleles
		    $FST{$breed_A}{$breed_B}=$FST{$breed_A}{$breed_B}/$allele_count;
		    if ($pavg_weigh>0) {
		    $FST_EW{$breed_A}{$breed_B}=$FST_EW{$breed_A}{$breed_B}/$pavg_weigh;
		    }else{
			$FST_EW{$breed_A}{$breed_B}=0;
		    }
		    unless ((exists $FST_EW{$breed_B}) && (defined $FST_EW{$breed_B}{$breed_A})){
			$all_fst_comparisons.="$FST_EW{$breed_A}{$breed_B}\t";
			if ($all_pair_FST_flag==0) {
			    $all_pair_FST_header.="$breed_A.$breed_B\t";
			}

		    }
		    
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
#		    $Second_alt_fst=$Weir_top/$Weir_base;
#		    $SEC_ALT_AVG_FST{$breed_A}+=$Second_alt_fst;
		    
		    ### Average pairwise FST for breed_A vs all breeds and alleles ###
		    $AVG_FST{$breed_A}+=$FST{$breed_A}{$breed_B};
		    $AVG_FST_EW{$breed_A}+=$FST_EW{$breed_A}{$breed_B};
		    $breed_comparisons++;
#		    print "$breed_B\t$FST{$breed_A}{$breed_B}\t";
		}
		$AVG_FST{$breed_A}=$AVG_FST{$breed_A}/$breed_comparisons;		
		$AVG_FST_EW{$breed_A}=$AVG_FST_EW{$breed_A}/$breed_comparisons;	
		$ALL_AVG_FST_EW{$chr}{$pos}{$breed_A}=$AVG_FST_EW{$breed_A};
		$ALT_AVG_FST{$breed_A}=$ALT_AVG_FST{$breed_A}/$breed_comparisons;		
#		$SEC_ALT_AVG_FST{$breed_A}=$SEC_ALT_AVG_FST{$breed_A}/$breed_comparisons;		
		print OUT "$AVG_FST_EW{$breed_A}\t";
	    }
	    ## Cont. of total FST variables
	    $Tdev= sqrt($Tdev/($Tr-1)); #standard deviation of sample sizes
	    $Tc = $Tdev/$Tn;
	    $Tc_to_2 = $Tc**2; #the squared coefficient of variation of sample sizes
#	    print STDERR "$pos\t$Tfst = ($Ts-(1/(2*$Tn-1))*(($Tp_avg*(1-$Tp_avg))-((($Tr-1)/$Tr)*$Ts)))   /    ((1-(2*$Tn*$Tc_to_2)/((2*$Tn-1)*$r))*($Tp_avg*(1-$Tp_avg))+(1+(2*$Tn*($r-1)*$Tc_to_2)/((2*$Tn-1)*$Tr))*($Ts/$r));\n";
	    unless (($Ts==0) || ($Tn<5)) {
		$Tfst = ($Ts-(1/(2*$Tn-1))*(($Tp_avg*(1-$Tp_avg))-((($Tr-1)/$Tr)*$Ts)))   /    ((1-(2*$Tn*$Tc_to_2)/((2*$Tn-1)*$r))*($Tp_avg*(1-$Tp_avg))+(1+(2*$Tn*($r-1)*$Tc_to_2)/((2*$Tn-1)*$Tr))*($Ts/$r));
	    }else{
		$Tfst=0;
	    }
	}
    

	if ($header[6] eq PASS) { # Passed filtering
	    
	    print OUT "$Tfst\t$go_annotation\t$text\n";
	    if 	($all_pair_FST_flag==0) {
		print SOUT "chr\tpos\t$all_pair_FST_header\n";
		$all_pair_FST_flag=1;
	    }
	    print SOUT "$chr\t$pos\t$all_fst_comparisons\n";
	}
	
	
	
	
#	print "total FST:$Tfst\n";
#		print "$pos\t$length_1\t$length_2\t$A\t$B\t$chromosomes\tminor allele frequency=$freq\t$allele_freq_vcf\t$tal\n";
    }
#    die if ($pos==69688517);
#	    }
#	}
}

close OUT;
close SOUT;

print "Done with site wise FST estimates\n";

$win=10000;
$slide=5000;
foreach $chr (%ALL_AVG_FST_EW) {
  snp:    foreach $pos (sort numerically keys %{$ALL_AVG_FST_EW{$chr}}) {
### activate in order to skip sites that are not segr in current pool
#      next snp unless (exists ($SEGREGATING{$chr}{$pos}));
      $begin = $win*int($pos/$win); #jump to window close to snp
      for ($i=0; $i<(130000000-$win); $i+=$slide){ #slide through a few windows
#print STDERR "      for ($i=$begin; $i<($MAX_POS{$chr}-$win); $i+=$slide){ #slide through a few windows\n";
#         print "$i\t$pos\n";
#         print STDERR " if ($pos>$i && ($pos<($i+$win))){\n";
          if ($pos>$i && ($pos<($i+$win))){
#             print "HALLO\n";
              $win_pos = $i+($win/2);
	      foreach $breed (sort keys %{$ALL_AVG_FST_EW{$chr}{$pos}}) { #Each bree
		  $WINDOW_FST{$chr}{$win_pos}{$breed}+=$ALL_AVG_FST_EW{$chr}{$pos}{$breed};#Add up FST
		  $VARIANT_COUNT{$chr}{$win_pos}{$breed}++;
#	    print "$start_pos\t$pos\t$breed\t$WINDOW_FST{$breed}\t$variants\n";
	      }	      
          }
          #####HOPPAR TILL NYTT FONSTER INNAN SNPEN KOMMIT INNANFOR GRANSEN!!!
          next snp if ($i>$pos); #go to next snp if window start site has gone past snp position
      }
  }
}

=pod
$win_size=10000;
$start_flag=$variants=0;
foreach $pos (sort numerically keys %ALL_AVG_FST_EW) { #Each variant
    if ($start_flag==0) {
	$start_pos=$pos;
    }
    $start_flag=1;
    if (($pos-$start_pos)<$win_size) { #Which window
	$variants++;
	foreach $breed (sort keys %{$ALL_AVG_FST_EW{$pos}}) { #Each bree
	    $WINDOW_FST{$breed}+=$ALL_AVG_FST_EW{$pos}{$breed};#Add up FST
#	    print "$start_pos\t$pos\t$breed\t$WINDOW_FST{$breed}\t$variants\n";
	}
    }else{
	foreach $breed (sort keys %WINDOW_FST) {
	    $middle=$start_pos+($win_size/2);
	    $MOVING_FST{$middle}{$breed}=$WINDOW_FST{$breed}/$variants;
#	    print "!!!!!!!!!$middle\t$breed\t$MOVING_FST{$middle}{$breed}\n";
	    $WINDOW_FST{$breed}=0;
	    $WINDOW_FST{$breed}+=$ALL_AVG_FST_EW{$pos}{$breed};
	}
	$variants=1;
	$start_pos=$pos;
    }
}
=cut

=pod



$win_size=10000;
$start_flag=$variants=$first_end_passed=0;

foreach $pos (sort numerically keys %ALL_AVG_FST_EW) { #Each variant
    if ($start_flag==0) {
	$start_pos=$pos;
    }
    $start_flag=1;
    if (($pos-$start_pos)<$win_size) { #Which window
	$variants_A++;
	foreach $breed (sort keys %{$ALL_AVG_FST_EW{$pos}}) { #Each bree
	    $WINDOW_FST_A{$breed}+=$ALL_AVG_FST_EW{$pos}{$breed};#Add up FST
#	    print "$start_pos\t$pos\t$breed\t$WINDOW_FST{$breed}\t$variants\n";
	    if ($first_end_passed==0){ 
		if (($pos-$start_pos)>($win_size/2)) {
		    $WINDOW_FST_B{$breed}+=$ALL_AVG_FST_EW{$pos}{$breed};#Add up FST for overlapping window
		    $variants_B++;
		}
	    }else{
		if (($pos-$start_pos)<($win_size/2)){
		    $WINDOW_FST_B{$breed}+=$ALL_AVG_FST_EW{$pos}{$breed};#Add up FST for overlapping window
		    $variants_B++;
		}else{
		    foreach $breed (sort keys %WINDOW_FST_B) {
			$middle_B=$start_pos;
			$first_end_passed=0;
			$MOVING_FST{$middle_B}{$breed}=$WINDOW_FST_B{$breed}/$variants_B;
#	    print "!!!!!!!!!$middle\t$breed\t$MOVING_FST{$middle}{$breed}\n";
			$WINDOW_FST_B{$breed}=0;
			$WINDOW_FST_B{$breed}+=$ALL_AVG_FST_EW{$pos}{$breed};
			$variants_B=1;
			$first_end_passed=0;
		    }
		    
		}
	    }
	}
    }else{
	foreach $breed (sort keys %WINDOW_FST_A) {
	    $middle=$start_pos+($win_size/2);
	    $first_end_passed=1;
	    $MOVING_FST{$middle}{$breed}=$WINDOW_FST_A{$breed}/$variants_A;
#	    print "!!!!!!!!!$middle\t$breed\t$MOVING_FST{$middle}{$breed}\n";
	    $WINDOW_FST_A{$breed}=0;
	    $WINDOW_FST_A{$breed}+=$ALL_AVG_FST_EW{$pos}{$breed};
	}
	$variants_A=1;
	$start_pos=$pos;
    }
}

=cut


$average_file="/proj/snic2020-6-127/nobackup/private/private/160DG_MEM_variantcalling/ANNOTATED_VARIANTS/".$chromosome."_moving_pairwise_average_FST.txt";
open OUT, ">>$average_file";
print OUT "chr\tposition\t";
$stop=0;
foreach $chr (sort numerically keys %WINDOW_FST) {
    foreach $position (sort numerically keys %{$WINDOW_FST{$chr}}) {
	foreach $breed (sort keys %{$WINDOW_FST{$chr}{$position}}) {
	    print OUT "$breed\t";
	}
	$stop=1;
	last if($stop==1);
    }
    print OUT "\n";
}

foreach $chr (sort numerically keys %WINDOW_FST) {
    foreach $position (sort numerically keys %{$WINDOW_FST{$chr}}) {
	print OUT "$chr\t$position\t";
	foreach $breed (sort keys %{$WINDOW_FST{$chr}{$position}}) {
	    if ($VARIANT_COUNT{$chr}{$position}{$breed}>4){
		$moving_average_fst=$WINDOW_FST{$chr}{$position}{$breed}/$VARIANT_COUNT{$chr}{$position}{$breed};
		print OUT "$moving_average_fst\t";
	    }else{
		print OUT "NA\t";
	    }
	}
	print OUT "\n";
    }
}	    
close OUT;

#foreach $breed (sort keys %SAMPLE_ORDER) {
#    print "$breed\t@{$SAMPLE_ORDER{$breed}}\n";
#}
=pod
open OUT, ">>individual_depth_per_variable_site.txt";
foreach $depth (sort numerically keys %DEPTH) {
    print OUT "$depth\t$DEPTH{$depth}\n";
}

foreach $sample (@BREED_ORDER) {
    print "$sample\t";
}


foreach $breed (sort keys %ALLELE_FREQS_BIN) {
    print "$breed\n";
    $outfile = "/proj/b2013119/private/Analyses/Chr1_allele_freqs_".$breed.".txt"; 
    open OUT, ">>$outfile" or die;
    foreach $freq (sort numerically keys %{$ALLELE_FREQS_BIN{$breed}}) {
          print OUT "$breed\t$freq\t$ALLELE_FREQS_BIN{$breed}{$freq}\n";
    }
    close OUT or die;
}
=cut
sub numerically {
    $a<=>$b;
}







#9 first columns are info, 10 and onwards are samples
#GT:AD:DP:GQ:PL
#0/0:6,0:6:18:0,18,202
#1/1:0,3:3:9:105,9,0 
#./.
#Genotype:Allelic depths:Read deapth:Genotype quality:Genotype likelihood


#foreach $column (@header) {
#    print "$column\n";
#}

##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">


##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency, for each ALT allele, in the same order as listed">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
##INFO=<ID=BaseQRankSum,Number=1,Type=Float,Description="Z-score from Wilcoxon rank sum test of Alt Vs

#. Ref base qualities">
##INFO=<ID=ClippingRankSum,Number=1,Type=Float,Description="Z-score From Wilcoxon rank sum test of Alt vs. Ref number of hard clipped bases">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">
##INFO=<ID=DS,Number=0,Type=Flag,Description="Were any of the samples downsampled?">
##INFO=<ID=EFF,Number=.,Type=String,Description="Predicted effects for this variant.Format: 'Effect ( Effect_Impact | Functional_Class | Codon_Change | Amino_Acid_Change| Amino_Acid_length | Gene_Name | Transcript_BioType | Gene_Coding | Transcript_ID | Exon_Rank  | Genotype_Number [ | ERRORS | WARNINGS ] )'">
##INFO=<ID=END,Number=1,Type=Integer,Description="Stop position of the interval">
##INFO=<ID=FS,Number=1,Type=Float,Description="Phred-scaled p-value using Fisher's exact test to detect strand bias">
##INFO=<ID=HaplotypeScore,Number=1,Type=Float,Description="Consistency of the site with at most two segregating haplotypes">
##INFO=<ID=InbreedingCoeff,Number=1,Type=Float,Description="Inbreeding coefficient as estimated from the genotype likelihoods per-sample when compared against the Hardy-Weinberg expectation">
##INFO=<ID=MLEAC,Number=A,Type=Integer,Description="Maximum likelihood expectation (MLE) for the allele counts (not necessarily the same as the AC), for each ALT allele, in the same order as listed">
##INFO=<ID=MLEAF,Number=A,Type=Float,Description="Maximum likelihood expectation (MLE) for the allele frequency (not necessarily the same as the AF), for each ALT allele, in the same order as listed">
##INFO=<ID=MQ,Number=1,Type=Float,Description="RMS Mapping Quality">
##INFO=<ID=MQ0,Number=1,Type=Integer,Description="Total Mapping Quality Zero Reads">
##INFO=<ID=MQRankSum,Number=1,Type=Float,Description="Z-score From Wilcoxon rank sum test of Alt vs. Ref read mapping qualities">
##INFO=<ID=NEGATIVE_TRAIN_SITE,Number=0,Type=Flag,Description="This variant was used to build the negative training set of bad variants">
##INFO=<ID=POSITIVE_TRAIN_SITE,Number=0,Type=Flag,Description="This variant was used to build the positive training set of good variants">
##INFO=<ID=QD,Number=1,Type=Float,Description="Variant Confidence/Quality by Depth">
##INFO=<ID=ReadPosRankSum,Number=1,Type=Float,Description="Z-score from Wilcoxon rank sum test of Alt vs. Ref read position bias">
##INFO=<ID=VQSLOD,Number=1,Type=Float,Description="Log odds ratio of being a true variant versus being false under the trained gaussian mixture model">
##INFO=<ID=culprit,Number=1,Type=String,Description="The annotation which was the worst performing in the Gaussian mixture model, likely the reason why the variant was filtered out">

