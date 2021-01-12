#! /usr/bin/perl -w


use Data::Dumper;

$min_sz=20; #minimum number of chromosomes sequenced per breed to include site
#$cons_pp100=1;
#$noncons_pp100=-1;
$blocks=50; #number of jack knif blocks
$outgroup='andean'; #set out group (cat 'cat' or andean fox 'andean')
$variant_type='SNPs_and_INDELs'; #'SNPs' for SNP based analysis and 'SNPs and INDELs' for analysis based on both types of sites

@breeds=('beagle','ckcs','gs','gr','lr','poodle','rotw','whwt');
#@breeds=('beagle','ckcs');

#results for pairwise (breeds) relative quantity of different classses of derived sites. 
$results_file="Purging_and_negative_selection_all_breeds_outgroup_".$outgroup."_jackknifewindows_".$blocks."_sites_".$variant_type."_pp0_pp7_autosames_210112.txt";

#results for pairwise (breeds) relative quantity of different classses of derived sites corrected for mutation rate differences. 
$rate_cor_results_file="Purging_and_negative_selection_all_breeds_outgroup_".$outgroup."_jackknifewindows_".$blocks."_sites_".$variant_type."_pp0_pp7_autosames_rate_corrected_210112.txt";

#results for pairwise (breeds) relative homozygosity of different classses of derived sites. 
$rresults_file="Derived_homozygosity_all_breeds_outgroup_".$outgroup."_jackknifewindows_".$blocks."_sites_".$variant_type."_pp0_pp7_autosames_210112.txt";

### Estimate genome size ###

$prev_pos=$prev_chr=$accum_genome_size=0;
for (</proj/snic2020-6-127/nobackup/private/private/160DG_MEM_variantcalling/ANNOTATED_VARIANTS/FST*all_variants_sz0_annotated.txt>) {
#for (</proj/b2013119/nobackup/private/160DG_MEM_variantcalling/FST*all_variants_sz0_annotated.txt>) {
#for (</proj/uppstore2017236/b2013119_nobackup/private/160DG_MEM_variantcalling/ANNOTATED_VARIANTS/FST_chr38_all_variants_sz0_annotated.txt>) {
    next if (/chrM/);
    next if (/chrX/);
    next if (/chrUn/);
    if (/chr(\d+)/) {
	$chr_test=$1;
    }
#print STDERR    "next if ($chr_test<37);\n";
    next if ($chr_test<0);
    $variant_file=$_;
    print "$variant_file\n";
    open IN, $variant_file;
  line:  while (<IN>) {
      $count++;
      next line unless ($count>1);
      $gene = $eff = $high_eff =  $gene_name= "NA";
      @sv = split, /\t/;
#      $koll=@sv;
#      print "$koll\n";
#      for ($i=0; $i<$koll; $i++) {
#       print "$i\t$sv[$i]\n";
#      }
      $chr=$sv[0];
#    $chr =~ s/chr//;
      $pos=$sv[1];
      next line if ($pos eq 'pos'); #ship header
      if ($prev_pos==0) {
	  $accum_genome_size=$pos;
      }else{
	  if ($chr == $prev_chr) {
	      $accum_genome_size=$accum_genome_size+($pos-$prev_pos);
	  }
      }
      $prev_pos=$pos;
      $prev_chr=$chr;
  }
#    print "$accum_genome_size\t$chr\t$pos\n";

}
print "$accum_genome_size\t$chr\t$pos\n";
#die;






### Start jack knifing: pre define blocks ###

$prev_pos=$prev_chr=$accum_pos=$leave_flag=$old_block=0;

for ($jack_block=1; $jack_block<=$blocks; $jack_block++) {
    $block_start=($jack_block-1)*$accum_genome_size/$blocks;
    $block_end=($jack_block*$accum_genome_size/$blocks);
    $JACK_BLOCK_STARTS{$jack_block}=$block_start;
    $JACK_BLOCK_ENDS{$jack_block}=$block_end;
#    $l_c_b_syn{$jack_block}=$l_b_c_syn{$jack_block}=$l_c_b_nonsyn{$jack_block}=$l_b_c_nonsyn{$jack_block}=$l_c_b_lof{$jack_block}=$l_b_c_lof{$jack_block}=0;
}



 #  %L_C_B_CONS = %L_B_C_CONS = ();
 #   $l_c_b_syn=$l_b_c_syn=$l_c_b_nonsyn=$l_b_c_nonsyn=$l_c_b_lof=$l_b_c_lof=0;
 #   $total_sites=0;
    


### load variants  ###
$count=0;
for (</proj/snic2020-6-127/nobackup/private/private/160DG_MEM_variantcalling/ANNOTATED_VARIANTS/FST*all_variants_sz0_annotated.txt>) {
#for (</proj/b2013119/nobackup/private/160DG_MEM_variantcalling/FST*all_variants_sz0_annotated.txt>) {
#for (</proj/uppstore2017236/b2013119_nobackup/private/160DG_MEM_variantcalling/ANNOTATED_VARIANTS/FST_chr38_all_variants_sz0_annotated.txt>) {
    next if (/chrM/);
    next if (/chrX/);
     next if (/chrUn/);
  next if (/chrUn/);
    if (/chr(\d+)/) {
	$chr_test=$1;
    }
#print STDERR    "next if ($chr_test<37);\n";
    next if ($chr_test<0);
    $variant_file=$_;
    print "$variant_file\n";
    open IN, $variant_file;
  line:  while (<IN>) {
      $count++;
      next line unless ($count>1);
      $gene = $eff = $high_eff =  $gene_name= "NA";
      @sv = split, /\t/;
#      $koll=@sv;
#      print "$koll\n";
#      for ($i=0; $i<$koll; $i++) {
#       print "$i\t$sv[$i]\n";
#      }
      $chr=$sv[0];
#    $chr =~ s/chr//;
      $pos=$sv[1];
 
     ### Potentially insert code that skips positions in regions that show evidence of selection in CKCS here ###

      next line if ($pos eq 'pos'); #ship header
      if ($prev_pos==0) {
	  $accum_pos=0;
      }else{
	  if ($chr == $prev_chr) {
	      $accum_pos=$accum_pos+($pos-$prev_pos);
	  }
      }
      $prev_pos=$pos;
      $prev_chr=$chr;
 #     print "$accum_pos\t$chr\t$pos\t";
      foreach $block (sort numerically keys %JACK_BLOCK_STARTS) {
	  if (($JACK_BLOCK_STARTS{$block}<$accum_pos) && ($accum_pos<$JACK_BLOCK_ENDS{$block})) {
	      $surrent_block=$block;
	      unless($surrent_block==$old_block) {
		  print "NEW_BLOCK:$surrent_block\t$chr\t$pos\n";
	      }
#	      print "$chr\t$pos\t$accum_pos\t$surrent_block\n"; 
	  }
      }
      $old_block=$surrent_block;
#      print "$surrent_block\n";

#      ###jack knife proceedure leaves one segment of genome out for every iteration ###
#      if (($accum_pos>$block_start) && ($accum_pos<$block_end)){
#	  if ($leave_flag==0) {
#	      print "start blocking $chr $pos $accum_pos\n";
#	      $leave_flag=1;
#	  }
#	  next line;
#     }
#      if ($leave_flag==1) {
#	  print "end blocking $chr $pos $accum_pos\n";
#	  $leave_flag=0;
#      }

      $ref=$sv[4];
      $alt=$sv[5];
      $cat=$sv[6];
#      print "$pos\t$ref\t$alt\t$cat\n";

 ###### only use biallelic SNPS #####
      if ($variant_type eq 'SNPs') {
      next line if (length($ref) >1);
      next line if (length($alt) >1);
      }
      
 
##### include all biallelic variants (including indels) ####
      if ($variant_type eq 'SNPs_and_INDELs') {
      next line if ($ref=~/\,/);
      next line if ($alt=~/\,/);
      }

      
      $andean=$sv[7];
      $wolf_alt=$sv[8];
      $eff=$sv[10];
      $phylo_fortysix=$sv[13];
      $phylo_hundred=$sv[12];
      #  $gene=$sv[5];
      $beagle_fst=$sv[14];
      $ckcs_fst=$sv[15];
      $gs_fst=$sv[16];
      $gr_fst=$sv[17];
      $lr_fst=$sv[18];
      $poodle_fst=$sv[19];
      $rotw_fst=$sv[20];
      $whwt_fst=$sv[21];
      $beagle_sz=$sv[26];
      $beagle_A=$sv[24];
      $beagle_B=$sv[25];
      $ckcs_sz=$sv[29];
      $ckcs_A=$sv[27];
      $ckcs_B=$sv[28];
      $gs_sz=$sv[32];
      $gs_A=$sv[30];
      $gs_B=$sv[31];
      $gr_sz=$sv[35];
      $gr_A=$sv[33];
      $gr_B=$sv[34];
      $lr_sz=$sv[38];
      $lr_A=$sv[36];
      $lr_B=$sv[37];
      $poodle_sz=$sv[41];
      $poodle_A=$sv[39];
      $poodle_B=$sv[40];
      $rotw_sz=$sv[44];
      $rotw_A=$sv[42];
      $rotw_B=$sv[43];
      
      $whwt_sz=$sv[47];
      $whwt_A=$sv[45];
      $whwt_B=$sv[46];
      $total_fst=$sv[22];
      $PHYLOP{$chr}{$pos}=$phylo_fortysix;
      
      ### only include sites with know ancestral state ###

      # use cat as outgroup
      if ($outgroup eq 'cat') {
	  next line unless ((lc $cat eq lc $alt) || (lc $cat eq lc $ref));
	  if (lc $cat eq lc $alt) { # dog reference is derived state
	      $DERIVED_AF{'beagle'}=$sv[24];
	      $DERIVED_AF{'ckcs'}=$sv[27];
	      $DERIVED_AF{'gs'}=$sv[30];
	      $DERIVED_AF{'gr'}=$sv[33];
	      $DERIVED_AF{'lr'}=$sv[36];
	      $DERIVED_AF{'poodle'}=$sv[39];
	      $DERIVED_AF{'rotw'}=$sv[42];
	      $DERIVED_AF{'whwt'}=$sv[45];
	  }
	  
	  if (lc $cat eq lc $ref) { #dog reference is ancestral state
	      $DERIVED_AF{'beagle'}=$sv[25];
	      $DERIVED_AF{'ckcs'}=$sv[28];
	      $DERIVED_AF{'gs'}=$sv[31];
	      $DERIVED_AF{'gr'}=$sv[34];
	      $DERIVED_AF{'lr'}=$sv[37];
	      $DERIVED_AF{'poodle'}=$sv[40];
	      $DERIVED_AF{'rotw'}=$sv[43];
	      $DERIVED_AF{'whwt'}=$sv[46];
	  }

      }


      # use andean fox as out group
      if ($outgroup eq 'andean') {
	  next line unless ((lc $andean eq lc $alt) || (lc $andean eq lc $ref));
	  if (lc $andean eq lc $alt) { # dog reference is derived state
	      $DERIVED_AF{'beagle'}=$sv[24];
	      $DERIVED_AF{'ckcs'}=$sv[27];
	      $DERIVED_AF{'gs'}=$sv[30];
	      $DERIVED_AF{'gr'}=$sv[33];
	      $DERIVED_AF{'lr'}=$sv[36];
	      $DERIVED_AF{'poodle'}=$sv[39];
	      $DERIVED_AF{'rotw'}=$sv[42];
	      $DERIVED_AF{'whwt'}=$sv[45];
	  }
	  
	  if (lc $andean eq lc $ref) { #dog reference is ancestral state
	      $DERIVED_AF{'beagle'}=$sv[25];
	      $DERIVED_AF{'ckcs'}=$sv[28];
	      $DERIVED_AF{'gs'}=$sv[31];
	      $DERIVED_AF{'gr'}=$sv[34];
	      $DERIVED_AF{'lr'}=$sv[37];
	      $DERIVED_AF{'poodle'}=$sv[40];
	      $DERIVED_AF{'rotw'}=$sv[43];
	      $DERIVED_AF{'whwt'}=$sv[46];
	  }
      }


      ### filter for variants with high sample size ###
      #  if ($total_fst>$min_fst || $beagle_fst>$min_fst || $ckcs_fst>$min_fst || $gs_fst>$min_fst || $gr_fst>$min_fst || $lr_fst>$min_fst || $poodle_fst>$min_fst || $rotw_fst>$min_fst || $whwt_fst>$min_fst) {
      if ($beagle_sz>$min_sz && $ckcs_sz>$min_sz && $gs_sz>$min_sz && $gr_sz>$min_sz && $lr_sz>$min_sz && $poodle_sz>$min_sz && $rotw_sz>$min_sz && $whwt_sz>$min_sz) {
	  $total_sites{$surrent_block}++;
	  
	      $SZ{'beagle'}=$beagle_sz;
	      $SZ{'ckcs'}=$ckcs_sz;
	      $SZ{'gs'}=$gs_sz;
	      $SZ{'gr'}=$gr_sz;
	      $SZ{'lr'}=$lr_sz;
	      $SZ{'poodle'}=$poodle_sz;
	      $SZ{'rotw'}=$rotw_sz;
	      $SZ{'whwt'}=$whwt_sz;

###### Test for purging of strongly deleterious mutation in CKCS and Beagle #####
	  
	  
# Contrast selection and purging at synonymous, nonsynonymous and LOF sites 
	  
	  
#	  if (lc $cat eq lc $alt) { # dog reference is derived state

	  foreach $breed_1 (sort @breeds) {
	      foreach $breed_2 (sort @breeds) {
		  unless ($breed_1 eq $breed_2) {
		      if ($DERIVED_AF{$breed_1}>0) { #only process sites that are variable in breed 1 as these are the only informative sites
#		      print "$breed_1\t$breed_2";
#		      print "\t$DERIVED_AF{$breed_1}\t$DERIVED_AF{$breed_2}\t";
			  $l=$DERIVED_AF{$breed_1}*(1-$DERIVED_AF{$breed_2});

			  #prob homozygous breed 1 not breed2
			  $ll=((($DERIVED_AF{$breed_1}*$SZ{$breed_1})**2-($DERIVED_AF{$breed_1}*$SZ{$breed_1}))/($SZ{$breed_1}*($SZ{$breed_1}-1))*(1-(($DERIVED_AF{$breed_2}*$SZ{$breed_2})**2-($DERIVED_AF{$breed_2}*$SZ{$breed_2}))/($SZ{$breed_2}*($SZ{$breed_2}-1))));

			  ### Nonsynonymous mutations
			  if ($eff eq "NON_SYNONYMOUS_CODING") {
			      $L_NONSYN{$surrent_block}{$breed_1}{$breed_2}+=$l;
			      $TOTAL_L_NONSYN{$breed_1}{$breed_2}+=$l;
			      $BLOCK_NONSYN_INF_SITES{$surrent_block}{$breed_1}{$breed_2}++;
			      $TOTAL_NONSYN_INF_SITES{$breed_1}{$breed_2}++;

			      #prob homozygous breed 1 not breed2
			      $LL_NONSYN{$surrent_block}{$breed_1}{$breed_2}+=$ll;
			      $TOTAL_LL_NONSYN{$breed_1}{$breed_2}+=$ll;
			  }
			  
			  ### Synonymous mutations
			  if ($eff eq "SYNONYMOUS_CODING") {
			      $L_SYN{$surrent_block}{$breed_1}{$breed_2}+=$l;
			      $TOTAL_L_SYN{$breed_1}{$breed_2}+=$l;
			      $BLOCK_SYN_INF_SITES{$surrent_block}{$breed_1}{$breed_2}++;
			      $TOTAL_SYN_INF_SITES{$breed_1}{$breed_2}++;

			      #prob homozygous breed 1 not breed2
			      $LL_SYN{$surrent_block}{$breed_1}{$breed_2}+=$ll;
			      $TOTAL_LL_SYN{$breed_1}{$breed_2}+=$ll;

			  }
			  
			  ### Loss of function mutations
			  if (($eff eq 'CODON_CHANGE_PLUS_CODON_DELETION') or ($eff eq 'CODON_CHANGE_PLUS_CODON_INSERTION') or ($eff eq 'CODON_DELETION') or ($eff eq 'CODON_INSERTION') or ($eff eq 'EXON_DELETED') or ($eff eq "FRAME_SHIFT") or ($eff eq 'SPLICE_SITE_ACCEPTOR') or ($eff eq 'SPLICE_SITE_DONOR') or ($eff eq 'START_LOST') or ($eff eq 'STOP_GAINED')) {
			      
#'CODON_CHANGE_PLUS_CODON_DELETION','CODON_CHANGE_PLUS_CODON_INSERTION','CODON_DELETION','CODON_INSERTION','EXON_DELETED','FRAME_SHIFT','SPLICE_SITE_ACCEPTOR','SPLICE_SITE_DONOR','START_LOST','STOP_GAINED',
#			  print STDERR "$breed_1\t$breed_2\t			  $L_LOF{$surrent_block}{$breed_1}{$breed_2}+=$l;
#\n";
			      $L_LOF{$surrent_block}{$breed_1}{$breed_2}+=$l;
			      $TOTAL_L_LOF{$breed_1}{$breed_2}+=$l;
			      $BLOCK_LOF_INF_SITES{$surrent_block}{$breed_1}{$breed_2}++;
			      $TOTAL_LOF_INF_SITES{$breed_1}{$breed_2}++;

			      #prob homozygous breed 1 not breed2
			      $LL_LOF{$surrent_block}{$breed_1}{$breed_2}+=$ll;
			      $TOTAL_LL_LOF{$breed_1}{$breed_2}+=$ll;
			  }
			  
			  unless ($phylo_hundred eq 'NA') {
			      if ($phylo_hundred>0) {
				  $L_CONS{$surrent_block}{$breed_1}{$breed_2}{cons_bin($phylo_hundred)}+=$l;
				  $TOTAL_L_CONS{$breed_1}{$breed_2}{cons_bin($phylo_hundred)}+=$l;
				  $BLOCK_CONS_INF_SITES{$surrent_block}{$breed_1}{$breed_2}{cons_bin($phylo_hundred)}++;
				  $TOTAL_CONS_INF_SITES{$breed_1}{$breed_2}{cons_bin($phylo_hundred)}++;
				  #prob homozygous breed 1 not breed2
				  $LL_CONS{$surrent_block}{$breed_1}{$breed_2}{cons_bin($phylo_hundred)}+=$ll;
				  $TOTAL_LL_CONS{$breed_1}{$breed_2}{cons_bin($phylo_hundred)}+=$ll;
				  
			      }
			  }
		      }
		  }
	      }
	  }
      }
  }			  
}



#Total (across entitre genome) comparison of purging and selection across breeds 
foreach $breed_1 (sort @breeds) {
    foreach $breed_2 (sort @breeds) {
	unless ($breed_1 eq $breed_2) {
#		unless (exists $TOTAL_R_SYN{$breed_2}{$breed_1}) {
	    print "LOF\t$breed_1\t$breed_2\ttotal informative sites: $TOTAL_LOF_INF_SITES{$breed_1}{$breed_2}\n";
	    $TOTAL_R_LOF{$breed_1}{$breed_2}=$TOTAL_L_LOF{$breed_1}{$breed_2}/$TOTAL_L_LOF{$breed_2}{$breed_1};
	    print STDERR "$TOTAL_R_LOF{$breed_1}{$breed_2}=$TOTAL_L_LOF{$breed_1}{$breed_2}/$TOTAL_L_LOF{$breed_2}{$breed_1};\n";
	    $TOTAL_R_SYN{$breed_1}{$breed_2}=$TOTAL_L_SYN{$breed_1}{$breed_2}/$TOTAL_L_SYN{$breed_2}{$breed_1};
	    $TOTAL_R_NONSYN{$breed_1}{$breed_2}=$TOTAL_L_NONSYN{$breed_1}{$breed_2}/$TOTAL_L_NONSYN{$breed_2}{$breed_1};


	    #prob homozygous breed 1 not breed2
	    $TOTAL_RR_LOF{$breed_1}{$breed_2}=$TOTAL_LL_LOF{$breed_1}{$breed_2}/$TOTAL_LL_LOF{$breed_2}{$breed_1};
	    $TOTAL_RR_SYN{$breed_1}{$breed_2}=$TOTAL_LL_SYN{$breed_1}{$breed_2}/$TOTAL_LL_SYN{$breed_2}{$breed_1};
	    $TOTAL_RR_NONSYN{$breed_1}{$breed_2}=$TOTAL_LL_NONSYN{$breed_1}{$breed_2}/$TOTAL_LL_NONSYN{$breed_2}{$breed_1};
	    foreach $bin (sort numerically keys %{$TOTAL_L_CONS{$breed_1}{$breed_2}}) {
		if (($TOTAL_L_CONS{$breed_1}{$breed_2}{$bin}>0) and ($TOTAL_L_CONS{$breed_2}{$breed_1}{$bin}>0)){
		    print "CONS\t$breed_1\t$breed_2\t$bin\ttotal informative sites: $TOTAL_CONS_INF_SITES{$breed_1}{$breed_2}{$bin}\n";
		    
		    $TOTAL_R_CONS{$breed_1}{$breed_2}{$bin}=$TOTAL_L_CONS{$breed_1}{$breed_2}{$bin}/$TOTAL_L_CONS{$breed_2}{$breed_1}{$bin};
		    print STDERR "$TOTAL_R_CONS{$breed_1}{$breed_2}{$bin}=$TOTAL_L_CONS{$breed_1}{$breed_2}{$bin}/$TOTAL_L_CONS{$breed_2}{$breed_1}{$bin};\n";

		    #prob homozygous breed 1 not breed2
		    $TOTAL_RR_CONS{$breed_1}{$breed_2}{$bin}=$TOTAL_LL_CONS{$breed_1}{$breed_2}{$bin}/$TOTAL_LL_CONS{$breed_2}{$breed_1}{$bin};
		}
	    }
	}
    }
}



print "TOTAL genome\n";
foreach $breed_1 (sort @breeds) {
    foreach $breed_2 (sort @breeds) {
	unless ($breed_1 eq $breed_2) {
#	    unless (exists $TOTAL_R_SYN{$breed_2}{$breed_1}) {
	    print "$breed_1\t$breed_2:\n";
	    print "LOF\tR\t$TOTAL_R_LOF{$breed_1}{$breed_2}\tINF_SITES\t$TOTAL_LOF_INF_SITES{$breed_1}{$breed_2}\n";
	    print "SYN\t$TOTAL_R_SYN{$breed_1}{$breed_2}\tINF_SITES\t$TOTAL_SYN_INF_SITES{$breed_1}{$breed_2}\n";
	    print "NONSYN\t$TOTAL_R_NONSYN{$breed_1}{$breed_2}\tINF_SITES\t$TOTAL_NONSYN_INF_SITES{$breed_1}{$breed_2}\n";
	    foreach $bin (sort numerically keys %{$TOTAL_L_CONS{$breed_1}{$breed_2}}) {
		print "CONS\t$bin\tR\t$TOTAL_R_CONS{$breed_1}{$breed_2}{$bin}\tINF_SITES\t$TOTAL_CONS_INF_SITES{$breed_1}{$breed_2}{$bin}\n";
	    }
	}
    }
}
#}


#estimate total L (%JACK_TOTAL_L) for all blocks when leaving out the ith observation
#estimate average L (%JACK_AVG_L) for all blocks when leaving out the ith observation
for ($pseudo_replicate=1; $pseudo_replicate<=$blocks; $pseudo_replicate++) {
    foreach $block (sort numerically keys %JACK_BLOCK_STARTS) {
	unless ($block == $pseudo_replicate) {
	    foreach $breed_1 (sort @breeds) {
		foreach $breed_2 (sort @breeds) {
		    unless ($breed_1 eq $breed_2) {
#			unless (exists $OLD_JACK_AVG_R_SYN{$pseudo_replicate}{$breed_2}{$breed_1}) {
#			$OLD_JACK_AVG_R_LOF{$pseudo_replicate}{$breed_1}{$breed_2}+=$R_LOF{$block}{$breed_1}{$breed_2}/($blocks-1); 
#			$JACK_AVG_L_LOF{$pseudo_replicate}{$breed_1}{$breed_2}+=$L_LOF{$block}{$breed_1}{$breed_2}/($blocks-1);
			$JACK_TOTAL_L_LOF{$pseudo_replicate}{$breed_1}{$breed_2}+=$L_LOF{$block}{$breed_1}{$breed_2};
			
#			$WEIGHTED_JACK_AVG_R_LOF{$pseudo_replicate}{$breed_1}{$breed_2}+=(1-$BLOCK_LOF_INF_SITES{$block}{$breed_1}{$breed_2}/$TOTAL_LOF_INF_SITES{$breed_1}{$breed_2})*$R_LOF{$block}{$breed_1}{$breed_2}/($blocks-1); 
			
#			if (($breed_1 eq 'ckcs') && ($breed_2 eq 'beagle')) {
#			    print "leave out\t$pseudo_replicate\tblock\t$block\t$breed_1\t$breed_2\t$R_LOF{$block}{$breed_1}{$breed_2}\t$OLD_JACK_AVG_R_LOF{$pseudo_replicate}{$breed_1}{$breed_2}\n";
#			    print STDERR "			$OLD_JACK_AVG_R_LOF{$pseudo_replicate}{$breed_1}{$breed_2}+=$R_LOF{$block}{$breed_1}{$breed_2}/($blocks-1); \n";
#			}
#			$AVG_OF_JACK_AVG_R_LOF{$breed_1}{$breed_2}+=($R_LOF{$block}{$breed_1}{$breed_2}/($blocks-1))/$blocks;
			
#			$OLD_JACK_AVG_R_SYN{$pseudo_replicate}{$breed_1}{$breed_2}+=$R_SYN{$block}{$breed_1}{$breed_2}/($blocks-1);
#			$JACK_AVG_L_SYN{$pseudo_replicate}{$breed_1}{$breed_2}+=$L_SYN{$block}{$breed_1}{$breed_2}/($blocks-1);
			$JACK_TOTAL_L_SYN{$pseudo_replicate}{$breed_1}{$breed_2}+=$L_SYN{$block}{$breed_1}{$breed_2};
#			$WEIGHTED_JACK_AVG_R_SYN{$pseudo_replicate}{$breed_1}{$breed_2}+=(1-$BLOCK_SYN_INF_SITES{$block}{$breed_1}{$breed_2}/$TOTAL_SYN_INF_SITES{$breed_1}{$breed_2})*$R_SYN{$block}{$breed_1}{$breed_2}/($blocks-1); 
			
#			$AVG_OF_JACK_AVG_R_SYN{$breed_1}{$breed_2}+=($R_SYN{$block}{$breed_1}{$breed_2}/($blocks-1))/$blocks;
			
#			$OLD_JACK_AVG_R_NONSYN{$pseudo_replicate}{$breed_1}{$breed_2}+=$R_NONSYN{$block}{$breed_1}{$breed_2}/($blocks-1);
#			$JACK_AVG_L_NONSYN{$pseudo_replicate}{$breed_1}{$breed_2}+=$L_NONSYN{$block}{$breed_1}{$breed_2}/($blocks-1);
			$JACK_TOTAL_L_NONSYN{$pseudo_replicate}{$breed_1}{$breed_2}+=$L_NONSYN{$block}{$breed_1}{$breed_2};
#			$WEIGHTED_JACK_AVG_R_NONSYN{$pseudo_replicate}{$breed_1}{$breed_2}+=(1-$BLOCK_NONSYN_INF_SITES{$block}{$breed_1}{$breed_2}/$TOTAL_NONSYN_INF_SITES{$breed_1}{$breed_2})*$R_NONSYN{$block}{$breed_1}{$breed_2}/($blocks-1); 
			
#			$AVG_OF_JACK_AVG_R_NONSYN{$breed_1}{$breed_2}+=($R_NONSYN{$block}{$breed_1}{$breed_2}/($blocks-1))/$blocks;

			#prob homozygous breed 1 not breed2
			$JACK_TOTAL_LL_LOF{$pseudo_replicate}{$breed_1}{$breed_2}+=$LL_LOF{$block}{$breed_1}{$breed_2};
			$JACK_TOTAL_LL_SYN{$pseudo_replicate}{$breed_1}{$breed_2}+=$LL_SYN{$block}{$breed_1}{$breed_2};
			$JACK_TOTAL_LL_NONSYN{$pseudo_replicate}{$breed_1}{$breed_2}+=$LL_NONSYN{$block}{$breed_1}{$breed_2};

			
			foreach $bin (sort numerically keys %{$L_CONS{$block}{$breed_1}{$breed_2}}) {
#			    print "$leave_out\t$block\t$breed_1\t$breed_2\t$bin\t";
#			    print STDERR "			    $JACK_AVG_R_CONS{$breed_1}{$breed_2}{$bin}+=($R_CONS{$block}{$breed_1}{$breed_2}{$bin}/($blocks-1));\n";
#			    $OLD_JACK_AVG_R_CONS{$pseudo_replicate}{$breed_1}{$breed_2}{$bin}+=($R_CONS{$block}{$breed_1}{$breed_2}{$bin}/($blocks-1));
#			    $JACK_AVG_L_CONS{$pseudo_replicate}{$breed_1}{$breed_2}{$bin}+=$L_CONS{$block}{$breed_1}{$breed_2}{$bin}/($blocks-1);
			    $JACK_TOTAL_L_CONS{$pseudo_replicate}{$breed_1}{$breed_2}{$bin}+=$L_CONS{$block}{$breed_1}{$breed_2}{$bin};

			    if ($L_CONS{$block}{$breed_1}{$breed_2}{$bin}>0) {
				$COUNT_AVG_R_CONS_NONZERO_BLOCKS{$breed_1}{$breed_2}{$bin}++;
			    }

			    #prob homozygous breed 1 not breed2
			    $JACK_TOTAL_LL_CONS{$pseudo_replicate}{$breed_1}{$breed_2}{$bin}+=$LL_CONS{$block}{$breed_1}{$breed_2}{$bin};
			   

#			    $WEIGHTED_JACK_AVG_R_CONS{$pseudo_replicate}{$breed_1}{$breed_2}{$bin}+=(1-$BLOCK_CONS_INF_SITES{$block}{$breed_1}{$breed_2}{$bin}/$TOTAL_CONS_INF_SITES{$breed_1}{$breed_2}{$bin})*$R_CONS{$block}{$breed_1}{$breed_2}{$bin}/($blocks-1); 
			    
#			    $AVG_OF_JACK_AVG_R_CONS{$breed_1}{$breed_2}{$bin}+=($R_CONS{$block}{$breed_1}{$breed_2}{$bin}/($blocks-1))/$blocks;
			    
#			    print "$JACK_AVG_R_CONS{$leave_out}{$breed_1}{$breed_2}{$bin}\n";
			}
		    }
		}
	    }
	}
    }
}




#estimate  average R for (%JACK_AVG_R) each individual pseudoreplicate by taking ratio of %JACK_AVG_L_breed1_breed2/%JACK_AVG_L_breed2_breed1

#estimate  R for (%JACK_R) each pseudoreplicate by taking ratio of %JACK_TOTAL_L_breed1_breed2/%JACK_TOTAL_L_breed2_breed1

#estimate R (%WEIGHTED_JACK_R) weighted by sample size for each pseudoreplicate 

foreach $pseudo_replicate (sort numerically keys %JACK_TOTAL_L_SYN) {  
    foreach $breed_1 (sort @breeds) {
	foreach $breed_2 (sort @breeds) {
	    unless ($breed_1 eq $breed_2) {
#		unless (exists $JACK_AVG_R_LOF{$breed_2}{$breed_1}) {
		    $JACK_R_LOF{$pseudo_replicate}{$breed_1}{$breed_2}=$JACK_TOTAL_L_LOF{$pseudo_replicate}{$breed_1}{$breed_2}/$JACK_TOTAL_L_LOF{$pseudo_replicate}{$breed_2}{$breed_1};
		    print STDERR "		    $JACK_R_LOF{$pseudo_replicate}{$breed_1}{$breed_2}=$JACK_TOTAL_L_LOF{$pseudo_replicate}{$breed_1}{$breed_2}/$JACK_TOTAL_L_LOF{$pseudo_replicate}{$breed_2}{$breed_1};\n";
		    $JACK_R_SYN{$pseudo_replicate}{$breed_1}{$breed_2}=$JACK_TOTAL_L_SYN{$pseudo_replicate}{$breed_1}{$breed_2}/$JACK_TOTAL_L_SYN{$pseudo_replicate}{$breed_2}{$breed_1};
		    $JACK_R_NONSYN{$pseudo_replicate}{$breed_1}{$breed_2}=$JACK_TOTAL_L_NONSYN{$pseudo_replicate}{$breed_1}{$breed_2}/$JACK_TOTAL_L_NONSYN{$pseudo_replicate}{$breed_2}{$breed_1};
		    
		    $h_LOF=$TOTAL_LOF_INF_SITES{$breed_1}{$breed_2}/$BLOCK_LOF_INF_SITES{$pseudo_replicate}{$breed_1}{$breed_2};
		    $h_SYN=$TOTAL_SYN_INF_SITES{$breed_1}{$breed_2}/$BLOCK_SYN_INF_SITES{$pseudo_replicate}{$breed_1}{$breed_2};
		    $h_NONSYN=$TOTAL_NONSYN_INF_SITES{$breed_1}{$breed_2}/$BLOCK_NONSYN_INF_SITES{$pseudo_replicate}{$breed_1}{$breed_2};

		    $WEIGHTED_JACK_R_LOF{$pseudo_replicate}{$breed_1}{$breed_2}=$h_LOF*$TOTAL_R_LOF{$breed_1}{$breed_2}-($h_LOF-1)*$JACK_R_LOF{$pseudo_replicate}{$breed_1}{$breed_2};
print STDERR "		    $WEIGHTED_JACK_R_LOF{$pseudo_replicate}{$breed_1}{$breed_2}=$h_LOF*$TOTAL_R_LOF{$breed_1}{$breed_2}-($h_LOF-1)*$JACK_R_LOF{$pseudo_replicate}{$breed_1}{$breed_2};\n";
		    $WEIGHTED_JACK_R_SYN{$pseudo_replicate}{$breed_1}{$breed_2}=$h_SYN*$TOTAL_R_SYN{$breed_1}{$breed_2}-($h_SYN-1)*$JACK_R_SYN{$pseudo_replicate}{$breed_1}{$breed_2};
		    $WEIGHTED_JACK_R_NONSYN{$pseudo_replicate}{$breed_1}{$breed_2}=$h_NONSYN*$TOTAL_R_NONSYN{$breed_1}{$breed_2}-($h_NONSYN-1)*$JACK_R_NONSYN{$pseudo_replicate}{$breed_1}{$breed_2};
		    
#		    $JACK_AVG_R_LOF{$pseudo_replicate}{$breed_1}{$breed_2}=$JACK_AVG_L_LOF{$pseudo_replicate}{$breed_1}{$breed_2}/$JACK_AVG_L_LOF{$pseudo_replicate}{$breed_2}{$breed_1};
#		    $JACK_AVG_R_SYN{$pseudo_replicate}{$breed_1}{$breed_2}=$JACK_AVG_L_SYN{$pseudo_replicate}{$breed_1}{$breed_2}/$JACK_AVG_L_SYN{$pseudo_replicate}{$breed_2}{$breed_1};
#		    $JACK_AVG_R_NONSYN{$pseudo_replicate}{$breed_1}{$breed_2}=$JACK_AVG_L_NONSYN{$pseudo_replicate}{$breed_1}{$breed_2}/$JACK_AVG_L_NONSYN{$pseudo_replicate}{$breed_2}{$breed_1};


		    #prob homozygous breed 1 not breed2
		    $JACK_RR_LOF{$pseudo_replicate}{$breed_1}{$breed_2}=$JACK_TOTAL_LL_LOF{$pseudo_replicate}{$breed_1}{$breed_2}/$JACK_TOTAL_LL_LOF{$pseudo_replicate}{$breed_2}{$breed_1};
		    $JACK_RR_SYN{$pseudo_replicate}{$breed_1}{$breed_2}=$JACK_TOTAL_LL_SYN{$pseudo_replicate}{$breed_1}{$breed_2}/$JACK_TOTAL_LL_SYN{$pseudo_replicate}{$breed_2}{$breed_1};
		    $JACK_RR_NONSYN{$pseudo_replicate}{$breed_1}{$breed_2}=$JACK_TOTAL_LL_NONSYN{$pseudo_replicate}{$breed_1}{$breed_2}/$JACK_TOTAL_LL_NONSYN{$pseudo_replicate}{$breed_2}{$breed_1};
		    
		    $WEIGHTED_JACK_RR_LOF{$pseudo_replicate}{$breed_1}{$breed_2}=$h_LOF*$TOTAL_RR_LOF{$breed_1}{$breed_2}-($h_LOF-1)*$JACK_RR_LOF{$pseudo_replicate}{$breed_1}{$breed_2};
		    $WEIGHTED_JACK_RR_SYN{$pseudo_replicate}{$breed_1}{$breed_2}=$h_SYN*$TOTAL_RR_SYN{$breed_1}{$breed_2}-($h_SYN-1)*$JACK_RR_SYN{$pseudo_replicate}{$breed_1}{$breed_2};
		    $WEIGHTED_JACK_RR_NONSYN{$pseudo_replicate}{$breed_1}{$breed_2}=$h_NONSYN*$TOTAL_RR_NONSYN{$breed_1}{$breed_2}-($h_NONSYN-1)*$JACK_RR_NONSYN{$pseudo_replicate}{$breed_1}{$breed_2};



		    
		    foreach $bin (sort numerically keys %{$JACK_TOTAL_L_CONS{$pseudo_replicate}{$breed_1}{$breed_2}}) {
			if ((exists $JACK_TOTAL_L_CONS{$pseudo_replicate}{$breed_1}{$breed_2}{$bin}) and (exists $JACK_TOTAL_L_CONS{$pseudo_replicate}{$breed_2}{$breed_1}{$bin})) { 
			    if ($BLOCK_CONS_INF_SITES{$pseudo_replicate}{$breed_1}{$breed_2}{$bin}>0) {
				$JACK_R_CONS{$pseudo_replicate}{$breed_1}{$breed_2}{$bin}=$JACK_TOTAL_L_CONS{$pseudo_replicate}{$breed_1}{$breed_2}{$bin}/$JACK_TOTAL_L_CONS{$pseudo_replicate}{$breed_2}{$breed_1}{$bin};
				print STDERR "			$JACK_R_CONS{$pseudo_replicate}{$breed_1}{$breed_2}{$bin}=$JACK_TOTAL_L_CONS{$pseudo_replicate}{$breed_1}{$breed_2}{$bin}/$JACK_TOTAL_L_CONS{$pseudo_replicate}{$breed_2}{$breed_1}{$bin};\n";
	
				#prob homozygous breed 1 not breed2		
				$JACK_RR_CONS{$pseudo_replicate}{$breed_1}{$breed_2}{$bin}=$JACK_TOTAL_LL_CONS{$pseudo_replicate}{$breed_1}{$breed_2}{$bin}/$JACK_TOTAL_LL_CONS{$pseudo_replicate}{$breed_2}{$breed_1}{$bin};


				$h_CONS="";
				$h_CONS=$TOTAL_CONS_INF_SITES{$breed_1}{$breed_2}{$bin}/$BLOCK_CONS_INF_SITES{$pseudo_replicate}{$breed_1}{$breed_2}{$bin};
				$WEIGHTED_JACK_R_CONS{$pseudo_replicate}{$breed_1}{$breed_2}{$bin}=$h_CONS*$TOTAL_R_CONS{$breed_1}{$breed_2}{$bin}-($h_CONS-1)*$JACK_R_CONS{$pseudo_replicate}{$breed_1}{$breed_2}{$bin};
				
				
				
				print STDERR "			$JACK_R_CONS{$pseudo_replicate}{$breed_1}{$breed_2}{$bin}=$JACK_TOTAL_L_CONS{$pseudo_replicate}{$breed_1}{$breed_2}{$bin}/$JACK_TOTAL_L_CONS{$pseudo_replicate}{$breed_2}{$breed_1}{$bin};\n";
#				$JACK_AVG_R_CONS{$pseudo_replicate}{$breed_1}{$breed_2}=$JACK_AVG_L_CONS{$pseudo_replicate}{$breed_1}{$breed_2}{$bin}/$JACK_AVG_L_CONS{$pseudo_replicate}{$breed_2}{$breed_1}{$bin};

				#prob homozygous breed 1 not breed2		
				$WEIGHTED_JACK_RR_CONS{$pseudo_replicate}{$breed_1}{$breed_2}{$bin}=$h_CONS*$TOTAL_RR_CONS{$breed_1}{$breed_2}{$bin}-($h_CONS-1)*$JACK_RR_CONS{$pseudo_replicate}{$breed_1}{$breed_2}{$bin};

			    }
			}
		    }
		}
	    }
	}
    }
#}


#calculate averages and weighted averages of the pseudo replicates (for which 1 block was left out each time)
foreach $pseudo_replicate (sort numerically keys %JACK_R_SYN) {
    foreach $breed_1 (sort @breeds) {
	foreach $breed_2 (sort @breeds) {
	    unless ($breed_1 eq $breed_2) {
#		unless (exists $AVG_OF_JACK_AVG_R_SYN{$breed_2}{$breed_1}) {
		    
		    
		    $AVG_OF_JACK_AVG_R_LOF{$breed_1}{$breed_2}+=$JACK_R_LOF{$pseudo_replicate}{$breed_1}{$breed_2}/$blocks;
		    $AVG_OF_JACK_AVG_R_SYN{$breed_1}{$breed_2}+=$JACK_R_SYN{$pseudo_replicate}{$breed_1}{$breed_2}/$blocks;
		    $AVG_OF_JACK_AVG_R_NONSYN{$breed_1}{$breed_2}+=$JACK_R_NONSYN{$pseudo_replicate}{$breed_1}{$breed_2}/$blocks;
		    $WEIGHTED_AVG_OF_JACK_AVG_R_LOF{$breed_1}{$breed_2}+=((1-$BLOCK_LOF_INF_SITES{$pseudo_replicate}{$breed_1}{$breed_2}/$TOTAL_LOF_INF_SITES{$breed_1}{$breed_2})*$JACK_R_LOF{$pseudo_replicate}{$breed_1}{$breed_2})/$blocks;
		    $WEIGHTED_AVG_OF_JACK_AVG_R_SYN{$breed_1}{$breed_2}+=((1-$BLOCK_SYN_INF_SITES{$pseudo_replicate}{$breed_1}{$breed_2}/$TOTAL_SYN_INF_SITES{$breed_1}{$breed_2})*$JACK_R_SYN{$pseudo_replicate}{$breed_1}{$breed_2})/$blocks;
		    $WEIGHTED_AVG_OF_JACK_AVG_R_NONSYN{$breed_1}{$breed_2}+=((1-$BLOCK_NONSYN_INF_SITES{$pseudo_replicate}{$breed_1}{$breed_2}/$TOTAL_NONSYN_INF_SITES{$breed_1}{$breed_2})*$JACK_R_NONSYN{$pseudo_replicate}{$breed_1}{$breed_2})/$blocks;
		    if (($breed_1 eq 'ckcs') && ($breed_2 eq 'beagle')) {
			print STDERR "\n$AVG_OF_JACK_AVG_R_LOF{$breed_1}{$breed_2}+=$JACK_R_LOF{$pseudo_replicate}{$breed_1}{$breed_2}/$blocks;\n";
			print "pseudo replicate: $pseudo_replicate\t$breed_1\t$breed_2\tAvg of jack avg LOF:$AVG_OF_JACK_AVG_R_LOF{$breed_1}{$breed_2}\n";
		    }


		    #prob homozygous breed 1 not breed2		
		    $AVG_OF_JACK_AVG_RR_LOF{$breed_1}{$breed_2}+=$JACK_RR_LOF{$pseudo_replicate}{$breed_1}{$breed_2}/$blocks;
		    $AVG_OF_JACK_AVG_RR_SYN{$breed_1}{$breed_2}+=$JACK_RR_SYN{$pseudo_replicate}{$breed_1}{$breed_2}/$blocks;
		    $AVG_OF_JACK_AVG_RR_NONSYN{$breed_1}{$breed_2}+=$JACK_RR_NONSYN{$pseudo_replicate}{$breed_1}{$breed_2}/$blocks;
		    $WEIGHTED_AVG_OF_JACK_AVG_RR_LOF{$breed_1}{$breed_2}+=((1-$BLOCK_LOF_INF_SITES{$pseudo_replicate}{$breed_1}{$breed_2}/$TOTAL_LOF_INF_SITES{$breed_1}{$breed_2})*$JACK_RR_LOF{$pseudo_replicate}{$breed_1}{$breed_2})/$blocks;
		    $WEIGHTED_AVG_OF_JACK_AVG_RR_SYN{$breed_1}{$breed_2}+=((1-$BLOCK_SYN_INF_SITES{$pseudo_replicate}{$breed_1}{$breed_2}/$TOTAL_SYN_INF_SITES{$breed_1}{$breed_2})*$JACK_RR_SYN{$pseudo_replicate}{$breed_1}{$breed_2})/$blocks;
		    $WEIGHTED_AVG_OF_JACK_AVG_RR_NONSYN{$breed_1}{$breed_2}+=((1-$BLOCK_NONSYN_INF_SITES{$pseudo_replicate}{$breed_1}{$breed_2}/$TOTAL_NONSYN_INF_SITES{$breed_1}{$breed_2})*$JACK_RR_NONSYN{$pseudo_replicate}{$breed_1}{$breed_2})/$blocks;



		    foreach $bin (sort numerically keys %{$TOTAL_L_CONS{$breed_1}{$breed_2}}) {
		    print "$breed_1\t$breed_2\t$bin\t";
		    $AVG_OF_JACK_AVG_R_CONS{$breed_1}{$breed_2}{$bin}+=$JACK_R_CONS{$pseudo_replicate}{$breed_1}{$breed_2}{$bin}/$blocks;
		    $WEIGHTED_AVG_OF_JACK_AVG_R_CONS{$breed_1}{$breed_2}{$bin}+=((1-$BLOCK_CONS_INF_SITES{$pseudo_replicate}{$breed_1}{$breed_2}{$bin}/$TOTAL_CONS_INF_SITES{$breed_1}{$breed_2}{$bin})*$JACK_R_CONS{$pseudo_replicate}{$breed_1}{$breed_2}{$bin})/$blocks;
#			print "KOLLA HAR!!!!!!!!\n";
#			print STDERR "$WEIGHTED_AVG_OF_JACK_AVG_R_CONS{$breed_1}{$breed_2}{$bin}+=((1-$BLOCK_CONS_INF_SITES{$pseudo_replicate}{$breed_1}{$breed_2}{$bin}/$TOTAL_CONS_INF_SITES{$breed_1}{$breed_2}{$bin})*$JACK_R_CONS{$pseudo_replicate}{$breed_1}{$breed_2}{$bin})/$blocks\n";
		    
#		    print STDERR "		    $AVG_OF_JACK_AVG_R_CONS{$breed_1}{$breed_2}{$bin}+=$JACK_AVG_R_CONS{$pseudo_replicate}{$breed_1}{$breed_2}{$bin}/$blocks;\n";

		    #prob homozygous breed 1 not breed2		
		    $AVG_OF_JACK_AVG_RR_CONS{$breed_1}{$breed_2}{$bin}+=$JACK_RR_CONS{$pseudo_replicate}{$breed_1}{$breed_2}{$bin}/$blocks;
		    $WEIGHTED_AVG_OF_JACK_AVG_RR_CONS{$breed_1}{$breed_2}{$bin}+=((1-$BLOCK_CONS_INF_SITES{$pseudo_replicate}{$breed_1}{$breed_2}{$bin}/$TOTAL_CONS_INF_SITES{$breed_1}{$breed_2}{$bin})*$JACK_RR_CONS{$pseudo_replicate}{$breed_1}{$breed_2}{$bin})/$blocks;

		    }
		}
	    }
	}
    }
#}


#unweigthed and weighted bias corrected jackknife estimator of R
foreach $breed_1 (sort @breeds) {
    foreach $breed_2 (sort @breeds) {
	unless ($breed_1 eq $breed_2) {
#	    unless (exists $BIAS_CORRECTED_R_SYN{$breed_2}{$breed_1}) {
		$BIAS_CORRECTED_R_LOF{$breed_1}{$breed_2}=$blocks*$TOTAL_R_LOF{$breed_1}{$breed_2}-($blocks-1)*$AVG_OF_JACK_AVG_R_LOF{$breed_1}{$breed_2};

print STDERR "		$BIAS_CORRECTED_R_LOF{$breed_1}{$breed_2}=$blocks*$TOTAL_R_LOF{$breed_1}{$breed_2}-($blocks-1)*$AVG_OF_JACK_AVG_R_LOF{$breed_1}{$breed_2};
\n";		$BIAS_CORRECTED_R_SYN{$breed_1}{$breed_2}=$blocks*$TOTAL_R_SYN{$breed_1}{$breed_2}-($blocks-1)*$AVG_OF_JACK_AVG_R_SYN{$breed_1}{$breed_2};
		$BIAS_CORRECTED_R_NONSYN{$breed_1}{$breed_2}=$blocks*$TOTAL_R_NONSYN{$breed_1}{$breed_2}-($blocks-1)*$AVG_OF_JACK_AVG_R_NONSYN{$breed_1}{$breed_2};
		
		$WEIGHTED_BIAS_CORRECTED_R_LOF{$breed_1}{$breed_2}=$blocks*$TOTAL_R_LOF{$breed_1}{$breed_2}-($blocks*$WEIGHTED_AVG_OF_JACK_AVG_R_LOF{$breed_1}{$breed_2});
		$WEIGHTED_BIAS_CORRECTED_R_SYN{$breed_1}{$breed_2}=$blocks*$TOTAL_R_SYN{$breed_1}{$breed_2}-($blocks*$WEIGHTED_AVG_OF_JACK_AVG_R_SYN{$breed_1}{$breed_2});
		$WEIGHTED_BIAS_CORRECTED_R_NONSYN{$breed_1}{$breed_2}=$blocks*$TOTAL_R_NONSYN{$breed_1}{$breed_2}-($blocks*$WEIGHTED_AVG_OF_JACK_AVG_R_NONSYN{$breed_1}{$breed_2});


		#prob homozygous breed 1 not breed2		
		$BIAS_CORRECTED_RR_LOF{$breed_1}{$breed_2}=$blocks*$TOTAL_RR_LOF{$breed_1}{$breed_2}-($blocks-1)*$AVG_OF_JACK_AVG_RR_LOF{$breed_1}{$breed_2};
		$BIAS_CORRECTED_RR_SYN{$breed_1}{$breed_2}=$blocks*$TOTAL_RR_SYN{$breed_1}{$breed_2}-($blocks-1)*$AVG_OF_JACK_AVG_RR_SYN{$breed_1}{$breed_2};
		$BIAS_CORRECTED_RR_NONSYN{$breed_1}{$breed_2}=$blocks*$TOTAL_RR_NONSYN{$breed_1}{$breed_2}-($blocks-1)*$AVG_OF_JACK_AVG_RR_NONSYN{$breed_1}{$breed_2};
		
		$WEIGHTED_BIAS_CORRECTED_RR_LOF{$breed_1}{$breed_2}=$blocks*$TOTAL_RR_LOF{$breed_1}{$breed_2}-($blocks*$WEIGHTED_AVG_OF_JACK_AVG_RR_LOF{$breed_1}{$breed_2});
		$WEIGHTED_BIAS_CORRECTED_RR_SYN{$breed_1}{$breed_2}=$blocks*$TOTAL_RR_SYN{$breed_1}{$breed_2}-($blocks*$WEIGHTED_AVG_OF_JACK_AVG_RR_SYN{$breed_1}{$breed_2});
		$WEIGHTED_BIAS_CORRECTED_RR_NONSYN{$breed_1}{$breed_2}=$blocks*$TOTAL_RR_NONSYN{$breed_1}{$breed_2}-($blocks*$WEIGHTED_AVG_OF_JACK_AVG_RR_NONSYN{$breed_1}{$breed_2});


		
		if (($breed_1 eq 'ckcs') && ($breed_2 eq 'beagle')) {
		    print "bias corrected jackknife estimator LOF:$breed_1\t$breed_2\t$BIAS_CORRECTED_R_LOF{$breed_1}{$breed_2}\n"; 
		}
		foreach $bin (sort numerically keys %{$TOTAL_L_CONS{$breed_1}{$breed_2}}) {
		    $BIAS_CORRECTED_R_CONS{$breed_1}{$breed_2}{$bin}=$blocks*$TOTAL_R_CONS{$breed_1}{$breed_2}{$bin}-($blocks-1)*$AVG_OF_JACK_AVG_R_CONS{$breed_1}{$breed_2}{$bin};
		    $WEIGHTED_BIAS_CORRECTED_R_CONS{$breed_1}{$breed_2}{$bin}=$blocks*$TOTAL_R_CONS{$breed_1}{$breed_2}{$bin}-($blocks*$WEIGHTED_AVG_OF_JACK_AVG_R_CONS{$breed_1}{$breed_2}{$bin});
#		    print "$breed_1\t$breed_2\t$bin\n";
#		    print STDERR "$BIAS_CORRECTED_R_CONS{$breed_1}{$breed_2}{$bin}=$blocks*$TOTAL_R_CONS{$breed_1}{$breed_2}{$bin}-($blocks-1)*$AVG_OF_JACK_AVG_R_CONS{$breed_1}{$breed_2}{$bin};\n";
#		    print STDERR "$WEIGHTED_BIAS_CORRECTED_R_CONS{$breed_1}{$breed_2}{$bin}=$blocks*$TOTAL_R_CONS{$breed_1}{$breed_2}{$bin}-($blocks*$WEIGHTED_AVG_OF_JACK_AVG_R_CONS{$breed_1}{$breed_2}{$bin});\n";
		   

		#prob homozygous breed 1 not breed2		
 		    $BIAS_CORRECTED_RR_CONS{$breed_1}{$breed_2}{$bin}=$blocks*$TOTAL_RR_CONS{$breed_1}{$breed_2}{$bin}-($blocks-1)*$AVG_OF_JACK_AVG_RR_CONS{$breed_1}{$breed_2}{$bin};
		    $WEIGHTED_BIAS_CORRECTED_RR_CONS{$breed_1}{$breed_2}{$bin}=$blocks*$TOTAL_RR_CONS{$breed_1}{$breed_2}{$bin}-($blocks*$WEIGHTED_AVG_OF_JACK_AVG_RR_CONS{$breed_1}{$breed_2}{$bin});

		}	    
	    }
	}
    }
#}



foreach $breed_1 (sort @breeds) {
    foreach $breed_2 (sort @breeds) {
	unless ($breed_1 eq $breed_2) {
#	    if (($breed_1 eq 'ckcs') && ($breed_2 eq 'beagle')) {
	    if (exists $TOTAL_R_LOF{$breed_1}{$breed_2}) {
		print "$breed_1\t$breed_2:\n";
		print "LOF\ttotal R\t$TOTAL_R_LOF{$breed_1}{$breed_2}\tbias corrected jackknife estimate\t$BIAS_CORRECTED_R_LOF{$breed_1}{$breed_2}\twieghted bias corrected jackknife estimate\t$WEIGHTED_BIAS_CORRECTED_R_LOF{$breed_1}{$breed_2}\tINF_SITES\t$TOTAL_LOF_INF_SITES{$breed_1}{$breed_2}\n";
		print "SYN\ttotal R\t$TOTAL_R_SYN{$breed_1}{$breed_2}\tbias corrected jackknife estimate\t$BIAS_CORRECTED_R_SYN{$breed_1}{$breed_2}\twieghted bias corrected jackknife estimate\t$WEIGHTED_BIAS_CORRECTED_R_SYN{$breed_1}{$breed_2}\tINF_SITES\t$TOTAL_SYN_INF_SITES{$breed_1}{$breed_2}\n";
		print "NONSYN\ttotal R\t$TOTAL_R_NONSYN{$breed_1}{$breed_2}\tbias corrected jackknife estimate\t$BIAS_CORRECTED_R_NONSYN{$breed_1}{$breed_2}\t\twieghted bias corrected jackknife estimate\t$WEIGHTED_BIAS_CORRECTED_R_NONSYN{$breed_1}{$breed_2}INF_SITES\t$TOTAL_NONSYN_INF_SITES{$breed_1}{$breed_2}\n";
		foreach $bin (sort numerically keys %{$AVG_OF_JACK_AVG_R_CONS{$breed_1}{$breed_2}}) {
		    print "CONS\t$bin\ttotal R\t$TOTAL_R_CONS{$breed_1}{$breed_2}{$bin}\taverage of all pseudo replicates\t$BIAS_CORRECTED_R_CONS{$breed_1}{$breed_2}{$bin}\twighted bias corrected jackknife estimate\t$WEIGHTED_BIAS_CORRECTED_R_CONS{$breed_1}{$breed_2}{$bin}\tINF_SITES\t$TOTAL_CONS_INF_SITES{$breed_1}{$breed_2}{$bin}\n";
		    
		}
	    }
	}
    }
}



#jackknife variance estimation 
foreach $block (sort numerically keys %JACK_R_SYN) {
    foreach $breed_1 (sort @breeds) {
	foreach $breed_2 (sort @breeds) {
	    unless ($breed_1 eq $breed_2) {
#		unless (exists $var_R_SYN{$breed_2}{$breed_1}) {

		$h_LOF=$h_SYN=$h_NONSYN="";
#		print "$block\t$breed_1\t$breed_2\n";
#		print "LOF sites in excluded block: $BLOCK_LOF_INF_SITES{$block}{$breed_1}{$breed_2}\n"; 
#		print "total LOF sites in all genome: $TOTAL_LOF_INF_SITES{$breed_1}{$breed_2}\n";
		$var_R_LOF{$breed_1}{$breed_2}+=((($JACK_R_LOF{$block}{$breed_1}{$breed_2}-$AVG_OF_JACK_AVG_R_LOF{$breed_1}{$breed_2})**2)*(($blocks-1)/$blocks));
#		print STDERR "		$var_R_LOF{$breed_1}{$breed_2}+=((($JACK_AVG_R_LOF{$block}{$breed_1}{$breed_2}-$AVG_OF_JACK_AVG_R_LOF{$breed_1}{$breed_2})**2)*(($blocks-1)/$blocks));\n";
		$var_R_SYN{$breed_1}{$breed_2}+=((($JACK_R_SYN{$block}{$breed_1}{$breed_2}-$AVG_OF_JACK_AVG_R_SYN{$breed_1}{$breed_2})**2)*(($blocks-1)/$blocks));
		$var_R_NONSYN{$breed_1}{$breed_2}+=((($JACK_R_NONSYN{$block}{$breed_1}{$breed_2}-$AVG_OF_JACK_AVG_R_NONSYN{$breed_1}{$breed_2})**2)*(($blocks-1)/$blocks));


		$h_LOF=$TOTAL_LOF_INF_SITES{$breed_1}{$breed_2}/$BLOCK_LOF_INF_SITES{$block}{$breed_1}{$breed_2};
		$h_SYN=$TOTAL_SYN_INF_SITES{$breed_1}{$breed_2}/$BLOCK_SYN_INF_SITES{$block}{$breed_1}{$breed_2};
		$h_NONSYN=$TOTAL_NONSYN_INF_SITES{$breed_1}{$breed_2}/$BLOCK_NONSYN_INF_SITES{$block}{$breed_1}{$breed_2};

		$weighted_var_R_LOF{$breed_1}{$breed_2}+=((1/($h_LOF-1))*($WEIGHTED_JACK_R_LOF{$block}{$breed_1}{$breed_2}-$WEIGHTED_BIAS_CORRECTED_R_LOF{$breed_1}{$breed_2})**2)/$blocks;
		$weighted_var_R_SYN{$breed_1}{$breed_2}+=((1/($h_SYN-1))*($WEIGHTED_JACK_R_SYN{$block}{$breed_1}{$breed_2}-$WEIGHTED_BIAS_CORRECTED_R_SYN{$breed_1}{$breed_2})**2)/$blocks;
		$weighted_var_R_NONSYN{$breed_1}{$breed_2}+=((1/($h_NONSYN-1))*($WEIGHTED_JACK_R_NONSYN{$block}{$breed_1}{$breed_2}-$WEIGHTED_BIAS_CORRECTED_R_NONSYN{$breed_1}{$breed_2})**2)/$blocks;
		##################################################


		#prob homozygous breed 1 not breed2		
		$var_RR_LOF{$breed_1}{$breed_2}+=((($JACK_RR_LOF{$block}{$breed_1}{$breed_2}-$AVG_OF_JACK_AVG_RR_LOF{$breed_1}{$breed_2})**2)*(($blocks-1)/$blocks));
		$var_RR_SYN{$breed_1}{$breed_2}+=((($JACK_RR_SYN{$block}{$breed_1}{$breed_2}-$AVG_OF_JACK_AVG_RR_SYN{$breed_1}{$breed_2})**2)*(($blocks-1)/$blocks));
		$var_RR_NONSYN{$breed_1}{$breed_2}+=((($JACK_RR_NONSYN{$block}{$breed_1}{$breed_2}-$AVG_OF_JACK_AVG_RR_NONSYN{$breed_1}{$breed_2})**2)*(($blocks-1)/$blocks));

		$weighted_var_RR_LOF{$breed_1}{$breed_2}+=((1/($h_LOF-1))*($WEIGHTED_JACK_RR_LOF{$block}{$breed_1}{$breed_2}-$WEIGHTED_BIAS_CORRECTED_RR_LOF{$breed_1}{$breed_2})**2)/$blocks;
		$weighted_var_RR_SYN{$breed_1}{$breed_2}+=((1/($h_SYN-1))*($WEIGHTED_JACK_RR_SYN{$block}{$breed_1}{$breed_2}-$WEIGHTED_BIAS_CORRECTED_RR_SYN{$breed_1}{$breed_2})**2)/$blocks;
		$weighted_var_RR_NONSYN{$breed_1}{$breed_2}+=((1/($h_NONSYN-1))*($WEIGHTED_JACK_RR_NONSYN{$block}{$breed_1}{$breed_2}-$WEIGHTED_BIAS_CORRECTED_RR_NONSYN{$breed_1}{$breed_2})**2)/$blocks;



		foreach $bin (sort numerically keys %{$L_CONS{$block}{$breed_1}{$breed_2}}) {
		    print "$block\t$breed_1\t$breed_2\n";
		    $h_CONS="";
print STDERR "	    $h_CONS=$TOTAL_CONS_INF_SITES{$breed_1}{$breed_2}{$bin}/$BLOCK_CONS_INF_SITES{$block}{$breed_1}{$breed_2}{$bin};\n";

		    $h_CONS=$TOTAL_CONS_INF_SITES{$breed_1}{$breed_2}{$bin}/$BLOCK_CONS_INF_SITES{$block}{$breed_1}{$breed_2}{$bin};
		    if ($h_CONS>1) {
			$weighted_var_R_CONS{$breed_1}{$breed_2}{$bin}+=((1/($h_CONS-1))*($WEIGHTED_JACK_R_CONS{$block}{$breed_1}{$breed_2}{$bin}-$WEIGHTED_BIAS_CORRECTED_R_CONS{$breed_1}{$breed_2}{$bin})**2)/$blocks;
			
			$VAR_R_CONS{$breed_1}{$breed_2}{$bin}+=((($JACK_R_CONS{$block}{$breed_1}{$breed_2}{$bin}-$AVG_OF_JACK_AVG_R_CONS{$breed_1}{$breed_2}{$bin})**2)*(($blocks-1)/$blocks));

			#prob homozygous breed 1 not breed2		
			$weighted_var_RR_CONS{$breed_1}{$breed_2}{$bin}+=((1/($h_CONS-1))*($WEIGHTED_JACK_RR_CONS{$block}{$breed_1}{$breed_2}{$bin}-$WEIGHTED_BIAS_CORRECTED_RR_CONS{$breed_1}{$breed_2}{$bin})**2)/$blocks;
			
			$VAR_RR_CONS{$breed_1}{$breed_2}{$bin}+=((($JACK_RR_CONS{$block}{$breed_1}{$breed_2}{$bin}-$AVG_OF_JACK_AVG_RR_CONS{$breed_1}{$breed_2}{$bin})**2)*(($blocks-1)/$blocks));

		    }
		    
		}
	    }
	}
    }
}
#}

open OUT, ">>$results_file";
print OUT "Mutation\tbreed_1\tinformative_sites_breed_1\tinformative_sites_breed_2\tbreed_2\tL_breed1_breed2\tL_breed2_breed1\tR\tbias corrected jackknife estimator of R\tweighted bias corrected jackknife estimator of R\tvariance\tSD\tmin_conf\tmax_conf\tweighted variance\tweighted SD\tweighted min_conf\tweighted max_conf\tjackknife blocks containing bin snps\n";
foreach $breed_1 (sort keys %var_R_SYN) {
    foreach $breed_2 (sort keys %{$var_R_SYN{$breed_1}}) {
	    $SD_R_LOF=sqrt($var_R_LOF{$breed_1}{$breed_2});
	    $min_conf_R_LOF=$TOTAL_R_LOF{$breed_1}{$breed_2}-1.96*$SD_R_LOF;
	    $max_conf_R_LOF=$TOTAL_R_LOF{$breed_1}{$breed_2}+1.96*$SD_R_LOF;

	    $weighted_SD_R_LOF=sqrt($weighted_var_R_LOF{$breed_1}{$breed_2});
	    $weighted_min_conf_R_LOF=$WEIGHTED_BIAS_CORRECTED_R_LOF{$breed_1}{$breed_2}-1.96*$weighted_SD_R_LOF;
	    $weighted_max_conf_R_LOF=$WEIGHTED_BIAS_CORRECTED_R_LOF{$breed_1}{$breed_2}+1.96*$weighted_SD_R_LOF;

	    $SD_R_SYN=sqrt($var_R_SYN{$breed_1}{$breed_2});
	    $min_conf_R_SYN=$TOTAL_R_SYN{$breed_1}{$breed_2}-1.96*$SD_R_SYN;
	    $max_conf_R_SYN=$TOTAL_R_SYN{$breed_1}{$breed_2}+1.96*$SD_R_SYN;

	    $weighted_SD_R_SYN=sqrt($var_R_SYN{$breed_1}{$breed_2});
	    $weighted_min_conf_R_SYN=$WEIGHTED_BIAS_CORRECTED_R_SYN{$breed_1}{$breed_2}-1.96*$weighted_SD_R_SYN;
	    $weighted_max_conf_R_SYN=$WEIGHTED_BIAS_CORRECTED_R_SYN{$breed_1}{$breed_2}+1.96*$weighted_SD_R_SYN;

	    $SD_R_NONSYN=sqrt($var_R_NONSYN{$breed_1}{$breed_2});
	    $min_conf_R_NONSYN=$TOTAL_R_NONSYN{$breed_1}{$breed_2}-1.96*$SD_R_NONSYN;
	    $max_conf_R_NONSYN=$TOTAL_R_NONSYN{$breed_1}{$breed_2}+1.96*$SD_R_NONSYN;

	    $weighted_SD_R_NONSYN=sqrt($var_R_NONSYN{$breed_1}{$breed_2});
	    $weighted_min_conf_R_NONSYN=$WEIGHTED_BIAS_CORRECTED_R_NONSYN{$breed_1}{$breed_2}-1.96*$weighted_SD_R_NONSYN;
	    $weighted_max_conf_R_NONSYN=$WEIGHTED_BIAS_CORRECTED_R_NONSYN{$breed_1}{$breed_2}+1.96*$weighted_SD_R_NONSYN;

	    print OUT "LOF\t$breed_1\t$breed_2\t$TOTAL_LOF_INF_SITES{$breed_1}{$breed_2}\t$TOTAL_LOF_INF_SITES{$breed_2}{$breed_1}\t$TOTAL_L_LOF{$breed_1}{$breed_2}\t$TOTAL_L_LOF{$breed_2}{$breed_1}\t$TOTAL_R_LOF{$breed_1}{$breed_2}\t$BIAS_CORRECTED_R_LOF{$breed_1}{$breed_2}\t$WEIGHTED_BIAS_CORRECTED_R_LOF{$breed_1}{$breed_2}\t$var_R_LOF{$breed_1}{$breed_2}\t$SD_R_LOF\t$min_conf_R_LOF\t$max_conf_R_LOF\t$weighted_var_R_LOF{$breed_1}{$breed_2}\t$weighted_SD_R_LOF\t$weighted_min_conf_R_LOF\t$weighted_max_conf_R_LOF\tNA\n";
	    print OUT "SYN\t$breed_1\t$breed_2\t$TOTAL_SYN_INF_SITES{$breed_1}{$breed_2}\t$TOTAL_SYN_INF_SITES{$breed_2}{$breed_1}\t$TOTAL_L_SYN{$breed_1}{$breed_2}\t$TOTAL_L_SYN{$breed_2}{$breed_1}\t$TOTAL_R_SYN{$breed_1}{$breed_2}\t$BIAS_CORRECTED_R_SYN{$breed_1}{$breed_2}\t$WEIGHTED_BIAS_CORRECTED_R_SYN{$breed_1}{$breed_2}\t$var_R_SYN{$breed_1}{$breed_2}\t$SD_R_SYN\t$min_conf_R_SYN\t$max_conf_R_SYN\t$weighted_var_R_SYN{$breed_1}{$breed_2}\t$weighted_SD_R_SYN\t$weighted_min_conf_R_SYN\t$weighted_max_conf_R_SYN\tNA\n";
	    print OUT "NONSYN\t$breed_1\t$breed_2\t$TOTAL_NONSYN_INF_SITES{$breed_1}{$breed_2}\t$TOTAL_NONSYN_INF_SITES{$breed_2}{$breed_1}\t$TOTAL_L_NONSYN{$breed_1}{$breed_2}\t$TOTAL_L_NONSYN{$breed_2}{$breed_1}\t$TOTAL_R_NONSYN{$breed_1}{$breed_2}\t$BIAS_CORRECTED_R_NONSYN{$breed_1}{$breed_2}\t$WEIGHTED_BIAS_CORRECTED_R_NONSYN{$breed_1}{$breed_2}\t$var_R_NONSYN{$breed_1}{$breed_2}\t$SD_R_NONSYN\t$min_conf_R_NONSYN\t$max_conf_R_NONSYN\t$weighted_var_R_NONSYN{$breed_1}{$breed_2}\t$weighted_SD_R_NONSYN\t$weighted_min_conf_R_NONSYN\t$weighted_max_conf_R_NONSYN\tNA\n";


	foreach $bin (sort numerically keys %{$VAR_R_CONS{$breed_1}{$breed_2}}) {

	    print "$breed_1\t$breed_2\t$bin\n";
	    if (exists  $VAR_R_CONS{$breed_1}{$breed_2}{$bin}) {
	    $SD_R_CONS=sqrt($VAR_R_CONS{$breed_1}{$breed_2}{$bin});
	    $min_conf_R_CONS=$TOTAL_R_CONS{$breed_1}{$breed_2}{$bin}-1.96*$SD_R_CONS;
	    $max_conf_R_CONS=$TOTAL_R_CONS{$breed_1}{$breed_2}{$bin}+1.96*$SD_R_CONS;

	    $weighted_SD_R_CONS=sqrt($weighted_var_R_CONS{$breed_1}{$breed_2}{$bin});
	    $weighted_min_conf_R_CONS=$WEIGHTED_BIAS_CORRECTED_R_CONS{$breed_1}{$breed_2}{$bin}-1.96*$weighted_SD_R_CONS;
	    $weighted_max_conf_R_CONS=$WEIGHTED_BIAS_CORRECTED_R_CONS{$breed_1}{$breed_2}{$bin}+1.96*$weighted_SD_R_CONS;
	    print OUT "$bin\t$breed_1\t$breed_2\t$TOTAL_CONS_INF_SITES{$breed_1}{$breed_2}{$bin}\t$TOTAL_CONS_INF_SITES{$breed_2}{$breed_1}{$bin}\t$TOTAL_L_CONS{$breed_1}{$breed_2}{$bin}\t$TOTAL_L_CONS{$breed_2}{$breed_1}{$bin}\t$TOTAL_R_CONS{$breed_1}{$breed_2}{$bin}\t$BIAS_CORRECTED_R_CONS{$breed_1}{$breed_2}{$bin}\t$WEIGHTED_BIAS_CORRECTED_R_CONS{$breed_1}{$breed_2}{$bin}\t$VAR_R_CONS{$breed_1}{$breed_2}{$bin}\t$SD_R_CONS\t$min_conf_R_CONS\t$max_conf_R_CONS\t$weighted_var_R_CONS{$breed_1}{$breed_2}{$bin}\t$weighted_SD_R_CONS\t$weighted_min_conf_R_CONS\t$weighted_max_conf_R_CONS\t$COUNT_AVG_R_CONS_NONZERO_BLOCKS{$breed_1}{$breed_2}{$bin}\n";

	    }
	}
    }
}
close OUT;



#Print results of derived homozygos analyses
open OUT, ">>$rresults_file";
print OUT "Mutation\tbreed_1\tinformative_sites_breed_1\tinformative_sites_breed_2\tbreed_2\tLL_breed1_breed2\tLL_breed2_breed1\tRR\tbias corrected jackknife estimator of RR\tweighted bias corrected jackknife estimator of RR\tvariance\tSD\tmin_conf\tmax_conf\tweighted variance\tweighted SD\tweighted min_conf\tweighted max_conf\tjackknife blocks containing bin snps\n";
foreach $breed_1 (sort keys %var_RR_SYN) {
    foreach $breed_2 (sort keys %{$var_RR_SYN{$breed_1}}) {
	    $SD_RR_LOF=sqrt($var_RR_LOF{$breed_1}{$breed_2});
	    $min_conf_RR_LOF=$TOTAL_RR_LOF{$breed_1}{$breed_2}-1.96*$SD_RR_LOF;
	    $max_conf_RR_LOF=$TOTAL_RR_LOF{$breed_1}{$breed_2}+1.96*$SD_RR_LOF;

	    $weighted_SD_RR_LOF=sqrt($weighted_var_RR_LOF{$breed_1}{$breed_2});
	    $weighted_min_conf_RR_LOF=$WEIGHTED_BIAS_CORRECTED_RR_LOF{$breed_1}{$breed_2}-1.96*$weighted_SD_RR_LOF;
	    $weighted_max_conf_RR_LOF=$WEIGHTED_BIAS_CORRECTED_RR_LOF{$breed_1}{$breed_2}+1.96*$weighted_SD_RR_LOF;

	    $SD_RR_SYN=sqrt($var_RR_SYN{$breed_1}{$breed_2});
	    $min_conf_RR_SYN=$TOTAL_RR_SYN{$breed_1}{$breed_2}-1.96*$SD_RR_SYN;
	    $max_conf_RR_SYN=$TOTAL_RR_SYN{$breed_1}{$breed_2}+1.96*$SD_RR_SYN;

	    $weighted_SD_RR_SYN=sqrt($var_RR_SYN{$breed_1}{$breed_2});
	    $weighted_min_conf_RR_SYN=$WEIGHTED_BIAS_CORRECTED_RR_SYN{$breed_1}{$breed_2}-1.96*$weighted_SD_RR_SYN;
	    $weighted_max_conf_RR_SYN=$WEIGHTED_BIAS_CORRECTED_RR_SYN{$breed_1}{$breed_2}+1.96*$weighted_SD_RR_SYN;

	    $SD_RR_NONSYN=sqrt($var_RR_NONSYN{$breed_1}{$breed_2});
	    $min_conf_RR_NONSYN=$TOTAL_RR_NONSYN{$breed_1}{$breed_2}-1.96*$SD_RR_NONSYN;
	    $max_conf_RR_NONSYN=$TOTAL_RR_NONSYN{$breed_1}{$breed_2}+1.96*$SD_RR_NONSYN;

	    $weighted_SD_RR_NONSYN=sqrt($var_RR_NONSYN{$breed_1}{$breed_2});
	    $weighted_min_conf_RR_NONSYN=$WEIGHTED_BIAS_CORRECTED_RR_NONSYN{$breed_1}{$breed_2}-1.96*$weighted_SD_RR_NONSYN;
	    $weighted_max_conf_RR_NONSYN=$WEIGHTED_BIAS_CORRECTED_RR_NONSYN{$breed_1}{$breed_2}+1.96*$weighted_SD_RR_NONSYN;

	    print OUT "LOF\t$breed_1\t$breed_2\t$TOTAL_LOF_INF_SITES{$breed_1}{$breed_2}\t$TOTAL_LOF_INF_SITES{$breed_2}{$breed_1}\t$TOTAL_LL_LOF{$breed_1}{$breed_2}\t$TOTAL_LL_LOF{$breed_2}{$breed_1}\t$TOTAL_RR_LOF{$breed_1}{$breed_2}\t$BIAS_CORRECTED_RR_LOF{$breed_1}{$breed_2}\t$WEIGHTED_BIAS_CORRECTED_RR_LOF{$breed_1}{$breed_2}\t$var_RR_LOF{$breed_1}{$breed_2}\t$SD_RR_LOF\t$min_conf_RR_LOF\t$max_conf_RR_LOF\t$weighted_var_RR_LOF{$breed_1}{$breed_2}\t$weighted_SD_RR_LOF\t$weighted_min_conf_RR_LOF\t$weighted_max_conf_RR_LOF\tNA\n";
	    print OUT "SYN\t$breed_1\t$breed_2\t$TOTAL_SYN_INF_SITES{$breed_1}{$breed_2}\t$TOTAL_SYN_INF_SITES{$breed_2}{$breed_1}\t$TOTAL_LL_SYN{$breed_1}{$breed_2}\t$TOTAL_LL_SYN{$breed_2}{$breed_1}\t$TOTAL_RR_SYN{$breed_1}{$breed_2}\t$BIAS_CORRECTED_RR_SYN{$breed_1}{$breed_2}\t$WEIGHTED_BIAS_CORRECTED_RR_SYN{$breed_1}{$breed_2}\t$var_RR_SYN{$breed_1}{$breed_2}\t$SD_RR_SYN\t$min_conf_RR_SYN\t$max_conf_RR_SYN\t$weighted_var_RR_SYN{$breed_1}{$breed_2}\t$weighted_SD_RR_SYN\t$weighted_min_conf_RR_SYN\t$weighted_max_conf_RR_SYN\tNA\n";
	    print OUT "NONSYN\t$breed_1\t$breed_2\t$TOTAL_NONSYN_INF_SITES{$breed_1}{$breed_2}\t$TOTAL_NONSYN_INF_SITES{$breed_2}{$breed_1}\t$TOTAL_LL_NONSYN{$breed_1}{$breed_2}\t$TOTAL_LL_NONSYN{$breed_2}{$breed_1}\t$TOTAL_RR_NONSYN{$breed_1}{$breed_2}\t$BIAS_CORRECTED_RR_NONSYN{$breed_1}{$breed_2}\t$WEIGHTED_BIAS_CORRECTED_RR_NONSYN{$breed_1}{$breed_2}\t$var_RR_NONSYN{$breed_1}{$breed_2}\t$SD_RR_NONSYN\t$min_conf_RR_NONSYN\t$max_conf_RR_NONSYN\t$weighted_var_RR_NONSYN{$breed_1}{$breed_2}\t$weighted_SD_RR_NONSYN\t$weighted_min_conf_RR_NONSYN\t$weighted_max_conf_RR_NONSYN\tNA\n";


	foreach $bin (sort numerically keys %{$VAR_RR_CONS{$breed_1}{$breed_2}}) {

	    print "$breed_1\t$breed_2\t$bin\n";
	    if (exists  $VAR_RR_CONS{$breed_1}{$breed_2}{$bin}) {
	    $SD_RR_CONS=sqrt($VAR_RR_CONS{$breed_1}{$breed_2}{$bin});
	    $min_conf_RR_CONS=$TOTAL_RR_CONS{$breed_1}{$breed_2}{$bin}-1.96*$SD_RR_CONS;
	    $max_conf_RR_CONS=$TOTAL_RR_CONS{$breed_1}{$breed_2}{$bin}+1.96*$SD_RR_CONS;

	    $weighted_SD_RR_CONS=sqrt($weighted_var_RR_CONS{$breed_1}{$breed_2}{$bin});
	    $weighted_min_conf_RR_CONS=$WEIGHTED_BIAS_CORRECTED_RR_CONS{$breed_1}{$breed_2}{$bin}-1.96*$weighted_SD_RR_CONS;
	    $weighted_max_conf_RR_CONS=$WEIGHTED_BIAS_CORRECTED_RR_CONS{$breed_1}{$breed_2}{$bin}+1.96*$weighted_SD_RR_CONS;
	    print OUT "$bin\t$breed_1\t$breed_2\t$TOTAL_CONS_INF_SITES{$breed_1}{$breed_2}{$bin}\t$TOTAL_CONS_INF_SITES{$breed_2}{$breed_1}{$bin}\t$TOTAL_LL_CONS{$breed_1}{$breed_2}{$bin}\t$TOTAL_LL_CONS{$breed_2}{$breed_1}{$bin}\t$TOTAL_RR_CONS{$breed_1}{$breed_2}{$bin}\t$BIAS_CORRECTED_RR_CONS{$breed_1}{$breed_2}{$bin}\t$WEIGHTED_BIAS_CORRECTED_RR_CONS{$breed_1}{$breed_2}{$bin}\t$VAR_RR_CONS{$breed_1}{$breed_2}{$bin}\t$SD_RR_CONS\t$min_conf_RR_CONS\t$max_conf_RR_CONS\t$weighted_var_RR_CONS{$breed_1}{$breed_2}{$bin}\t$weighted_SD_RR_CONS\t$weighted_min_conf_RR_CONS\t$weighted_max_conf_RR_CONS\t$COUNT_AVG_RR_CONS_NONZERO_BLOCKS{$breed_1}{$breed_2}{$bin}\n";

	    }
	}
    }
}
close OUT;



sub cons_bin {
    $cons_bin='';
    $cons=$_[0];
#    if (($cons>0) && ($cons<=1)){
#	$cons_bin=0;
#    }
#    if (($cons>1) && ($cons<=2)) {
#	$cons_bin=1;
#    }
#    if (($cons>2) && ($cons<=3)) {
#	$cons_bin=2;
#    }
#    if (($cons>3) && ($cons<=4)) {
#	$cons_bin=3;
#    }
#    if (($cons>4) && ($cons<=5)) {
#	$cons_bin=4;
#    }
    if (($cons>2) && ($cons<=5)) {
	$cons_bin=2;
    }
#    if (($cons>5) && ($cons<=6)) {
#	$cons_bin=5;
#    }
#    if (($cons>6) && ($cons<=7)) {
#	$cons_bin=6;
#    }
    if ($cons>5) {
	$cons_bin=5;
    }
#    if ($cons>7) {
#	$cons_bin=7;
#   }
#    if (($cons>8) && ($cons<20)) {
#	$cons_bin=8;                   
#    }
#    if (($cons>8) && ($cons<=9)) {
#	$cons_bin=9;
#    }
#    if (($cons>9) && ($cons<=10)) {
#	$cons_bin=10;
#    }
#    if ($cons>10) {
#	$cons_bin=11;
#   }
    return $cons_bin;
}


 

sub numerically {
  $a<=>$b;
}







      ### Distribution of avg. pairwise FST per breed ###
      $FST_BIN{'beagle'}{bin($beagle_fst)}++;
      $FST_BIN{'ckcs'}{bin($ckcs_fst)}++;
      $FST_BIN{'gs'}{bin($gs_fst)}++;
      $FST_BIN{'gr'}{bin($gr_fst)}++;
      $FST_BIN{'lr'}{bin($lr_fst)}++;
      $FST_BIN{'poodle'}{bin($poodle_fst)}++;
      $FST_BIN{'rotw'}{bin($rotw_fst)}++;
      $FST_BIN{'whwt'}{bin($whwt_fst)}++;
      


      sub bin {
	$freq=$_[0];
	if ($freq<=0.05) {
	  $bin=0.05;
	}
	if (($freq>0.05) && ($freq<=0.1)) {
	  $bin=0.1;
	}
	if (($freq>0.1) && ($freq<=0.15)) {
	  $bin=0.15;
	}
	if (($freq>0.15) && ($freq<=0.2)) {
	  $bin=0.2;
	}
	if (($freq>0.2) && ($freq<=0.25)) {
	  $bin=0.25;
	}
	if (($freq>0.25) && ($freq<=0.3)) {
	  $bin=0.3;
	}
	if (($freq>0.3) && ($freq<=0.35)) {
	  $bin=0.35;
	}
	if (($freq>0.35) && ($freq<=0.4)) {
	  $bin=0.4;
	}
	if (($freq>0.4) && ($freq<=0.45)) {
	  $bin=0.45;
	}
	if (($freq>0.45) && ($freq<=0.5)) {
	  $bin=0.5;
	}
	if (($freq>0.5) && ($freq<=0.55)) {
	  $bin=0.55;
	}
	if (($freq>0.55) && ($freq<=0.6)) {
	  $bin=0.6;
	}
	if (($freq>0.6) && ($freq<=0.65)) {
	  $bin=0.65;
	}
	if (($freq>0.65) && ($freq<=0.7)) {
	  $bin=0.7;
	}
	if (($freq>0.7) && ($freq<=0.75)) {
	  $bin=0.75;
	}
	if (($freq>0.75) && ($freq<=0.8)) {
	  $bin=0.8;
	}
	if (($freq>0.8) && ($freq<=0.85)) {
	  $bin=0.85;
	}
	if (($freq>0.85) && ($freq<=0.9)) {
	  $bin=0.9;
	}
	if (($freq>0.9) && ($freq<=0.95)) {
	  $bin=0.95;
	}
	if ($freq>0.95) {
	  $bin=1;
	}
	return $bin;
      }
 


      sub zbin {
	$freq=$_[0];
	if ($freq<=-2) {
	  $zbin=-2;
	}
	if (($freq>-2) && ($freq<=-1.5)) {
	  $zbin=-1.5;
	}
	if (($freq>-1.5) && ($freq<=-1)) {
	  $zbin=-1;
	}
	if (($freq>-1) && ($freq<=-0.5)) {
	  $zbin=-0.5;
	}
	if (($freq>-0.5) && ($freq<=0)) {
	  $zbin=-0;
	}
	if (($freq>0) && ($freq<=0.5)) {
	  $zbin=0.5;
	}
	if (($freq>0.5) && ($freq<=1)) {
	  $zbin=1;
	}
	if (($freq>1) && ($freq<=1.5)) {
	  $zbin=1.5;
	}
	if (($freq>1.5) && ($freq<=2)) {
	  $zbin=2;
	}
	if (($freq>2) && ($freq<=2.5)) {
	  $zbin=2.5;
	}
	if (($freq>2.5) && ($freq<=3)) {
	  $zbin=3;
	}
	if (($freq>3) && ($freq<=3.5)) {
	  $zbin=3.5;
	}
	if (($freq>3.5) && ($freq<=4)) {
	  $zbin=4;
	}
	if (($freq>4) && ($freq<=4.5)) {
	  $zbin=4.5;
	}
	if (($freq>4.5) && ($freq<=5)) {
	  $zbin=5;
	}
	if (($freq>5) && ($freq<=5.5)) {
	  $zbin=5.5;
	}
	if (($freq>5.5) && ($freq<=6)) {
	  $zbin=6;
	}
	if (($freq>6) && ($freq<=6.5)) {
	  $zbin=6.5;
	}
	if (($freq>6.5) && ($freq<=7)) {
	  $zbin=7;
	}
	if (($freq>7) && ($freq<=7.5)) {
	  $zbin=7.5;
	}
	if (($freq>7.5) && ($freq<=8)) {
	  $zbin=8;
	}

	  return $zbin;
      }
