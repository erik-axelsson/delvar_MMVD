#! /usr/bin/perl -w

#$infile = "/Users/erikaxelsson/Documents/Projects/Dog_resequencing/Post_2012_publication/Novartis/Analyses/ANNOTATIONS/phylop/chrM.phyloP100way.wigFix";
$outfile = "/proj/uppstore2017236/b2013119_nobackup/Liftover/PHYLOP/agrarian_arctic_canfam3_SNPs_with_100WAY_phylop.txt";
$phylop_stats = "/proj/uppstore2017236/b2013119_nobackup/Liftover/PHYLOP/agrarian_arctic_liftover_stats.txt";
open STATS, ">>$phylop_stats";

### Get phylop values for human hg38 positions that are segregating in dog. ### 
$phylo_count=0;
for (</proj/uppstore2017236/b2013119_nobackup/Liftover/Agrarian_Arctic_SNPs/phylop_100way_for_hg38*>){
  $infile=$_;
#  $outfile = $infile;
  print "$infile\n";
  $count=0;
  open IN, $infile;
  while (<IN>) {
    @line = split (/\s/);
#    print "@line\n";
    $chr=$line[0];
    $chr=~s/chr//;
    $PHYLOP_HUM_POS{$chr}{$line[1]}=$line[3]; ###MANISKO POS!!!
    $phylo_count++;
  }
  close IN;
}
print STATS "$phylo_count dog SNPs with phylop values\n";


### Get Canfam3 and hg38 position of all dog SNPs ###
$lifover=0;
for (</proj/uppstore2017236/b2013119_nobackup/Liftover/Agrarian_Arctic_SNPs/*canfam3_coordinates_for_all_SNPs.bed>) {
  $infile=$_;
  print "$infile\n";
  open IN, $infile;
  while (<IN>) {
    @line = split (/\s/);
    #    print "@line\n";
    $chr=$line[0]; #human chromosome
    $chr=~s/chr//;
    $pos=$line[1]; #human pos
    $dog_chr=$line[3]; #dog chromosome
    $dog_chr=~s/chr//;
    $dog_pos=$line[4]; #human pos
    $CONVERT_CHR{$chr}{$pos} = $dog_chr; #human chr, human pos -> dog chr. ##MANNISKO POS!!!!
    $CONVERT_POS{$chr}{$pos} = $dog_pos; #human chr, human pos -> dog pos.
    $liftover++;
  }
  close IN;
}
print STATS "dog SNPs with human coordinates $liftover\n";


### Convert phylop values from human to dog coordinates ###
foreach $chr (sort keys %PHYLOP_HUM_POS) {
  foreach $pos ( sort keys %{$PHYLOP_HUM_POS{$chr}}) {
    if (exists $CONVERT_CHR{$chr}{$pos}) {
      $dog_chr=$CONVERT_CHR{$chr}{$pos};
      $dog_pos=$CONVERT_POS{$chr}{$pos};
      $PHYLOP_DOG_POS{$dog_chr}{$dog_pos}=$PHYLOP_HUM_POS{$chr}{$pos};
      $DOG_TO_HUMAN_CHR{$dog_chr}{$dog_pos}=$chr;
      $DOG_TO_HUMAN_POS{$dog_chr}{$dog_pos}=$pos;
    }
  }
}

### Get Canfam3 SNPs that map back to original position after liftover to human and back to dog again ###
$reciprocal=0;
for (</proj/uppstore2017236/b2013119_nobackup/Liftover/Agrarian_Arctic_SNPs/and_back_to_human/*all_SNPs_canfam3_coordinates.bed>) {
  $infile=$_;
  print "$infile\n";
  open IN, $infile;
  while (<IN>) {
    chomp;
    @line = split (/\s/);
    $chr=$line[0]; #dog chromosome
    $chr=~s/chr//;
    $dog_chr=$line[3]; #dog chromosome
    $dog_chr=~s/chr//;
    $pos=$line[1]; #dog chromosome
    $pos=~s/\s//;
    $dog_pos=$line[4]; #dog chromosome
    $dog_pos=~s/\s//;
    if ($chr eq $dog_chr) {
      if ($pos==$dog_pos) {
	$TRUE{$chr}{$pos}=1;  ###HUND pos!!!
	$reciprocal++;
      }
    }
  }
}
print STATS "SNPs that map reciprocally in dog: $reciprocal\n";

### Filter and print phylop values for those with canfam3 coordinate that map back to  original dog position after liftover back and fort to human genome ###
open OUT, ">>$outfile";
foreach $chr (sort numerically keys %PHYLOP_DOG_POS) {
  foreach $pos (sort numerically keys %{$PHYLOP_DOG_POS{$chr}}) {
    #   print "$chr\t$pos\n";
    #print STDERR "    if (exists $CONVERT{$chr}{$pos}) {\n";
    #print STDERR "      if (exists $TRUE{$chr}{$pos}) \n";
    if (exists $TRUE{$chr}{$pos}) {
      $dog_pos=$pos+1;  ### BED files have 0-based coordinates, that are turned back to 1-based coordinates here.
      $hum_pos=$DOG_TO_HUMAN_POS{$chr}{$pos}+1;
      print OUT "$chr\t$dog_pos\t$DOG_TO_HUMAN_CHR{$chr}{$pos}\t$hum_pos\t$PHYLOP_DOG_POS{$chr}{$pos}\n"; ## dog chr, dog pos, human chr, human pos, phylop value
    }
  }
}


=pod
#open OUT, ">>$outfile";
foreach $chr (sort keys %CONVERT) {
  foreach $pos (sort keys %{$CONVERT{$chr}}) {
    print "$chr\t$pos\t@{$CONVERT{$chr}{$pos}}\n";
  }
}
=cut

sub numerically {
    $a<=>$b;
}
