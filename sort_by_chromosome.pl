#! /usr/bin/perl -w


#for (</Users/erikaxelsson/Documents/Projects/Dog_resequencing/Post_2012_publication/Novartis/Analyses/SNP_INDEL/canfam3_chr38_all_SNPs.bed>) {
for (</proj/uppstore2017236/b2013119_nobackup/Liftover/Agrarian_Arctic_SNPs/*coordinates.bed>) {
  $dog_bed_file = $_;
  print "$dog_bed_file\n";
  if ($dog_bed_file=~/chr(.+?)\_/) {
    $chromosome=$1;
  }
#  print "$chromosome\n";
  open IN, $dog_bed_file;
  while (<IN>) {
    chomp;
    @line=split("\t", $_);
#    print "$line[0]\n";
    $chr=$line[0];
    $pos=$line[1];
#    print "@line\n";
    push (@{$COORDINATES{$chr}{$pos}}, @line);
  }
}
close IN;
#foreach $chr (sort keys %COORDINATES) {
#  foreach $pos (sort numerically keys %{$COORDINATES{$chr}}) {
#    print "$chr\t$pos\t@{$COORDINATES{$chr}{$pos}}\n";
#  }
#}
#die;

foreach $chr (sort keys %COORDINATES) {
  $outfile = "/proj/uppstore2017236/b2013119_nobackup/Liftover/Agrarian_Arctic_SNPs/hg38_".$chr."_canfam3_coordinates_for_all_SNPs.bed";
  $grep_pattern_file = "/proj/uppstore2017236/b2013119_nobackup/Liftover/Agrarian_Arctic_SNPs/hg38_".$chr."_canfam3_coordinates_for_all_SNPs.grep.pattern.txt";
  open OUT, ">>$outfile";
  open UT, ">>$grep_pattern_file";
  foreach $pos (sort numerically keys %{$COORDINATES{$chr}}) {
    print OUT "@{$COORDINATES{$chr}{$pos}}\n";
    print UT "$chr\t$pos\t\n";
  }
}

sub numerically {
  $a<=>$b;
}

