#! /usr/bin/perl -w



#$infile = "/proj/b2013119/private/Analyses/PHYLOP/46WAY/chrUn_gl000244.phyloP46way.placental.wigFix";
#$outfile = "/Users/erikaxelsson/Documents/Projects/Dog_resequencing/Post_2012_publication/Novartis/Analyses/ANNOTATIONS/phylop/chr2_b.phyloP100way.bed";

for (</proj/uppstore2017236/b2013119/private/Analyses/PHYLOP/100WAY/*.wigFix>){
  $infile=$_;
  print "$infile\n\n";
#  $count=0;
  open IN, $infile;
  while (<IN>) {
    $line=$_;
    if ($line=~/fixed/) {
      @line = split (/\s/,$line);
#      print "@line\n";
      if ($line[1]=~ /chrom\=(.+)/) {
	$chr=$1;
	$outfile = "/proj/uppstore2017236/b2013119/private/Analyses/PHYLOP/100WAY/".$chr."_hg38.phyloP100way.bed";
	$outfile=~s/\.wigFix/\.bed/;
	print "$outfile\n";
	close OUT;
	open OUT, ">>$outfile";
      }
      if ($line[2]=~ /start\=(.+)/) {
	$pos=$1-1; ### Wig fix format is 1-based, BED format is 0-based
	$end=$pos+1;
      }
#      print "$chr\t$pos\n";
    } else {
      print OUT "$chr\t$pos\t$end\t$line";
      $pos++;
      $end++;
    }
#    $count++; 
  }
  close IN;
  close OUT;
}

