#! /usr/bin/perl -w

@chromosomes=(4..9);
#push @chromosomes, 'X';
#push @chromosomes, 'M';
push @chromosomes, 'Un';

foreach $chr (sort @chromosomes){
    $CHROMOSOMES{$chr}=1;
}

for (</proj/uppstore2017236/b2013119/private/dmytro/VCF/chrFiltered/*.vcf>) {
  $sv_file = $_;
  print "$sv_file\n";
  if ($sv_file=~/proj\/uppstore2017236\/b2013119\/private\/dmytro\/VCF\/chrFiltered\/Arctic_vs_Agrarian_merged_markDupl_BQSR_chrAll_SNPs_VQSRfilterPASSED_chr(.+?)\.vcf/) {
    $chromosome=$1;
  }
  if (exists $CHROMOSOMES{$chromosome}){ 
      print "$chromosome\n";
      $out_file="/proj/uppstore2017236/b2013119_nobackup/Liftover/Agrarian_Arctic_SNPs_BED/canfam3_chr".$chromosome."_all_SNPs.bed";
      print "$out_file\n";
      
      ### load variants  ###
      open OUT, ">>$out_file";
      open IN, $sv_file;
      while (<IN>) {
	  if (/^chr/) {
	      $count++;
	      @sv = split, /\t/;
	      $chr=$sv[0];
	      $pos=$sv[1]-1; ### vcf format is 1-based, BED format is 0-based 
	      $end=$pos+1;
	      print OUT "$chr\t$pos\t$end\t$chr\t$pos\n";
	  }
      }
      close OUT; 
  }
}

