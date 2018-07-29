#!/usr/bin/perl
# input  .tsv  tar_dir_path
my $infile = $ARGV[0];
open(INPUT, "<$infile");
my $outdir = $ARGV[1];

while(<INPUT>)
{
  my $line = $_;
  chomp($line);
  my @columns = split /\t/, $line;
  if ($line =~ /contig/) {
  } else {
      my $i = $columns[0];
      my $outfile = "$i.TARGET.tbl";
      open(OUTPUT, ">$outdir/$outfile");
      print OUTPUT "$columns[0]|$columns[1]|$columns[2]|$columns[3]|";
      close(OUTPUT);
  }
}
close(INPUT);
