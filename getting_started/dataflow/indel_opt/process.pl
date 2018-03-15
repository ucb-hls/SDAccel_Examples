#!/usr/bin/perl
#use strict;
#use warnings;
use File::Path qw( make_path );

  my $path = "tmp.txt";
	my $success = 0;
  my @newlines;
  for (my $i=0; $i < 4; $i++){
	if(!open FILE,"<", $path){
		print $!;
	} else {
    push(@newlines, "//INST_$i============================================//\n");
		while(<FILE>){
			$line = $_;
      #chomp ($line);  
      $line =~ s/_X/_$i/g;
      $line =~ s/_\|X/_X/g;
      $line =~ s/_\~X/$i/g;
      push(@newlines, $line);
		}
		close FILE;
	}

  }
  open(FILE, ">tmp_out.txt") || die "Error!";
  print FILE @newlines;
	close(FILE)


