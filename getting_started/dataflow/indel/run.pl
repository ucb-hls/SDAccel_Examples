#!/usr/bin/perl
#use strict;
#use warnings;
use Data::Dumper qw(Dumper);
local $/ = undef;
open FILE, "index.txt" or die "Couldn't open file: $!";
binmode FILE;
my $string = <FILE>;
#print ($string);
my @indices = split / /, $string;
foreach my $index (@indices){
#	print $index;	
	system("./host $index >> ch22_out")
}
close FILE;
