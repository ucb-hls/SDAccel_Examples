#!/usr/bin/perl
#use strict;
#use warnings;
use File::Path qw( make_path );
use File::Spec;
use File::Copy;

use Data::Dumper qw(Dumper);

my $dir = "./ch22_out" ;
if ( !-d $dir) {
    make_path $dir or die "Failed to create path: $dir";
}
local $/ = undef;

my $fn = "index.txt";
open FILE, "$fn" or die "Couldn't open file: $!";
binmode FILE;
my $string = <FILE>;
#print ($string);
my @indices = split / /, $string;
foreach my $index (@indices){
	$index =~ s/\n+//g;	
	#$index =~ tr/\r\n//d;
	$index =~ s/\s+//g;	
	print $index;
#	print $index;	
	system("./host $index > ./ch22_out/$index\.out ");
	#print("./host ".$index." > ./ch22_out/".$index."\.out ");
	copy("./sdaccel_profile_summary.csv","./ch22_out/$index\.sdaccel_profile_summary\.csv") or die "Copy failed: $!";
}
close FILE;
