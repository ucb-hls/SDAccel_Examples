#!/usr/bin/perl
#use strict;
#use warnings;
use File::Path qw( make_path );
use File::Spec;
use File::Copy;

use Text::CSV;

#use Data::Dumper qw(Dumper);
#local $/ = undef;

scalar (@ARGV) == 2 or die "RUN perl parse.pl  ./ indel_simple_multi";

my $dir = $ARGV[0];
my $report = $ARGV[1];
#my $dir = "./ch22_out/" ;
#my $fn = "../indel_tests/ch22-schedule.txt";
#my $fn = "../indel_tests/short.txt";
#my $report = "indel_opt_1ddr";

#open FILE, "$fn" or die "Couldn't open file: $!";
#binmode FILE;
#my $string = <FILE>;
#print ($string);
#my @indices = split / /, $string;

my @indices = (1..22);
open RSLT, ">$report.csv";

my @failed = ();
my $print_title = 1;
foreach my $index (@indices){
	$index =~ s/\n+//g;	
	$index =~ s/\s+//g;	
	if (not $index =~ m/(\d+)/g) {next; }
        # time_10
	my $out_path = $dir."time_$index";
	#my ($success, $parse_time, $prep_time, $exe_time) = parseOutput($out_path);
	my ($success, $parse_time, $prep_time, $exe_time) = parseOutput($out_path);

	my $prof_path = $dir."sdaccel_profile_summary.csv.$index.full.16";
	(undef, undef, undef, my $exe_time) = parseOutput($prof_path);
	
	#print "$success, $parse_time, $prep_time, $exe_time\n";
	#my $total_time = $parse_time + $prep_time + $exe_time;
	print "$parse_time\n";
	print "$exe_time\n";
	
	my $total_time = $exe_time + $parse_time; 
	#my $prof_path = $dir.$index.".sdaccel_profile_summary.csv";
	#sdaccel_profile_summary.csv.21.full.16

	my ($cl_result_ptr , $ot_result_ptr) = parseReport($prof_path); 
	%cl_result = %$cl_result_ptr;
	%ot_result = %$ot_result_ptr;
	#print "Key: $_ and Value: $cl_result{$_}\n" foreach (keys%cl_result);
	#print "\n";
	#print "Key: $_ and Value: $ot_result{$_}\n" foreach (keys%ot_result);
	#print "\n";
	#print %cl_result;
	
	if ($print_title eq 1) {
		#print RSLT "Target #, Success, Total time, Parse Time, Preprocessing Time, Execution Time, "; 	
		print RSLT "Target #, Success, Total time, Parse Time, Execution Time, "; 	
		print RSLT "cl:clEnqueueMigrateMemObjects, Host-Global Read, Host-Global Write, "; 
		print RSLT "Kernel-Global Read, Kernel-Global Write";
		print RSLT "\n";
		$print_title = 0;
	}
	if ($success ne 1){ push(@failed, $index);}
	#else { 	
		#print RSLT "$index, $success, $total_time, $parse_time, $prep_time, $exe_time, ";
		print RSLT "$index, $success, $total_time, $parse_time, $exe_time, ";
		print RSLT $cl_result{"clEnqueueMigrateMemObjects"}->[1].", "; 

		print RSLT $ot_result{"4_2"}->[6].", ";
		print RSLT $ot_result{"4_3"}->[6].", ";

		print RSLT ($ot_result{"5_2"}->[6]/ 1000) .", ";
		print RSLT ($ot_result{"5_3"}->[6]/ 1000) .", ";
		print RSLT "\n";	
	#}
}

	print "Failed Tests: ";
	print "$_, " foreach (@failed);
	print "\n";
close RSLT;

#################################Parse Output###################################
sub parseOutput {

	my $path = $_[0];
	my $success = 0;
	my $parse_time = undef;
	my $prep_time = undef;
	my $exe_time = undef;
	if(!open FILE,"<:encoding(utf8)", $path){
		print $!;
	} else {
		while(<FILE>){
			$line = $_;
			chomp ($line);  
			#print $line;
			#print("==================\n");
			#if($line =~m /Parsing time is :(\d+) ms/){
			if($line =~m /Parse time is: (\d+) ms/){
				$parse_time = $1;
			}
			if($line =~m /Preprocess time is : (\d+) ms/){
				$prep_time = $1; 
			}
			if($line =~m /OpenCl Execution time is: (\d+\.?\d*) ms/){
				$exe_time = $1;
			}
			if($line =~m /TEST PASSED/){
				$success = 1;	
			}
			#DEVICE_EXEC_TIME,xilinx:aws-vu9p-f1:4ddr-xpr-2pr:4.0-0,6092875.350452
			if($line =~ m /DEVICE_EXEC_TIME,xilinx:aws-vu9p-f1:4ddr-xpr-2pr:4.0-0,(\S+),/){
				$exe_time = $1
			}
		}
		close FILE;
	}
	

	return ($success, $parse_time, $prep_time, $exe_time); 
}


#################################Parse Report###################################
sub parseReport {

	my $csv = Text::CSV->new({ sep_char => ',' });

	my $path = $_[0];
	my %cl_result; 
	my %ot_result;
	my $section = 0;
	my $local = -1; 
	my $space = 0;

	open(my $data, '<', $path) or die "Could not open '$path' $!\n";

		while (my $line = <$data>) {
			#$line = $_;
			chomp ($line);  
			#else {
			#	warn "Line could not be parsed: $line\n";
			#}


			if ($line =~ m/^\s*$/ and ($section ne 0)) {
				$local = 0;
				#print "set to local\n";
			}
			#Kernel,Number Of Enqueues,Total Time (ms),Minimum Time (ms),Average Time (ms),Maximum Time (ms),
			#if ($csv->parse($line) and ($section ne 0) and ($section ne 1)) {
			elsif (($section ne 0) and ($local >= 0) and $csv->parse($line)) {
				if ($local eq 0) { $section = $section + 1;}
				#print"ot $section _ $local $line\n";
				my @fields = $csv->fields();
				my $size = @fields;
				#print("size $size\n");
				if ($size < 1){ next; }
				$ot_result{$section."_".$local} = [@fields];
				$local = $local + 1;
			}
			#API Name,Number Of Calls,Total Time (ms),Minimum Time (ms),Average Time (ms),Maximum Time (ms),
			elsif ($csv->parse($line) and ($section eq 1)) {
				#print"cl $line\n";
			#if ($section eq 1) {
				my @fields = $csv->fields();
				my $size = @fields;
				if ($size < 1){ next; }
				my $cl_entry = shift @fields;
				$cl_result{$cl_entry} = [@fields];
			#my @value_array = @{$new_hash{$_}};
			}
			#
			elsif($line =~ m/OpenCL API Calls/){
				$section = 1; 	
			}
		}



		#close FILE;
	return (\%cl_result, \%ot_result);

}
