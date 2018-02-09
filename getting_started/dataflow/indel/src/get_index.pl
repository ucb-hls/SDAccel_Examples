my @files = glob "./ch22-ir/*.SEQ.tbl";
print scalar @files;
print "\n";
for (0..$#files){
  $files[$_] =~ s/\.SEQ\.tbl//;  
  $files[$_] =~ s/\.\/ch22-ir\///;  
}

foreach $file (@files){
	printf("$file ");
}
