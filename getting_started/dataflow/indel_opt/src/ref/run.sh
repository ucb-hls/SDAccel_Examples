./test 25 > out 
./test 125 >> out 
./test 525 >> out 
./test 625 >> out 
./test 725 >> out 
./test 825 >> out 
./test 925 >> out 
./test 1025 >> out 

vimdiff -c 'set diffopt+=iwhite' out ./ir_toy/GOLDEN.txt 
