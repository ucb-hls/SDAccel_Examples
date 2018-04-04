#!/bin/bash
echo "" > sw_out

#array=(1 25 125 525 625 725 825 925 1025)
array=(25)
for value in "${array[@]}"
do 
    echo $value
    make check TARGETS=sw_emu DEVICES=$AWS_PLATFORM all  >> sw_out 
    #make check TARGETS=sw_emu DEVICES=$AWS_PLATFORM check_ARGS=$value  all  >> out  
    #make TARGETS=hw DEVICES=$AWS_PLATFORM all > out 
done
echo finish
