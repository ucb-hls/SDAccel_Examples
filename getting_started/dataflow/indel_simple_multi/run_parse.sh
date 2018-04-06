
#TEST=22
for TEST in `seq 1 4`;do
    #SORTED=
    #SORTED=sorted_
    SCHE=/home/centos/target-dir/${TEST}valid.txt 
    #SCHE=../indel_tests/ch22-schedule.txt 
    #SCHE=/home/centos/aws-fpga/SDAccel/examples/xilinx/getting_started/dataflow/indel_simple_multi_sorted/schedule_${SORTED}${TEST}.txt 
    TARGET=/home/centos/target-dir/$TEST/
    #TARGET=../indel_tests/ch22-ir/
    ERR=out_${SORTED}$TEST
    OUT=time_${SORTED}$TEST
    ./host $SCHE $TARGET >$OUT 2>$ERR 
    #cp sdaccel_profile_summary.csv sdaccel_profile_summary.csv.${SORTED}${TEST}.full.16
done
