
for TEST in `seq 15 19`;do
    #TEST=19
    SORTED=
    #SORTED=sorted_
    SCHE=/home/centos/target-dir/${TEST}valid.txt 
    #SCHE=/home/centos/aws-fpga/SDAccel/examples/xilinx/getting_started/dataflow/indel_simple_multi_sorted/schedule_${SORTED}${TEST}.txt 
    TARGET=/home/centos/target-dir/$TEST/
    ERR=out_${SORTED}$TEST
    OUT=time_${SORTED}$TEST
    ./host $SCHE $TARGET >$OUT 2>$ERR 
    cp sdaccel_profile_summary.csv sdaccel_profile_summary.csv.${SORTED}${TEST}.full.16
done