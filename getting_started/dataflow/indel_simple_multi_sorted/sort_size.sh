#source $XILINX_SDX/settings64.sh

#TEST=20
#for TEST in `seq 15 21`;do
for TEST in `seq 1 7`;do

SCHE=/home/centos/target-dir/${TEST}valid.txt
TARGET=/home/centos/target-dir/$TEST/
ERR=count_${TEST}.txt
OUT=schedule_sorted_${TEST}.txt

#XCL_EMULATION_MODE=sw_emu ./host /home/centos/target-dir/21valid.txt /home/centos/target-dir/21/ 2>21schedule.txt
XCL_EMULATION_MODE=sw_emu ./host $SCHE $TARGET 2>$ERR
python sort.py $ERR > $OUT
done
