BM=Indel_Accel.hw.xilinx_aws-vu9p-f1_4ddr-xpr-2pr_4_0
XCLBIN=/home/centos/aws-fpga/SDAccel/examples/xilinx/getting_started/dataflow/indel/xclbin/

$SDACCEL_DIR/tools/create_sdaccel_afi.sh -xclbin=$XCLBIN/$BM.xclbin -o=$BM \
        -s3_bucket="hqjenny-100" -s3_dcp_key=dcp/$BM -s3_logs_key=logs/
