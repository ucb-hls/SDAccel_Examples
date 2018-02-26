COMMON_REPO := ../../../

include $(COMMON_REPO)/utility/boards.mk
include $(COMMON_REPO)/libs/xcl2/xcl2.mk
include $(COMMON_REPO)/libs/opencl/opencl.mk

# Kernel profile 
PROFILE:=yes

# Host Application
host_SRCS=./src/host.cpp $(xcl2_SRCS) 
host_HDRS=$(xcl2_HDRS) ./src/Indel_Accel.h
host_CXXFLAGS=-I./src/ $(xcl2_CXXFLAGS) $(opencl_CXXFLAGS) -I/opt/Xilinx/SDx/2017.1.op/Vivado_HLS/include/ -g --debug
host_LDFLAGS=$(opencl_LDFLAGS)

EXES=host

# Kernel
# OpenCL file
Indel_Accel_SRCS=./src/Indel_Accel.cl

#Indel_Accel_CLFLAGS= -I./src/   --report estimate -g --debug --max_memory_ports Indel_Accel
Indel_Accel_CLFLAGS= -I./src/  --debug --max_memory_ports Indel_Accel
#Indel_Accel_CLFLAGS= -I/usr/include/linux/ -I/opt/Xilinx/SDx/2017.1.op/ext/glibc-2.17/glibc-2.17/ -I/opt/Xilinx/SDx/2017.1.op/ext/glibc-2.17/glibc-2.17/include/ -I/opt/Xilinx/SDx/2017.1.op/Vivado_HLS/include/ --report estimate -g --debug --max_memory_ports Indel_Accel 
Indel_Accel_LDCLFLAGS= --xp misc:map_connect=add.kernel.Indel_Accel_1.M_AXI_GMEM0.core.OCL_REGION_0.M00_AXI --xp misc:map_connect=add.kernel.Indel_Accel_1.M_AXI_GMEM1.core.OCL_REGION_0.M01_AXI --xp misc:map_connect=add.kernel.Indel_Accel_1.M_AXI_GMEM2.core.OCL_REGION_0.M02_AXI --xp vivado_param:project.runs.noReportGeneration=0
#Indel_Accel_NDEVICES=xilinx:adm-pcie-7v3:1ddr



XOS=Indel_Accel

# xclbin
Indel_Accel_XOS=Indel_Accel

XCLBINS=Indel_Accel

# check
check_EXE=host
#check_ARGS=125
check_XCLBINS=Indel_Accel

CHECKS=check

include $(COMMON_REPO)/utility/rules.mk

