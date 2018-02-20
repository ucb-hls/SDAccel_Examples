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
Indel_Accel_SRCS=./src/Indel_Accel.cpp
Indel_Accel_CLFLAGS= -k Indel_Accel -I./src/ --report estimate -g --debug --sp  Indel_Accel.m_axi_gmem0:bank0 --sp  Indel_Accel.m_axi_gmem1:bank1 --sp  Indel_Accel.m_axi_gmem2:bank2


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

