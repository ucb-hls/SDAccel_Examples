COMMON_REPO := ../../../

include $(COMMON_REPO)/utility/boards.mk
include $(COMMON_REPO)/libs/xcl2/xcl2.mk
include $(COMMON_REPO)/libs/opencl/opencl.mk

# Host Application
host_SRCS=./src/host.cpp $(xcl2_SRCS) 
host_HDRS=$(xcl2_HDRS) ./src/Indel_Accel.h
host_CXXFLAGS=-I./src/ $(xcl2_CXXFLAGS) $(opencl_CXXFLAGS)
host_LDFLAGS=$(opencl_LDFLAGS)

EXES=host

# Kernel
Indel_Accel_SRCS=./src/Indel_Accel.cpp
Indel_Accel_CLFLAGS= -k Indel_Accel -I./src/

XOS=Indel_Accel

# xclbin
Indel_Accel_XOS=Indel_Accel

XCLBINS=Indel_Accel

# check
check_EXE=host
check_XCLBINS=Indel_Accel

CHECKS=check

include $(COMMON_REPO)/utility/rules.mk

