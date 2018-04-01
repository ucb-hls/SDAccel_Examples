COMMON_REPO := ../../../

include $(COMMON_REPO)/utility/boards.mk
include $(COMMON_REPO)/libs/xcl2/xcl2.mk
include $(COMMON_REPO)/libs/opencl/opencl.mk
include $(COMMON_REPO)/libs/oclHelper/oclHelper.mk

# Kernel profile 
PROFILE:=yes

# Host Application
host_SRCS=./src/host.cpp $(xcl2_SRCS) $(oclHelper_SRCS) 
host_HDRS=$(xcl2_HDRS) ./src/Indel_Accel.h
host_CXXFLAGS=-I./src/ $(xcl2_CXXFLAGS) $(opencl_CXXFLAGS)  $(oclHelper_CXXFLAGS)  -I/opt/Xilinx/SDx/2017.1.op/Vivado_HLS/include/ -g --debug
host_LDFLAGS=$(opencl_LDFLAGS)

EXES=host

# Kernel
Indel_Accel_SRCS=./src/Indel_Accel.cpp
Indel_Accel_CLFLAGS= --optimize 3 --kernel_frequency 150  -k Indel_Accel -I./src/ --report estimate -g --debug  #--max_memory_ports Indel_Accel

#Indel_Accel_CLFLAGS= -k Indel_Accel -I./src/ --report estimate -g --debug
#Indel_Accel_LDCLFLAGS= --nk Indel_Accel:10 --xp vivado_param:project.runs.noReportGeneration=0 --xp misc:map_connect=add.kernel.Indel_Accel_1.M_AXI_GMEM0.core.OCL_REGION_0.M00_AXI --xp misc:map_connect=add.kernel.Indel_Accel_1.M_AXI_GMEM1.core.OCL_REGION_0.M01_AXI --xp misc:map_connect=add.kernel.Indel_Accel_1.M_AXI_GMEM2.core.OCL_REGION_0.M02_AXI --xp misc:map_connect=add.kernel.Indel_Accel_1.M_AXI_GMEM3.core.OCL_REGION_0.M03_AXI  --xp misc:map_connect=add.kernel.Indel_Accel_2.M_AXI_GMEM0.core.OCL_REGION_0.M00_AXI --xp misc:map_connect=add.kernel.Indel_Accel_2.M_AXI_GMEM1.core.OCL_REGION_0.M01_AXI --xp misc:map_connect=add.kernel.Indel_Accel_2.M_AXI_GMEM2.core.OCL_REGION_0.M02_AXI --xp misc:map_connect=add.kernel.Indel_Accel_2.M_AXI_GMEM3.core.OCL_REGION_0.M03_AXI --xp misc:map_connect=add.kernel.Indel_Accel_3.M_AXI_GMEM0.core.OCL_REGION_0.M00_AXI --xp misc:map_connect=add.kernel.Indel_Accel_3.M_AXI_GMEM1.core.OCL_REGION_0.M01_AXI --xp misc:map_connect=add.kernel.Indel_Accel_3.M_AXI_GMEM2.core.OCL_REGION_0.M02_AXI --xp misc:map_connect=add.kernel.Indel_Accel_3.M_AXI_GMEM3.core.OCL_REGION_0.M03_AXI --xp misc:map_connect=add.kernel.Indel_Accel_4.M_AXI_GMEM0.core.OCL_REGION_0.M00_AXI --xp misc:map_connect=add.kernel.Indel_Accel_4.M_AXI_GMEM1.core.OCL_REGION_0.M01_AXI --xp misc:map_connect=add.kernel.Indel_Accel_4.M_AXI_GMEM2.core.OCL_REGION_0.M02_AXI --xp misc:map_connect=add.kernel.Indel_Accel_4.M_AXI_GMEM3.core.OCL_REGION_0.M03_AXI 
#Indel_Accel_LDCLFLAGS= --nk Indel_Accel:1 --xp vivado_param:project.runs.noReportGeneration=0  --xp misc:map_connect=add.kernel.Indel_Accel_1.M_AXI_GMEM0.core.OCL_REGION_0.M00_AXI #--xp misc:map_connect=add.kernel.Indel_Accel_2.M_AXI_GMEM0.core.OCL_REGION_0.M01_AXI --xp misc:map_connect=add.kernel.Indel_Accel_3.M_AXI_GMEM0.core.OCL_REGION_0.M02_AXI --xp misc:map_connect=add.kernel.Indel_Accel_4.M_AXI_GMEM0.core.OCL_REGION_0.M03_AXI  --xp misc:map_connect=add.kernel.Indel_Accel_5.M_AXI_GMEM0.core.OCL_REGION_0.M00_AXI --xp misc:map_connect=add.kernel.Indel_Accel_6.M_AXI_GMEM0.core.OCL_REGION_0.M01_AXI --xp misc:map_connect=add.kernel.Indel_Accel_7.M_AXI_GMEM0.core.OCL_REGION_0.M02_AXI --xp misc:map_connect=add.kernel.Indel_Accel_8.M_AXI_GMEM0.core.OCL_REGION_0.M03_AXI  
Indel_Accel_LDCLFLAGS= --optimize 3 --kernel_frequency  150 --nk Indel_Accel:16 --xp vivado_param:project.runs.noReportGeneration=0   #--xp misc:map_connect=add.kernel.Indel_Accel_1.M_AXI_GMEM0.core.OCL_REGION_0.M00_AXI --xp misc:map_connect=add.kernel.Indel_Accel_2.M_AXI_GMEM0.core.OCL_REGION_0.M01_AXI --xp misc:map_connect=add.kernel.Indel_Accel_3.M_AXI_GMEM0.core.OCL_REGION_0.M02_AXI --xp misc:map_connect=add.kernel.Indel_Accel_4.M_AXI_GMEM0.core.OCL_REGION_0.M03_AXI  --xp misc:map_connect=add.kernel.Indel_Accel_5.M_AXI_GMEM0.core.OCL_REGION_0.M00_AXI --xp misc:map_connect=add.kernel.Indel_Accel_6.M_AXI_GMEM0.core.OCL_REGION_0.M01_AXI --xp misc:map_connect=add.kernel.Indel_Accel_7.M_AXI_GMEM0.core.OCL_REGION_0.M02_AXI --xp misc:map_connect=add.kernel.Indel_Accel_8.M_AXI_GMEM0.core.OCL_REGION_0.M03_AXI  



Indel_Accel_NDEVICES=xilinx:adm-pcie-7v3:4ddr

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

