/**********
  Copyright (c) 2017, Xilinx, Inc.
  All rights reserved.

  Redistribution and use in source and binary forms, with or without modification,
  are permitted provided that the following conditions are met:

  1. Redistributions of source code must retain the above copyright notice,
  this list of conditions and the following disclaimer.

  2. Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

  3. Neither the name of the copyright holder nor the names of its contributors
  may be used to endorse or promote products derived from this software
  without specific prior written permission.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
  THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
  IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
  OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
  EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 **********/
//OpenCL utility layer include
#include "xcl2.hpp"
#include <vector>
#include <assert.h>
#include <map>

#include <pthread.h>
#include <chrono>
#include "ap_int.h"

#include "Indel_Accel.h"
#include "Indel_Accel_SW.cpp"

#define WORK_GROUP 4 
#define WORK_ITEM_PER_GROUP 1
#define PARALLEL_UNITS 4
#define NUM_KERNELS 8
// JENNY TODO
// Put input and output onto different memory banks 
// https://github.com/Xilinx/SDAccel_Examples/blob/master/getting_started/kernel_to_gmem/
//
//typedef std::map<char, ap_uint<8>> BasePairMap;
typedef std::map<char, char> BasePairMap;
BasePairMap m;

#include "Indel_Accel_Host.cpp"

void Indel_Rank (const int consensus_size, const int reads_size, int*  min_whd, int* __restrict new_ref_idx) {
    int new_ref[READS_SIZE];

    int min_score = 0x7fffffff;
    int min_idx = consensus_size + 1;
    int i, j;
    score: for (i = 1; i < consensus_size; i++) {
        int score = 0;
        for (j = 0; j < reads_size; j++) {
            printf("[%d,%d]-%d\t", i, j, min_whd[(i * reads_size + j) << 1]);
            int tmp = min_whd[(i * reads_size + j) << 1] - min_whd[j << 1];
            score += (tmp > 0) ? tmp: -tmp;
        }
        min_idx = (score < min_score) ? i : min_idx;
        //scores[i] = score;
    }
    //printf( "min_idx: %d\n", min_idx);
    //assert(min_idx < consensus_size);
    
    rank: for (j = 0; j < reads_size; j++) {
            new_ref[j] = min_whd[(min_idx * reads_size + j) << 1];
            new_ref_idx[j] = min_whd[((min_idx * reads_size +j) << 1) + 1];
            }

    print: for (j = 0; j < reads_size; j++) {
     //printf("Read %2d whd %2d index %2d\n", j, new_ref[j], new_ref_idx[j]);
        printf("Kernel: Read %2d whd %4d  index %2d\n", j, new_ref[j], new_ref_idx[j]);
    }
}


int main(int argc, char** argv)
{

    printf("Parse scheduele file\n");
    int num_tests = 0;
    int* test_indices = parse_schedule("../indel_tests/ir_toy-schedule.txt", &num_tests);

    //if (argc < 2){
    //    return 1;
    //}
    
    m['A'] = 0;
    m['T'] = 1;
    m['C'] = 2;
    m['G'] = 3;
    m['U'] = 4;
 
    //unsigned banks[4] = {XCL_MEM_DDR_BANK0, XCL_MEM_DDR_BANK1, XCL_MEM_DDR_BANK2, XCL_MEM_DDR_BANK3};
    //unsigned ddr_bank = XCL_MEM_DDR_BANK0 | XCL_MEM_DDR_BANK1 | XCL_MEM_DDR_BANK2 | XCL_MEM_DDR_BANK3;
    unsigned ddr_bank = XCL_MEM_DDR_BANK0;

    //OPENCL HOST CODE AREA START
    std::vector<cl::Device> devices = xcl::get_xil_devices();
    cl::Device device = devices[0];

    cl::Context context(device);
    //cl::CommandQueue q(context, device, CL_QUEUE_PROFILING_ENABLE | CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE);
    std::vector<cl::CommandQueue> qs( NUM_KERNELS, cl::CommandQueue(context, device, CL_QUEUE_PROFILING_ENABLE));
    //cl::CommandQueue q(context, device, CL_QUEUE_PROFILING_ENABLE);
    
    std::string device_name = device.getInfo<CL_DEVICE_NAME>(); 

    //Create Program and Kernel
    std::string binaryFile = xcl::find_binary_file(device_name, "Indel_Accel");
    cl::Program::Binaries bins = xcl::import_binary_file(binaryFile);
    devices.resize(1);
    cl::Program program(context, devices, bins);
    //cl::Kernel krnl_indel(program,"Indel_Accel");

    std::vector<cl::Kernel> krnl_indels;
    krnl_indels.push_back(cl::Kernel(program,"Indel_Accel"));
    krnl_indels.push_back(cl::Kernel(program,"Indel_Accel"));
    krnl_indels.push_back(cl::Kernel(program,"Indel_Accel"));
    krnl_indels.push_back(cl::Kernel(program,"Indel_Accel"));
    std::chrono::high_resolution_clock::time_point start, finish;
    std::chrono::milliseconds duration;

    // Used for printing results
    int* reads_size_arr = (int*) malloc(num_tests * sizeof(int));
    int ** new_ref_idx_ref_arr = (int**) malloc(num_tests * sizeof(int*));
    std::vector<std::vector<int,aligned_allocator<int>> *> new_ref_idx_arr(num_tests);

  for (int test_idx = 0; test_idx< num_tests; test_idx+= PARALLEL_UNITS) {
    int kernel_idx = ((test_idx / 4)) % NUM_KERNELS; 

    std::chrono::milliseconds parse_time[PARALLEL_UNITS];
    //int kernel_idx = test_idx % 2;
    //cl::Kernel  krnl_indel = krnl_indels[0];
    cl::Kernel &  krnl_indel = krnl_indels[kernel_idx];
    cl::CommandQueue & q = qs[kernel_idx];
    //char* test_num = arg;
    std::vector<cl::Memory> inBufVec, outBufVec;

    int * new_ref_idx_ref_0;
    int * new_ref_idx_ref_1;
    int * new_ref_idx_ref_2;
    int * new_ref_idx_ref_3;
    //int * new_ref_idx_0;
    //int * new_ref_idx_1;    
    //int * new_ref_idx_2;
    //int * new_ref_idx_3;
    int narg=0;
    int reads_size_0, reads_size_1, reads_size_2, reads_size_3;

    const char* file_prefix = "../indel_tests/ir_toy/";
    char test_num[5];
    char con[256] = "";
    char reads[256]="";

    //cl::CommandQueue q_X (context, device, CL_QUEUE_PROFILING_ENABLE);
    std::chrono::high_resolution_clock::time_point start, finish;
    std::chrono::milliseconds duration;
    char* con_arr, *reads_arr, *weights_arr; 
    int * con_len, con_size, *reads_len; //, reads_size; 
    int con_total_len, reads_total_len, weights_total_len; 
    int * min_whd, * min_whd_idx, * new_ref;  //* new_ref_idx_ref;
 
    //Run_Unit(0, test_idx, context, q, krnl_indel, inBufVec, outBufVec, new_ref_idx_ref_0, new_ref_idx_0, parse_time, narg, num_tests, test_indices, reads_size_0);
    //Run_Unit(1, test_idx, context, q, krnl_indel, inBufVec, outBufVec, new_ref_idx_ref_1, new_ref_idx_1, parse_time, narg, num_tests, test_indices, reads_size_1);
    //Run_Unit(2, test_idx, context, q, krnl_indel, inBufVec, outBufVec, new_ref_idx_ref_2, new_ref_idx_2, parse_time, narg, num_tests, test_indices, reads_size_2);
    //Run_Unit(3, test_idx, context, q, krnl_indel, inBufVec, outBufVec, new_ref_idx_ref_3, new_ref_idx_3, parse_time, narg, num_tests, test_indices, reads_size_3);


//INST_0============================================//
    int test_idx_0= test_idx + 0;

    printf("test_idx_0: %d\n", test_idx_0);
    //int reads_size_0;
    int* new_ref_0;

    if (test_idx_0 < num_tests){
        snprintf(test_num, sizeof(test_num), "%d", test_indices[test_idx_0]);
        con[0] = '\0'; 
        strcat(con, file_prefix);
        strcat(con, test_num);
        strcat(con, ".SEQ");

        reads[0] = '\0';
        strcat(reads, file_prefix);
        strcat(reads, test_num);
        strcat(reads, ".READS");

        printf("TARGET %s\n", test_num);

        start = std::chrono::high_resolution_clock::now();
        // Malloc the largest array
        con_arr = (char*) malloc( CON_SIZE * CON_LEN * sizeof(char));
        reads_arr = (char*) malloc( READS_SIZE * READS_LEN * sizeof(char));
        weights_arr = (char*) malloc( READS_SIZE * READS_LEN * sizeof(char));

        parse( con,1, con_arr, &con_len, &con_size, &con_total_len);
        parse( reads, 4, reads_arr, &reads_len, &reads_size_0 , &reads_total_len);
        parse( reads, 5, weights_arr, &reads_len, &reads_size_0, &weights_total_len);

        finish = std::chrono::high_resolution_clock::now();
        //std::chrono::duration<double> elapsed = finish - start;   
        duration = std::chrono::duration_cast<std::chrono::milliseconds>(finish - start);
        parse_time[0] = duration;
        //std::cout << "Parsing time is :" << duration.count() << " ms\n";
        //printf("Parse time is: %0.3f ms\n", duration.count());

        min_whd = (int*) malloc(con_size * reads_size_0 * sizeof(int));
        min_whd_idx = (int*) malloc(con_size * reads_size_0 * sizeof(int));
        new_ref_0 = (int*) malloc(reads_size_0 * sizeof(int));
        //int* new_ref_idx = (int*) malloc(reads_size_0 * sizeof(int));
        new_ref_idx_ref_0 = (int*) malloc(reads_size_0 * sizeof(int));

        whd(con_arr, con_size, con_len, reads_arr, reads_size_0, reads_len, weights_arr, min_whd, min_whd_idx);
        score_whd (min_whd,  min_whd_idx, con_size, reads_size_0, new_ref_0, new_ref_idx_ref_0);

        printf("Software Success!\n");
    } else {
        con_size=reads_size_0=0;
        con_total_len=reads_total_len=0;
        con_arr=reads_arr=weights_arr=NULL;
        con_len=reads_len=NULL;
        min_whd=min_whd_idx=new_ref_0=new_ref_idx_ref_0=NULL;   
    }

    //return 0;
    //Allocate Memory in Host Memory
    //std::vector<char,aligned_allocator<char>> con_arr_buffer     ( con_arr, con_arr + CON_SIZE * CON_LEN);
    //std::vector<char,aligned_allocator<char>> reads_arr_buffer     (reads_arr, reads_arr + READS_SIZE * READS_LEN);
    //std::vector<char,aligned_allocator<char>> weights_arr_buffer     (weights_arr, weights_arr + READS_SIZE * READS_LEN);

    start = std::chrono::high_resolution_clock::now();
    std::vector<ap_uint<4>,aligned_allocator<ap_uint<4>>> con_arr_buffer_0     ( CON_SIZE * CON_LEN);
    //std::vector<ap_uint<4>,aligned_allocator<ap_uint<4>>> con_arr_buffer_0     ( CON_SIZE * CON_LEN >> 1);
    //std::vector<char,aligned_allocator<char>> con_arr_buffer_0      ( con_total_len >> 1);
    std::vector<ap_uint<4>,aligned_allocator<ap_uint<4>>> reads_arr_buffer_0     ( READS_SIZE * READS_LEN);
    //std::vector<ap_uint<4>,aligned_allocator<ap_uint<4>>> reads_arr_buffer_0     ( READS_SIZE * READS_LEN >> 1);
    //std::vector<char,aligned_allocator<char>> reads_arr_buffer_0    ( reads_total_len >> 1);
    std::vector<char,aligned_allocator<char>> weights_arr_buffer_0 (READS_SIZE * READS_LEN);
    //std::vector<char,aligned_allocator<char>> weights_arr_buffer_0  (weights_arr, weights_arr + READS_SIZE * READS_LEN);
    //std::vector<char,aligned_allocator<char>> weights_arr_buffer_0  (weights_arr, weights_arr + reads_total_len);

    for(int i = 0 ; i < con_total_len; i++){
    //for(int i = 0 ; i < CON_SIZE * CON_LEN; i++){
        con_arr_buffer_0[i] = m[con_arr[i]];
    }

    //for(int i = 0 ; i < CON_SIZE * CON_LEN; i++){
   // for(int i = 0 ; i < con_total_len; i++){
   //     int idx = i >> 1;
   //     int offset = (i % 2) << 2;
   //     char c = m[con_arr[i]] & 0xf;
   //     con_arr_buffer_0[idx] = (con_arr_buffer_0[idx] & ~(0xf << offset)) | (c << offset);
   //     //printf("con_arr_buffer_0[%d + %d]: %hhX \n", idx, offset, con_arr_buffer_0[idx]);
   // }

    printf("Read Buffer:");     
    for(int i = 0 ; i < reads_total_len; i++){
    //for(int i = 0 ; i < READS_LEN * READS_SIZE; i++){
        reads_arr_buffer_0[i] = m[reads_arr[i]];
        weights_arr_buffer_0[i] = weights_arr[i];
        //unsigned char print_var = reads_arr_buffer[i];
        //printf("%x", print_var);     
    }
   // for(int i = 0 ; i < reads_total_len; i++){
   //     int idx = i >> 1;
   //     int offset = (i % 2) << 2;
   //     char c = m[reads_arr[i]] & 0xf;
   //     reads_arr_buffer_0[idx] = (reads_arr_buffer_0[idx] & ~(0xf << offset)) | (c << offset);
   //     //printf("reads_arr_buffer_0[%d]: %hhX \n", idx, (char)reads_arr_buffer_0[idx]);
   // }

    //std::vector<int,aligned_allocator<int>> con_len_buffer_0     (con_len, con_len + con_size);
    //std::vector<int,aligned_allocator<int>> con_len_buffer_0     (con_len, con_len + CON_SIZE);
    std::vector<int,aligned_allocator<int>> con_len_buffer_0     (CON_SIZE);
    
    for(int i = 0 ; i < con_size; i++){con_len_buffer_0[i] = con_len[i];}
    //std::vector<int,aligned_allocator<int>> reads_len_buffer_0     (reads_len, reads_len + reads_size_0);
    //std::vector<int,aligned_allocator<int>> reads_len_buffer_0     (reads_len, reads_len + READS_SIZE);
    std::vector<int,aligned_allocator<int>> reads_len_buffer_0     (READS_SIZE);
    for(int i = 0 ; i < reads_size_0; i++){reads_len_buffer_0[i] = reads_len[i];}

    printf("narg: %d\n", narg);
    printf("reads_size_0: %d\n", reads_size_0);
    printf("con_size: %d\n", con_size );
    printf("con_len: \n");
    for(int i =0; i < con_size; i++) {
      printf("%d\t", con_len_buffer_0[i]);
    }
    printf("\n");

    //std::vector<int,aligned_allocator<int>> whd_buffer(con_size * reads_size_0 << 1);
    //std::vector<int,aligned_allocator<int>> new_ref_idx_0(new_ref_idx_0 + reads_size_0);
    //std::vector<int,aligned_allocator<int>> new_ref_idx_0(READS_SIZE);
    std::vector<int,aligned_allocator<int>> * new_ref_idx_0 = new std::vector<int,aligned_allocator<int>> (READS_SIZE); 

    for(int i = 0 ; i < reads_size_0; i++){
        (*new_ref_idx_0)[i] = 0;
    }
    printf("\n");

    finish = std::chrono::high_resolution_clock::now();

    duration = std::chrono::duration_cast<std::chrono::milliseconds>(finish - start);
    std::cout << "Preprocess time is : " << duration.count() << " s\n";
    //printf("Preprocess time is: %0.3f ms\n", duration.count());


    cl_mem_ext_ptr_t con_arr_buffer_ptr_0, reads_arr_buffer_ptr_0, weights_arr_buffer_ptr_0, con_len_buffer_ptr_0, reads_len_buffer_ptr_0, whd_buffer_ptr_0,  new_ref_idx_ptr_0; 
    //ddr_bank = banks[kernel_idx % 4];
    con_arr_buffer_ptr_0.flags  = ddr_bank; 
    con_len_buffer_ptr_0.flags  = ddr_bank; 
    reads_arr_buffer_ptr_0.flags  = ddr_bank; 
    reads_len_buffer_ptr_0.flags  = ddr_bank; 
    weights_arr_buffer_ptr_0.flags  = ddr_bank; 
    new_ref_idx_ptr_0.flags  = ddr_bank; 
    //whd_buffer_ptr_0.flags = ddr_bank;
 
    std::cout << "Preprocess time is : " << duration.count() << " ms\n";

    // Setting input and output objects
    con_arr_buffer_ptr_0.obj = con_arr_buffer_0.data();
    //con_len_buffer_ptr_0.obj = con_len;  // this would cause extra memcpy since it is not aligned
    con_len_buffer_ptr_0.obj = con_len_buffer_0.data();
    reads_arr_buffer_ptr_0.obj = reads_arr_buffer_0.data();
    //reads_len_buffer_ptr_0.obj = reads_len;
    reads_len_buffer_ptr_0.obj = reads_len_buffer_0.data();
    weights_arr_buffer_ptr_0.obj = weights_arr_buffer_0.data();
    //whd_buffer_ptr_0.obj = whd_buffer_0.data();
    new_ref_idx_ptr_0.obj = new_ref_idx_0->data();

    std::cout << "Preprocess time is : " << duration.count() << " ms\n";
    // Setting param to zero 
    con_arr_buffer_ptr_0.param = 0; con_len_buffer_ptr_0.param = 0;
    reads_arr_buffer_ptr_0.param = 0; reads_len_buffer_ptr_0.param = 0; 
    weights_arr_buffer_ptr_0.param = 0; new_ref_idx_ptr_0.param = 0;
    whd_buffer_ptr_0.param = 0;

    printf("Finish creating buffers\n");
    //Allocate Buffer in Global Memory

    cl::Buffer con_arr_input_0(context,CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY | CL_MEM_EXT_PTR_XILINX ,\ 
            CON_SIZE * CON_LEN, &con_arr_buffer_ptr_0);
    cl::Buffer reads_arr_input_0(context,CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY | CL_MEM_EXT_PTR_XILINX ,\ 
            READS_SIZE * READS_LEN, &reads_arr_buffer_ptr_0);
    cl::Buffer weights_arr_input_0(context,CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY | CL_MEM_EXT_PTR_XILINX ,\ 
            READS_SIZE * READS_LEN, &weights_arr_buffer_ptr_0);

    cl::Buffer con_len_input_0(context,CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY | CL_MEM_EXT_PTR_XILINX , 
            CON_SIZE * sizeof(int), &con_len_buffer_ptr_0);
    cl::Buffer reads_len_input_0(context,CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY | CL_MEM_EXT_PTR_XILINX , 
            READS_SIZE * sizeof(int), &reads_len_buffer_ptr_0);

    //cl::Buffer whd_output_0(context, CL_MEM_USE_HOST_PTR | CL_MEM_WRITE_ONLY| CL_MEM_EXT_PTR_XILINX ,
    //        reads_size_0 * con_size * 2 * sizeof(int), &whd_buffer_ptr_0);
    cl::Buffer new_ref_idx_output_0(context, CL_MEM_USE_HOST_PTR | CL_MEM_WRITE_ONLY | CL_MEM_EXT_PTR_XILINX ,
            READS_SIZE * sizeof(int), &new_ref_idx_ptr_0);

    inBufVec.push_back(con_len_input_0);
    inBufVec.push_back(reads_len_input_0);
    inBufVec.push_back(con_arr_input_0);
    inBufVec.push_back(reads_arr_input_0);
    inBufVec.push_back(weights_arr_input_0);
    //outBufVec.push_back(whd_output_0);
    outBufVec.push_back(new_ref_idx_output_0);

    //Set the Kernel Arguments
    krnl_indel.setArg(narg++, con_arr_input_0);
    krnl_indel.setArg(narg++, con_size);
    krnl_indel.setArg(narg++, con_len_input_0);
    krnl_indel.setArg(narg++, reads_arr_input_0);
    krnl_indel.setArg(narg++, reads_size_0);
    krnl_indel.setArg(narg++, reads_len_input_0);
    krnl_indel.setArg(narg++, weights_arr_input_0);
    //krnl_indel.setArg(narg++, whd_output_0);
    krnl_indel.setArg(narg++, new_ref_idx_output_0);
   

//INST_1============================================//
    int test_idx_1= test_idx + 1;

    printf("test_idx_1: %d\n", test_idx_1);
    //int reads_size_1;
    int* new_ref_1;

    if (test_idx_1 < num_tests){
        snprintf(test_num, sizeof(test_num), "%d", test_indices[test_idx_1]);
        con[0] = '\0'; 
        strcat(con, file_prefix);
        strcat(con, test_num);
        strcat(con, ".SEQ");

        reads[0] = '\0';
        strcat(reads, file_prefix);
        strcat(reads, test_num);
        strcat(reads, ".READS");

        printf("TARGET %s\n", test_num);

        start = std::chrono::high_resolution_clock::now();
        // Malloc the largest array
        con_arr = (char*) malloc( CON_SIZE * CON_LEN * sizeof(char));
        reads_arr = (char*) malloc( READS_SIZE * READS_LEN * sizeof(char));
        weights_arr = (char*) malloc( READS_SIZE * READS_LEN * sizeof(char));

        parse( con,1, con_arr, &con_len, &con_size, &con_total_len);
        parse( reads, 4, reads_arr, &reads_len, &reads_size_1 , &reads_total_len);
        parse( reads, 5, weights_arr, &reads_len, &reads_size_1, &weights_total_len);

        finish = std::chrono::high_resolution_clock::now();
        //std::chrono::duration<double> elapsed = finish - start;   
        duration = std::chrono::duration_cast<std::chrono::milliseconds>(finish - start);
        parse_time[1] = duration;
        //std::cout << "Parsing time is :" << duration.count() << " ms\n";
        //printf("Parse time is: %0.3f ms\n", duration.count());

        min_whd = (int*) malloc(con_size * reads_size_1 * sizeof(int));
        min_whd_idx = (int*) malloc(con_size * reads_size_1 * sizeof(int));
        new_ref_1 = (int*) malloc(reads_size_1 * sizeof(int));
        //int* new_ref_idx = (int*) malloc(reads_size_1 * sizeof(int));
        new_ref_idx_ref_1 = (int*) malloc(reads_size_1 * sizeof(int));

        whd(con_arr, con_size, con_len, reads_arr, reads_size_1, reads_len, weights_arr, min_whd, min_whd_idx);
        score_whd (min_whd,  min_whd_idx, con_size, reads_size_1, new_ref_1, new_ref_idx_ref_1);

        printf("Software Success!\n");
    } else {
        con_size=reads_size_1=0;
        con_total_len=reads_total_len=0;
        con_arr=reads_arr=weights_arr=NULL;
        con_len=reads_len=NULL;
        min_whd=min_whd_idx=new_ref_1=new_ref_idx_ref_1=NULL;   
    }

    //return 0;
    //Allocate Memory in Host Memory
    //std::vector<char,aligned_allocator<char>> con_arr_buffer     ( con_arr, con_arr + CON_SIZE * CON_LEN);
    //std::vector<char,aligned_allocator<char>> reads_arr_buffer     (reads_arr, reads_arr + READS_SIZE * READS_LEN);
    //std::vector<char,aligned_allocator<char>> weights_arr_buffer     (weights_arr, weights_arr + READS_SIZE * READS_LEN);

    start = std::chrono::high_resolution_clock::now();
    std::vector<ap_uint<4>,aligned_allocator<ap_uint<4>>> con_arr_buffer_1     ( CON_SIZE * CON_LEN);
    //std::vector<ap_uint<4>,aligned_allocator<ap_uint<4>>> con_arr_buffer_1     ( CON_SIZE * CON_LEN >> 1);
    //std::vector<char,aligned_allocator<char>> con_arr_buffer_1      ( con_total_len >> 1);
    std::vector<ap_uint<4>,aligned_allocator<ap_uint<4>>> reads_arr_buffer_1     ( READS_SIZE * READS_LEN);
    //std::vector<ap_uint<4>,aligned_allocator<ap_uint<4>>> reads_arr_buffer_1     ( READS_SIZE * READS_LEN >> 1);
    //std::vector<char,aligned_allocator<char>> reads_arr_buffer_1    ( reads_total_len >> 1);
    std::vector<char,aligned_allocator<char>> weights_arr_buffer_1 (READS_SIZE * READS_LEN);
    //std::vector<char,aligned_allocator<char>> weights_arr_buffer_1  (weights_arr, weights_arr + READS_SIZE * READS_LEN);
    //std::vector<char,aligned_allocator<char>> weights_arr_buffer_1  (weights_arr, weights_arr + reads_total_len);

    for(int i = 0 ; i < con_total_len; i++){
    //for(int i = 0 ; i < CON_SIZE * CON_LEN; i++){
        con_arr_buffer_1[i] = m[con_arr[i]];
    }

    //for(int i = 0 ; i < CON_SIZE * CON_LEN; i++){
   // for(int i = 0 ; i < con_total_len; i++){
   //     int idx = i >> 1;
   //     int offset = (i % 2) << 2;
   //     char c = m[con_arr[i]] & 0xf;
   //     con_arr_buffer_1[idx] = (con_arr_buffer_1[idx] & ~(0xf << offset)) | (c << offset);
   //     //printf("con_arr_buffer_1[%d + %d]: %hhX \n", idx, offset, con_arr_buffer_1[idx]);
   // }

    printf("Read Buffer:");     
    for(int i = 0 ; i < reads_total_len; i++){
    //for(int i = 0 ; i < READS_LEN * READS_SIZE; i++){
        reads_arr_buffer_1[i] = m[reads_arr[i]];
        weights_arr_buffer_1[i] = weights_arr[i];
        //unsigned char print_var = reads_arr_buffer[i];
        //printf("%x", print_var);     
    }
   // for(int i = 0 ; i < reads_total_len; i++){
   //     int idx = i >> 1;
   //     int offset = (i % 2) << 2;
   //     char c = m[reads_arr[i]] & 0xf;
   //     reads_arr_buffer_1[idx] = (reads_arr_buffer_1[idx] & ~(0xf << offset)) | (c << offset);
   //     //printf("reads_arr_buffer_1[%d]: %hhX \n", idx, (char)reads_arr_buffer_1[idx]);
   // }

    //std::vector<int,aligned_allocator<int>> con_len_buffer_1     (con_len, con_len + con_size);
    //std::vector<int,aligned_allocator<int>> con_len_buffer_1     (con_len, con_len + CON_SIZE);
    std::vector<int,aligned_allocator<int>> con_len_buffer_1     (CON_SIZE);
    
    for(int i = 0 ; i < con_size; i++){con_len_buffer_1[i] = con_len[i];}
    //std::vector<int,aligned_allocator<int>> reads_len_buffer_1     (reads_len, reads_len + reads_size_1);
    //std::vector<int,aligned_allocator<int>> reads_len_buffer_1     (reads_len, reads_len + READS_SIZE);
    std::vector<int,aligned_allocator<int>> reads_len_buffer_1     (READS_SIZE);
    for(int i = 0 ; i < reads_size_1; i++){reads_len_buffer_1[i] = reads_len[i];}

    printf("narg: %d\n", narg);
    printf("reads_size_1: %d\n", reads_size_1);
    printf("con_size: %d\n", con_size );
    printf("con_len: \n");
    for(int i =0; i < con_size; i++) {
      printf("%d\t", con_len_buffer_1[i]);
    }
    printf("\n");

    //std::vector<int,aligned_allocator<int>> whd_buffer(con_size * reads_size_1 << 1);
    //std::vector<int,aligned_allocator<int>> new_ref_idx_1(new_ref_idx_1 + reads_size_1);
    //std::vector<int,aligned_allocator<int>> new_ref_idx_1(READS_SIZE);
    std::vector<int,aligned_allocator<int>> * new_ref_idx_1 = new std::vector<int,aligned_allocator<int>> (READS_SIZE); 

    for(int i = 0 ; i < reads_size_1; i++){
        (*new_ref_idx_1)[i] = 0;
    }
    printf("\n");

    finish = std::chrono::high_resolution_clock::now();

    duration = std::chrono::duration_cast<std::chrono::milliseconds>(finish - start);
    std::cout << "Preprocess time is : " << duration.count() << " s\n";
    //printf("Preprocess time is: %0.3f ms\n", duration.count());


    cl_mem_ext_ptr_t con_arr_buffer_ptr_1, reads_arr_buffer_ptr_1, weights_arr_buffer_ptr_1, con_len_buffer_ptr_1, reads_len_buffer_ptr_1, whd_buffer_ptr_1,  new_ref_idx_ptr_1; 
    //ddr_bank = banks[kernel_idx % 4];
    con_arr_buffer_ptr_1.flags  = ddr_bank; 
    con_len_buffer_ptr_1.flags  = ddr_bank; 
    reads_arr_buffer_ptr_1.flags  = ddr_bank; 
    reads_len_buffer_ptr_1.flags  = ddr_bank; 
    weights_arr_buffer_ptr_1.flags  = ddr_bank; 
    new_ref_idx_ptr_1.flags  = ddr_bank; 
    //whd_buffer_ptr_1.flags = ddr_bank;
 
    std::cout << "Preprocess time is : " << duration.count() << " ms\n";

    // Setting input and output objects
    con_arr_buffer_ptr_1.obj = con_arr_buffer_1.data();
    //con_len_buffer_ptr_1.obj = con_len;  // this would cause extra memcpy since it is not aligned
    con_len_buffer_ptr_1.obj = con_len_buffer_1.data();
    reads_arr_buffer_ptr_1.obj = reads_arr_buffer_1.data();
    //reads_len_buffer_ptr_1.obj = reads_len;
    reads_len_buffer_ptr_1.obj = reads_len_buffer_1.data();
    weights_arr_buffer_ptr_1.obj = weights_arr_buffer_1.data();
    //whd_buffer_ptr_1.obj = whd_buffer_1.data();
    new_ref_idx_ptr_1.obj = new_ref_idx_1->data();

    std::cout << "Preprocess time is : " << duration.count() << " ms\n";
    // Setting param to zero 
    con_arr_buffer_ptr_1.param = 0; con_len_buffer_ptr_1.param = 0;
    reads_arr_buffer_ptr_1.param = 0; reads_len_buffer_ptr_1.param = 0; 
    weights_arr_buffer_ptr_1.param = 0; new_ref_idx_ptr_1.param = 0;
    whd_buffer_ptr_1.param = 0;

    printf("Finish creating buffers\n");
    //Allocate Buffer in Global Memory

    cl::Buffer con_arr_input_1(context,CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY | CL_MEM_EXT_PTR_XILINX ,\ 
            CON_SIZE * CON_LEN, &con_arr_buffer_ptr_1);
    cl::Buffer reads_arr_input_1(context,CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY | CL_MEM_EXT_PTR_XILINX ,\ 
            READS_SIZE * READS_LEN, &reads_arr_buffer_ptr_1);
    cl::Buffer weights_arr_input_1(context,CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY | CL_MEM_EXT_PTR_XILINX ,\ 
            READS_SIZE * READS_LEN, &weights_arr_buffer_ptr_1);

    cl::Buffer con_len_input_1(context,CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY | CL_MEM_EXT_PTR_XILINX , 
            CON_SIZE * sizeof(int), &con_len_buffer_ptr_1);
    cl::Buffer reads_len_input_1(context,CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY | CL_MEM_EXT_PTR_XILINX , 
            READS_SIZE * sizeof(int), &reads_len_buffer_ptr_1);

    //cl::Buffer whd_output_1(context, CL_MEM_USE_HOST_PTR | CL_MEM_WRITE_ONLY| CL_MEM_EXT_PTR_XILINX ,
    //        reads_size_1 * con_size * 2 * sizeof(int), &whd_buffer_ptr_1);
    cl::Buffer new_ref_idx_output_1(context, CL_MEM_USE_HOST_PTR | CL_MEM_WRITE_ONLY | CL_MEM_EXT_PTR_XILINX ,
            READS_SIZE * sizeof(int), &new_ref_idx_ptr_1);

    inBufVec.push_back(con_len_input_1);
    inBufVec.push_back(reads_len_input_1);
    inBufVec.push_back(con_arr_input_1);
    inBufVec.push_back(reads_arr_input_1);
    inBufVec.push_back(weights_arr_input_1);
    //outBufVec.push_back(whd_output_1);
    outBufVec.push_back(new_ref_idx_output_1);

    //Set the Kernel Arguments
    krnl_indel.setArg(narg++, con_arr_input_1);
    krnl_indel.setArg(narg++, con_size);
    krnl_indel.setArg(narg++, con_len_input_1);
    krnl_indel.setArg(narg++, reads_arr_input_1);
    krnl_indel.setArg(narg++, reads_size_1);
    krnl_indel.setArg(narg++, reads_len_input_1);
    krnl_indel.setArg(narg++, weights_arr_input_1);
    //krnl_indel.setArg(narg++, whd_output_1);
    krnl_indel.setArg(narg++, new_ref_idx_output_1);
   

//INST_2============================================//
    int test_idx_2= test_idx + 2;

    printf("test_idx_2: %d\n", test_idx_2);
    //int reads_size_2;
    int* new_ref_2;

    if (test_idx_2 < num_tests){
        snprintf(test_num, sizeof(test_num), "%d", test_indices[test_idx_2]);
        con[0] = '\0'; 
        strcat(con, file_prefix);
        strcat(con, test_num);
        strcat(con, ".SEQ");

        reads[0] = '\0';
        strcat(reads, file_prefix);
        strcat(reads, test_num);
        strcat(reads, ".READS");

        printf("TARGET %s\n", test_num);

        start = std::chrono::high_resolution_clock::now();
        // Malloc the largest array
        con_arr = (char*) malloc( CON_SIZE * CON_LEN * sizeof(char));
        reads_arr = (char*) malloc( READS_SIZE * READS_LEN * sizeof(char));
        weights_arr = (char*) malloc( READS_SIZE * READS_LEN * sizeof(char));

        parse( con,1, con_arr, &con_len, &con_size, &con_total_len);
        parse( reads, 4, reads_arr, &reads_len, &reads_size_2 , &reads_total_len);
        parse( reads, 5, weights_arr, &reads_len, &reads_size_2, &weights_total_len);

        finish = std::chrono::high_resolution_clock::now();
        //std::chrono::duration<double> elapsed = finish - start;   
        duration = std::chrono::duration_cast<std::chrono::milliseconds>(finish - start);
        parse_time[2] = duration;
        //std::cout << "Parsing time is :" << duration.count() << " ms\n";
        //printf("Parse time is: %0.3f ms\n", duration.count());

        min_whd = (int*) malloc(con_size * reads_size_2 * sizeof(int));
        min_whd_idx = (int*) malloc(con_size * reads_size_2 * sizeof(int));
        new_ref_2 = (int*) malloc(reads_size_2 * sizeof(int));
        //int* new_ref_idx = (int*) malloc(reads_size_2 * sizeof(int));
        new_ref_idx_ref_2 = (int*) malloc(reads_size_2 * sizeof(int));

        whd(con_arr, con_size, con_len, reads_arr, reads_size_2, reads_len, weights_arr, min_whd, min_whd_idx);
        score_whd (min_whd,  min_whd_idx, con_size, reads_size_2, new_ref_2, new_ref_idx_ref_2);

        printf("Software Success!\n");
    } else {
        con_size=reads_size_2=0;
        con_total_len=reads_total_len=0;
        con_arr=reads_arr=weights_arr=NULL;
        con_len=reads_len=NULL;
        min_whd=min_whd_idx=new_ref_2=new_ref_idx_ref_2=NULL;   
    }

    //return 0;
    //Allocate Memory in Host Memory
    //std::vector<char,aligned_allocator<char>> con_arr_buffer     ( con_arr, con_arr + CON_SIZE * CON_LEN);
    //std::vector<char,aligned_allocator<char>> reads_arr_buffer     (reads_arr, reads_arr + READS_SIZE * READS_LEN);
    //std::vector<char,aligned_allocator<char>> weights_arr_buffer     (weights_arr, weights_arr + READS_SIZE * READS_LEN);

    start = std::chrono::high_resolution_clock::now();
    std::vector<ap_uint<4>,aligned_allocator<ap_uint<4>>> con_arr_buffer_2     ( CON_SIZE * CON_LEN);
    //std::vector<ap_uint<4>,aligned_allocator<ap_uint<4>>> con_arr_buffer_2     ( CON_SIZE * CON_LEN >> 1);
    //std::vector<char,aligned_allocator<char>> con_arr_buffer_2      ( con_total_len >> 1);
    std::vector<ap_uint<4>,aligned_allocator<ap_uint<4>>> reads_arr_buffer_2     ( READS_SIZE * READS_LEN);
    //std::vector<ap_uint<4>,aligned_allocator<ap_uint<4>>> reads_arr_buffer_2     ( READS_SIZE * READS_LEN >> 1);
    //std::vector<char,aligned_allocator<char>> reads_arr_buffer_2    ( reads_total_len >> 1);
    std::vector<char,aligned_allocator<char>> weights_arr_buffer_2 (READS_SIZE * READS_LEN);
    //std::vector<char,aligned_allocator<char>> weights_arr_buffer_2  (weights_arr, weights_arr + READS_SIZE * READS_LEN);
    //std::vector<char,aligned_allocator<char>> weights_arr_buffer_2  (weights_arr, weights_arr + reads_total_len);

    for(int i = 0 ; i < con_total_len; i++){
    //for(int i = 0 ; i < CON_SIZE * CON_LEN; i++){
        con_arr_buffer_2[i] = m[con_arr[i]];
    }

    //for(int i = 0 ; i < CON_SIZE * CON_LEN; i++){
   // for(int i = 0 ; i < con_total_len; i++){
   //     int idx = i >> 1;
   //     int offset = (i % 2) << 2;
   //     char c = m[con_arr[i]] & 0xf;
   //     con_arr_buffer_2[idx] = (con_arr_buffer_2[idx] & ~(0xf << offset)) | (c << offset);
   //     //printf("con_arr_buffer_2[%d + %d]: %hhX \n", idx, offset, con_arr_buffer_2[idx]);
   // }

    printf("Read Buffer:");     
    for(int i = 0 ; i < reads_total_len; i++){
    //for(int i = 0 ; i < READS_LEN * READS_SIZE; i++){
        reads_arr_buffer_2[i] = m[reads_arr[i]];
        weights_arr_buffer_2[i] = weights_arr[i];
        //unsigned char print_var = reads_arr_buffer[i];
        //printf("%x", print_var);     
    }
   // for(int i = 0 ; i < reads_total_len; i++){
   //     int idx = i >> 1;
   //     int offset = (i % 2) << 2;
   //     char c = m[reads_arr[i]] & 0xf;
   //     reads_arr_buffer_2[idx] = (reads_arr_buffer_2[idx] & ~(0xf << offset)) | (c << offset);
   //     //printf("reads_arr_buffer_2[%d]: %hhX \n", idx, (char)reads_arr_buffer_2[idx]);
   // }

    //std::vector<int,aligned_allocator<int>> con_len_buffer_2     (con_len, con_len + con_size);
    //std::vector<int,aligned_allocator<int>> con_len_buffer_2     (con_len, con_len + CON_SIZE);
    std::vector<int,aligned_allocator<int>> con_len_buffer_2     (CON_SIZE);
    
    for(int i = 0 ; i < con_size; i++){con_len_buffer_2[i] = con_len[i];}
    //std::vector<int,aligned_allocator<int>> reads_len_buffer_2     (reads_len, reads_len + reads_size_2);
    //std::vector<int,aligned_allocator<int>> reads_len_buffer_2     (reads_len, reads_len + READS_SIZE);
    std::vector<int,aligned_allocator<int>> reads_len_buffer_2     (READS_SIZE);
    for(int i = 0 ; i < reads_size_2; i++){reads_len_buffer_2[i] = reads_len[i];}

    printf("narg: %d\n", narg);
    printf("reads_size_2: %d\n", reads_size_2);
    printf("con_size: %d\n", con_size );
    printf("con_len: \n");
    for(int i =0; i < con_size; i++) {
      printf("%d\t", con_len_buffer_2[i]);
    }
    printf("\n");

    //std::vector<int,aligned_allocator<int>> whd_buffer(con_size * reads_size_2 << 1);
    //std::vector<int,aligned_allocator<int>> new_ref_idx_2(new_ref_idx_2 + reads_size_2);
    //std::vector<int,aligned_allocator<int>> new_ref_idx_2(READS_SIZE);
    std::vector<int,aligned_allocator<int>> * new_ref_idx_2 = new std::vector<int,aligned_allocator<int>> (READS_SIZE); 

    for(int i = 0 ; i < reads_size_2; i++){
        (*new_ref_idx_2)[i] = 0;
    }
    printf("\n");

    finish = std::chrono::high_resolution_clock::now();

    duration = std::chrono::duration_cast<std::chrono::milliseconds>(finish - start);
    std::cout << "Preprocess time is : " << duration.count() << " s\n";
    //printf("Preprocess time is: %0.3f ms\n", duration.count());


    cl_mem_ext_ptr_t con_arr_buffer_ptr_2, reads_arr_buffer_ptr_2, weights_arr_buffer_ptr_2, con_len_buffer_ptr_2, reads_len_buffer_ptr_2, whd_buffer_ptr_2,  new_ref_idx_ptr_2; 
    //ddr_bank = banks[kernel_idx % 4];
    con_arr_buffer_ptr_2.flags  = ddr_bank; 
    con_len_buffer_ptr_2.flags  = ddr_bank; 
    reads_arr_buffer_ptr_2.flags  = ddr_bank; 
    reads_len_buffer_ptr_2.flags  = ddr_bank; 
    weights_arr_buffer_ptr_2.flags  = ddr_bank; 
    new_ref_idx_ptr_2.flags  = ddr_bank; 
    //whd_buffer_ptr_2.flags = ddr_bank;
 
    std::cout << "Preprocess time is : " << duration.count() << " ms\n";

    // Setting input and output objects
    con_arr_buffer_ptr_2.obj = con_arr_buffer_2.data();
    //con_len_buffer_ptr_2.obj = con_len;  // this would cause extra memcpy since it is not aligned
    con_len_buffer_ptr_2.obj = con_len_buffer_2.data();
    reads_arr_buffer_ptr_2.obj = reads_arr_buffer_2.data();
    //reads_len_buffer_ptr_2.obj = reads_len;
    reads_len_buffer_ptr_2.obj = reads_len_buffer_2.data();
    weights_arr_buffer_ptr_2.obj = weights_arr_buffer_2.data();
    //whd_buffer_ptr_2.obj = whd_buffer_2.data();
    new_ref_idx_ptr_2.obj = new_ref_idx_2->data();

    std::cout << "Preprocess time is : " << duration.count() << " ms\n";
    // Setting param to zero 
    con_arr_buffer_ptr_2.param = 0; con_len_buffer_ptr_2.param = 0;
    reads_arr_buffer_ptr_2.param = 0; reads_len_buffer_ptr_2.param = 0; 
    weights_arr_buffer_ptr_2.param = 0; new_ref_idx_ptr_2.param = 0;
    whd_buffer_ptr_2.param = 0;

    printf("Finish creating buffers\n");
    //Allocate Buffer in Global Memory

    cl::Buffer con_arr_input_2(context,CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY | CL_MEM_EXT_PTR_XILINX ,\ 
            CON_SIZE * CON_LEN, &con_arr_buffer_ptr_2);
    cl::Buffer reads_arr_input_2(context,CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY | CL_MEM_EXT_PTR_XILINX ,\ 
            READS_SIZE * READS_LEN, &reads_arr_buffer_ptr_2);
    cl::Buffer weights_arr_input_2(context,CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY | CL_MEM_EXT_PTR_XILINX ,\ 
            READS_SIZE * READS_LEN, &weights_arr_buffer_ptr_2);

    cl::Buffer con_len_input_2(context,CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY | CL_MEM_EXT_PTR_XILINX , 
            CON_SIZE * sizeof(int), &con_len_buffer_ptr_2);
    cl::Buffer reads_len_input_2(context,CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY | CL_MEM_EXT_PTR_XILINX , 
            READS_SIZE * sizeof(int), &reads_len_buffer_ptr_2);

    //cl::Buffer whd_output_2(context, CL_MEM_USE_HOST_PTR | CL_MEM_WRITE_ONLY| CL_MEM_EXT_PTR_XILINX ,
    //        reads_size_2 * con_size * 2 * sizeof(int), &whd_buffer_ptr_2);
    cl::Buffer new_ref_idx_output_2(context, CL_MEM_USE_HOST_PTR | CL_MEM_WRITE_ONLY | CL_MEM_EXT_PTR_XILINX ,
            READS_SIZE * sizeof(int), &new_ref_idx_ptr_2);

    inBufVec.push_back(con_len_input_2);
    inBufVec.push_back(reads_len_input_2);
    inBufVec.push_back(con_arr_input_2);
    inBufVec.push_back(reads_arr_input_2);
    inBufVec.push_back(weights_arr_input_2);
    //outBufVec.push_back(whd_output_2);
    outBufVec.push_back(new_ref_idx_output_2);

    //Set the Kernel Arguments
    krnl_indel.setArg(narg++, con_arr_input_2);
    krnl_indel.setArg(narg++, con_size);
    krnl_indel.setArg(narg++, con_len_input_2);
    krnl_indel.setArg(narg++, reads_arr_input_2);
    krnl_indel.setArg(narg++, reads_size_2);
    krnl_indel.setArg(narg++, reads_len_input_2);
    krnl_indel.setArg(narg++, weights_arr_input_2);
    //krnl_indel.setArg(narg++, whd_output_2);
    krnl_indel.setArg(narg++, new_ref_idx_output_2);
   

//INST_3============================================//
    int test_idx_3= test_idx + 3;

    printf("test_idx_3: %d\n", test_idx_3);
    //int reads_size_3;
    int* new_ref_3;

    if (test_idx_3 < num_tests){
        snprintf(test_num, sizeof(test_num), "%d", test_indices[test_idx_3]);
        con[0] = '\0'; 
        strcat(con, file_prefix);
        strcat(con, test_num);
        strcat(con, ".SEQ");

        reads[0] = '\0';
        strcat(reads, file_prefix);
        strcat(reads, test_num);
        strcat(reads, ".READS");

        printf("TARGET %s\n", test_num);

        start = std::chrono::high_resolution_clock::now();
        // Malloc the largest array
        con_arr = (char*) malloc( CON_SIZE * CON_LEN * sizeof(char));
        reads_arr = (char*) malloc( READS_SIZE * READS_LEN * sizeof(char));
        weights_arr = (char*) malloc( READS_SIZE * READS_LEN * sizeof(char));

        parse( con,1, con_arr, &con_len, &con_size, &con_total_len);
        parse( reads, 4, reads_arr, &reads_len, &reads_size_3 , &reads_total_len);
        parse( reads, 5, weights_arr, &reads_len, &reads_size_3, &weights_total_len);

        finish = std::chrono::high_resolution_clock::now();
        //std::chrono::duration<double> elapsed = finish - start;   
        duration = std::chrono::duration_cast<std::chrono::milliseconds>(finish - start);
        parse_time[3] = duration;
        //std::cout << "Parsing time is :" << duration.count() << " ms\n";
        //printf("Parse time is: %0.3f ms\n", duration.count());

        min_whd = (int*) malloc(con_size * reads_size_3 * sizeof(int));
        min_whd_idx = (int*) malloc(con_size * reads_size_3 * sizeof(int));
        new_ref_3 = (int*) malloc(reads_size_3 * sizeof(int));
        //int* new_ref_idx = (int*) malloc(reads_size_3 * sizeof(int));
        new_ref_idx_ref_3 = (int*) malloc(reads_size_3 * sizeof(int));

        whd(con_arr, con_size, con_len, reads_arr, reads_size_3, reads_len, weights_arr, min_whd, min_whd_idx);
        score_whd (min_whd,  min_whd_idx, con_size, reads_size_3, new_ref_3, new_ref_idx_ref_3);

        printf("Software Success!\n");
    } else {
        con_size=reads_size_3=0;
        con_total_len=reads_total_len=0;
        con_arr=reads_arr=weights_arr=NULL;
        con_len=reads_len=NULL;
        min_whd=min_whd_idx=new_ref_3=new_ref_idx_ref_3=NULL;   
    }

    //return 0;
    //Allocate Memory in Host Memory
    //std::vector<char,aligned_allocator<char>> con_arr_buffer     ( con_arr, con_arr + CON_SIZE * CON_LEN);
    //std::vector<char,aligned_allocator<char>> reads_arr_buffer     (reads_arr, reads_arr + READS_SIZE * READS_LEN);
    //std::vector<char,aligned_allocator<char>> weights_arr_buffer     (weights_arr, weights_arr + READS_SIZE * READS_LEN);

    start = std::chrono::high_resolution_clock::now();
    std::vector<ap_uint<4>,aligned_allocator<ap_uint<4>>> con_arr_buffer_3     ( CON_SIZE * CON_LEN);
    //std::vector<ap_uint<4>,aligned_allocator<ap_uint<4>>> con_arr_buffer_3     ( CON_SIZE * CON_LEN >> 1);
    //std::vector<char,aligned_allocator<char>> con_arr_buffer_3      ( con_total_len >> 1);
    std::vector<ap_uint<4>,aligned_allocator<ap_uint<4>>> reads_arr_buffer_3     ( READS_SIZE * READS_LEN);
    //std::vector<ap_uint<4>,aligned_allocator<ap_uint<4>>> reads_arr_buffer_3     ( READS_SIZE * READS_LEN >> 1);
    //std::vector<char,aligned_allocator<char>> reads_arr_buffer_3    ( reads_total_len >> 1);
    std::vector<char,aligned_allocator<char>> weights_arr_buffer_3 (READS_SIZE * READS_LEN);
    //std::vector<char,aligned_allocator<char>> weights_arr_buffer_3  (weights_arr, weights_arr + READS_SIZE * READS_LEN);
    //std::vector<char,aligned_allocator<char>> weights_arr_buffer_3  (weights_arr, weights_arr + reads_total_len);

    for(int i = 0 ; i < con_total_len; i++){
    //for(int i = 0 ; i < CON_SIZE * CON_LEN; i++){
        con_arr_buffer_3[i] = m[con_arr[i]];
    }

    //for(int i = 0 ; i < CON_SIZE * CON_LEN; i++){
   // for(int i = 0 ; i < con_total_len; i++){
   //     int idx = i >> 1;
   //     int offset = (i % 2) << 2;
   //     char c = m[con_arr[i]] & 0xf;
   //     con_arr_buffer_3[idx] = (con_arr_buffer_3[idx] & ~(0xf << offset)) | (c << offset);
   //     //printf("con_arr_buffer_3[%d + %d]: %hhX \n", idx, offset, con_arr_buffer_3[idx]);
   // }

    printf("Read Buffer:");     
    for(int i = 0 ; i < reads_total_len; i++){
    //for(int i = 0 ; i < READS_LEN * READS_SIZE; i++){
        reads_arr_buffer_3[i] = m[reads_arr[i]];
        weights_arr_buffer_3[i] = weights_arr[i];
        //unsigned char print_var = reads_arr_buffer[i];
        //printf("%x", print_var);     
    }
   // for(int i = 0 ; i < reads_total_len; i++){
   //     int idx = i >> 1;
   //     int offset = (i % 2) << 2;
   //     char c = m[reads_arr[i]] & 0xf;
   //     reads_arr_buffer_3[idx] = (reads_arr_buffer_3[idx] & ~(0xf << offset)) | (c << offset);
   //     //printf("reads_arr_buffer_3[%d]: %hhX \n", idx, (char)reads_arr_buffer_3[idx]);
   // }

    //std::vector<int,aligned_allocator<int>> con_len_buffer_3     (con_len, con_len + con_size);
    //std::vector<int,aligned_allocator<int>> con_len_buffer_3     (con_len, con_len + CON_SIZE);
    std::vector<int,aligned_allocator<int>> con_len_buffer_3     (CON_SIZE);
    
    for(int i = 0 ; i < con_size; i++){con_len_buffer_3[i] = con_len[i];}
    //std::vector<int,aligned_allocator<int>> reads_len_buffer_3     (reads_len, reads_len + reads_size_3);
    //std::vector<int,aligned_allocator<int>> reads_len_buffer_3     (reads_len, reads_len + READS_SIZE);
    std::vector<int,aligned_allocator<int>> reads_len_buffer_3     (READS_SIZE);
    for(int i = 0 ; i < reads_size_3; i++){reads_len_buffer_3[i] = reads_len[i];}

    printf("narg: %d\n", narg);
    printf("reads_size_3: %d\n", reads_size_3);
    printf("con_size: %d\n", con_size );
    printf("con_len: \n");
    for(int i =0; i < con_size; i++) {
      printf("%d\t", con_len_buffer_3[i]);
    }
    printf("\n");

    //std::vector<int,aligned_allocator<int>> whd_buffer(con_size * reads_size_3 << 1);
    //std::vector<int,aligned_allocator<int>> new_ref_idx_3(new_ref_idx_3 + reads_size_3);
    //std::vector<int,aligned_allocator<int>> new_ref_idx_3(READS_SIZE);
    std::vector<int,aligned_allocator<int>> * new_ref_idx_3 = new std::vector<int,aligned_allocator<int>> (READS_SIZE); 

    for(int i = 0 ; i < reads_size_3; i++){
        (*new_ref_idx_3)[i] = 0;
    }
    printf("\n");

    finish = std::chrono::high_resolution_clock::now();

    duration = std::chrono::duration_cast<std::chrono::milliseconds>(finish - start);
    std::cout << "Preprocess time is : " << duration.count() << " s\n";
    //printf("Preprocess time is: %0.3f ms\n", duration.count());


    cl_mem_ext_ptr_t con_arr_buffer_ptr_3, reads_arr_buffer_ptr_3, weights_arr_buffer_ptr_3, con_len_buffer_ptr_3, reads_len_buffer_ptr_3, whd_buffer_ptr_3,  new_ref_idx_ptr_3; 
    //ddr_bank = banks[kernel_idx % 4];
    con_arr_buffer_ptr_3.flags  = ddr_bank; 
    con_len_buffer_ptr_3.flags  = ddr_bank; 
    reads_arr_buffer_ptr_3.flags  = ddr_bank; 
    reads_len_buffer_ptr_3.flags  = ddr_bank; 
    weights_arr_buffer_ptr_3.flags  = ddr_bank; 
    new_ref_idx_ptr_3.flags  = ddr_bank; 
    //whd_buffer_ptr_3.flags = ddr_bank;
 
    std::cout << "Preprocess time is : " << duration.count() << " ms\n";

    // Setting input and output objects
    con_arr_buffer_ptr_3.obj = con_arr_buffer_3.data();
    //con_len_buffer_ptr_3.obj = con_len;  // this would cause extra memcpy since it is not aligned
    con_len_buffer_ptr_3.obj = con_len_buffer_3.data();
    reads_arr_buffer_ptr_3.obj = reads_arr_buffer_3.data();
    //reads_len_buffer_ptr_3.obj = reads_len;
    reads_len_buffer_ptr_3.obj = reads_len_buffer_3.data();
    weights_arr_buffer_ptr_3.obj = weights_arr_buffer_3.data();
    //whd_buffer_ptr_3.obj = whd_buffer_3.data();
    new_ref_idx_ptr_3.obj = new_ref_idx_3->data();

    std::cout << "Preprocess time is : " << duration.count() << " ms\n";
    // Setting param to zero 
    con_arr_buffer_ptr_3.param = 0; con_len_buffer_ptr_3.param = 0;
    reads_arr_buffer_ptr_3.param = 0; reads_len_buffer_ptr_3.param = 0; 
    weights_arr_buffer_ptr_3.param = 0; new_ref_idx_ptr_3.param = 0;
    whd_buffer_ptr_3.param = 0;

    printf("Finish creating buffers\n");
    //Allocate Buffer in Global Memory

    cl::Buffer con_arr_input_3(context,CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY | CL_MEM_EXT_PTR_XILINX ,\ 
            CON_SIZE * CON_LEN, &con_arr_buffer_ptr_3);
    cl::Buffer reads_arr_input_3(context,CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY | CL_MEM_EXT_PTR_XILINX ,\ 
            READS_SIZE * READS_LEN, &reads_arr_buffer_ptr_3);
    cl::Buffer weights_arr_input_3(context,CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY | CL_MEM_EXT_PTR_XILINX ,\ 
            READS_SIZE * READS_LEN, &weights_arr_buffer_ptr_3);

    cl::Buffer con_len_input_3(context,CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY | CL_MEM_EXT_PTR_XILINX , 
            CON_SIZE * sizeof(int), &con_len_buffer_ptr_3);
    cl::Buffer reads_len_input_3(context,CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY | CL_MEM_EXT_PTR_XILINX , 
            READS_SIZE * sizeof(int), &reads_len_buffer_ptr_3);

    //cl::Buffer whd_output_3(context, CL_MEM_USE_HOST_PTR | CL_MEM_WRITE_ONLY| CL_MEM_EXT_PTR_XILINX ,
    //        reads_size_3 * con_size * 2 * sizeof(int), &whd_buffer_ptr_3);
    cl::Buffer new_ref_idx_output_3(context, CL_MEM_USE_HOST_PTR | CL_MEM_WRITE_ONLY | CL_MEM_EXT_PTR_XILINX ,
            READS_SIZE * sizeof(int), &new_ref_idx_ptr_3);

    inBufVec.push_back(con_len_input_3);
    inBufVec.push_back(reads_len_input_3);
    inBufVec.push_back(con_arr_input_3);
    inBufVec.push_back(reads_arr_input_3);
    inBufVec.push_back(weights_arr_input_3);
    //outBufVec.push_back(whd_output_3);
    outBufVec.push_back(new_ref_idx_output_3);

    //Set the Kernel Arguments
    krnl_indel.setArg(narg++, con_arr_input_3);
    krnl_indel.setArg(narg++, con_size);
    krnl_indel.setArg(narg++, con_len_input_3);
    krnl_indel.setArg(narg++, reads_arr_input_3);
    krnl_indel.setArg(narg++, reads_size_3);
    krnl_indel.setArg(narg++, reads_len_input_3);
    krnl_indel.setArg(narg++, weights_arr_input_3);
    //krnl_indel.setArg(narg++, whd_output_3);
    krnl_indel.setArg(narg++, new_ref_idx_output_3);
   


//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
    //Copy input data to device global memory
    q.enqueueMigrateMemObjects(inBufVec, 0/* 0 means from host*/);
    //q.finish();
 
    //Launch the Kernel
    //cl::Event event;
    //-------------------
    start = std::chrono::high_resolution_clock::now(); 
    //int work_group = WORK_GROUP;
    int work_group = 1;

    cl::Event events[work_group];

    //for(int i = 0; i < work_group; i++) {
    int i = 0; 
    krnl_indel.setArg(narg+0, i);
    krnl_indel.setArg(narg+1, work_group);

    q.enqueueTask(krnl_indel, NULL, &events[0]);
    //q.finish();

    //events[0].wait();
    //Copy Result from Device Global Memory to Host Local Memory
    q.enqueueMigrateMemObjects(outBufVec, CL_MIGRATE_MEM_OBJECT_HOST);
    //q.finish();
    //OPENCL HOST CODE AREA END   
    //int * whd_buffer_arr = &whd_buffer[0];
    //Indel_Rank(con_size, reads_size, whd_buffer_arr, new_ref_idx); 
    finish = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(finish - start);

    std::cout << "OpenCl Execution time is: " << duration.count() << " ms\n";  

//INST_0============================================//
   if (test_idx_0 < num_tests){ 
      reads_size_arr[test_idx_0] = reads_size_0;
      new_ref_idx_arr[test_idx_0] = new_ref_idx_0;
      new_ref_idx_ref_arr[test_idx_0] = new_ref_idx_ref_0;
    }


//INST_1============================================//
   if (test_idx_1 < num_tests){ 
      reads_size_arr[test_idx_1] = reads_size_1;
      new_ref_idx_arr[test_idx_1] = new_ref_idx_1;
      new_ref_idx_ref_arr[test_idx_1] = new_ref_idx_ref_1;
    }


//INST_2============================================//
   if (test_idx_2 < num_tests){ 
      reads_size_arr[test_idx_2] = reads_size_2;
      new_ref_idx_arr[test_idx_2] = new_ref_idx_2;
      new_ref_idx_ref_arr[test_idx_2] = new_ref_idx_ref_2;
    }


//INST_3============================================//
   if (test_idx_3 < num_tests){ 
      reads_size_arr[test_idx_3] = reads_size_3;
      new_ref_idx_arr[test_idx_3] = new_ref_idx_3;
      new_ref_idx_ref_arr[test_idx_3] = new_ref_idx_ref_3;
    }


    
   // Compare_Results(0, reads_size_0, new_ref_idx_0, new_ref_idx_ref_0);
   // Compare_Results(1, reads_size_1, new_ref_idx_1, new_ref_idx_ref_1);
   // Compare_Results(2, reads_size_2, new_ref_idx_2, new_ref_idx_ref_2);
   // Compare_Results(3, reads_size_3, new_ref_idx_3, new_ref_idx_ref_3);
  }

  for(int i; i < NUM_KERNELS; i ++) {
     qs[i].finish();
  }
  for (int test_idx = 0; test_idx< num_tests; test_idx++) {
      Compare_Results(test_idx, reads_size_arr[test_idx], *(new_ref_idx_arr[test_idx]) , new_ref_idx_ref_arr[test_idx]);
  }

  return EXIT_SUCCESS;
}
