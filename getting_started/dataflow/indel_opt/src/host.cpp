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

#include <chrono>
#include "ap_int.h"

#include "Indel_Accel.h"
#include "Indel_Accel_SW.cpp"

#define WORK_GROUP 4 
#define WORK_ITEM_PER_GROUP 1
// JENNY TODO
// Put input and output onto different memory banks 
// https://github.com/Xilinx/SDAccel_Examples/blob/master/getting_started/kernel_to_gmem/
//
//typedef std::map<char, ap_uint<4>> BasePairMap;
typedef std::map<char, char> BasePairMap;

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



int count_lines(char* filename){

    FILE *fp;
    fp = fopen(filename, "r");

    // Check if file exists
    if (fp == NULL)
    {
        printf("Could not open file %s\n", filename);
        return 0;
    }

    int count = 0;
    // Extract characters from file and store in character c
    char *line = NULL;
    size_t line_len = 0;
    ssize_t read;

    while ((read = getline(&line, &line_len, fp)) != -1) {
        //for (c = getc(fp); c != EOF; c = getc(fp))
        //  if (c == '\n') // Increment count if this character is newline
        count = count + 1;
    }
    // Close the file
    fclose(fp);
    return count;
}

void parse(const char* file_prefix, int col_num, char* con_arr, int** con_len_arr, int* con_size) {

    FILE *fp;
    char con[256]="";
    strcat(con, file_prefix);
    strcat(con, ".tbl");

    int num_lines = count_lines(con);
    //printf("num_lines %d\n", num_lines);
    fp=fopen(con, "r");

    if (fp == NULL)
        exit(EXIT_FAILURE);

    // Length of each segment 
    int * con_len = (int*) malloc(sizeof(int) * num_lines);
    assert(con_len != NULL);
    * con_len_arr = con_len;
    // Location of the end of each segment  
    int* pos = (int*) malloc(sizeof(int) * num_lines);
    assert(pos != NULL);

    char separators[] = "|";
    char *line = NULL;
    size_t line_len = 0;
    ssize_t read;

    int line_num = 0;
    int base = 0;
    while ((read = getline(&line, &line_len, fp)) != -1) {
        //printf("Retrieved line of length %zu :\n", read);
        //printf("%s", line);

        char *p = strtok(line, separators); 

        assert (p != NULL);

        int i = 1; 
        while (p != NULL){
            //printf("%s\n", p);
            //target[i++] = atoi(p); 
            p = strtok(NULL, separators);

            if(i == col_num){ 
                //if(line_num == 0){
                //  len = strlen(p); 
                //  con_arr = (char*) malloc(len * num_lines * sizeof(char));
                //}
                con_len[line_num] = strlen(p);
                //fprintf(stderr, "line_num %d strlen %d\n", line_num, con_len[line_num]);
                int k; 
                for (k = 0; k < con_len[line_num]; k++){
                    //printf("%c", p[k]);
                    strncpy(&con_arr[base + k], &p[k], 1);
                    //printf("%d\n", line_num* len + k);
                    //fprintf(stderr, "%c", con_arr[base + k]);
                }
               	//fprintf(stderr, "\n");
            }
            i++;
        }
        base += con_len[line_num];
        pos[line_num] = base; 
        line_num++;
    }
    //printf("%d\n", len);
    //int arr_size = len * num_lines;
    //printf("%d\n", arr_size);
    //for(int i= 0; i < arr_size; i++ ){
    //  printf("%d ", con_arr[i]);
    //}

    *con_size = num_lines;
    free(line);
    fclose(fp); 
}

int main(int argc, char** argv)
{

    if (argc < 2){
        return 1;
    }
    
    char* test_num = argv[1];
    const char* file_prefix = "../indel_tests/ir_toy/";
    char con[256] = "";
    strcat(con, file_prefix);
    strcat(con, test_num);
    strcat(con, ".SEQ");

    char reads[256]="";
    strcat(reads, file_prefix);
    strcat(reads, test_num);
    strcat(reads, ".READS");

    printf("TARGET %s\n", test_num);

    BasePairMap m;
    m['A'] = 0;
    m['T'] = 1;
    m['C'] = 2;
    m['G'] = 3;
    m['U'] = 4;

    char* con_arr, *reads_arr, *weights_arr; 
    int *con_len, con_size, *reads_len, reads_size; 

    // Malloc the largest array 
  auto start = std::chrono::high_resolution_clock::now();
    con_arr = (char*) malloc( CON_SIZE * CON_LEN * sizeof(char));
    reads_arr = (char*) malloc( READS_SIZE * READS_LEN * sizeof(char));
    weights_arr = (char*) malloc( READS_SIZE * READS_LEN * sizeof(char));

    parse( con,1, con_arr, &con_len, &con_size);
    parse( reads, 4, reads_arr, &reads_len, &reads_size);
    parse( reads, 5, weights_arr, &reads_len, &reads_size);

  auto finish = std::chrono::high_resolution_clock::now();
  //std::chrono::duration<double> elapsed = finish - start;   
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(finish - start);
  std::cout << "Parsing time is :" << duration.count() << " ms\n";
  //printf("Parse time is: %0.3f ms\n", duration.count());



    int* min_whd = (int*) malloc(con_size * reads_size * sizeof(int));
    int* min_whd_idx = (int*) malloc(con_size * reads_size * sizeof(int));
    int* new_ref = (int*) malloc(reads_size * sizeof(int));
    //int* new_ref_idx = (int*) malloc(reads_size * sizeof(int));
    int* new_ref_idx_ref = (int*) malloc(reads_size * sizeof(int));

    whd(con_arr, con_size, con_len, reads_arr, reads_size, reads_len, weights_arr, min_whd, min_whd_idx);
    score_whd (min_whd,  min_whd_idx, con_size, reads_size, new_ref, new_ref_idx_ref);

    printf("Software Success!\n");
    //return 0;
    //Allocate Memory in Host Memory
    //std::vector<char,aligned_allocator<char>> con_arr_buffer     ( con_arr, con_arr + CON_SIZE * CON_LEN);
    //std::vector<char,aligned_allocator<char>> reads_arr_buffer     (reads_arr, reads_arr + READS_SIZE * READS_LEN);
    //std::vector<char,aligned_allocator<char>> weights_arr_buffer     (weights_arr, weights_arr + READS_SIZE * READS_LEN);

    start = std::chrono::high_resolution_clock::now();
    //std::vector<ap_uint<4>,aligned_allocator<ap_uint<4>>> con_arr_buffer     ( CON_SIZE * CON_LEN);
    std::vector<char,aligned_allocator<char>> con_arr_buffer     ( CON_SIZE * CON_LEN >> 1);
    //std::vector<ap_uint<4>,aligned_allocator<ap_uint<4>>> reads_arr_buffer     ( READS_SIZE * READS_LEN);
    std::vector<char,aligned_allocator<char>> reads_arr_buffer     ( READS_SIZE * READS_LEN >> 1);
    std::vector<char,aligned_allocator<char>> weights_arr_buffer     (weights_arr, weights_arr + READS_SIZE * READS_LEN);

  //  for(int i = 0 ; i < CON_SIZE * CON_LEN; i++){
  //      con_arr_buffer[i] = m[con_arr[i]];
  //  }

    for(int i = 0 ; i < CON_SIZE * CON_LEN; i++){
        int idx = i >> 1;
        int offset = (i % 2) << 2;
        char c = m[con_arr[i]] & 0xf;
        con_arr_buffer[idx] = (con_arr_buffer[idx] & ~(0xf << offset)) | (c << offset);
        //printf("con_arr_buffer[%d + %d]: %hhX \n", idx, offset, con_arr_buffer[idx]);
    }

 //   printf("Read Buffer:");     
 //   for(int i = 0 ; i < READS_LEN * READS_SIZE; i++){
 //       reads_arr_buffer[i] = m[reads_arr[i]];
 //       //unsigned char print_var = reads_arr_buffer[i];
 //       //printf("%x", print_var);     
 //   }
    for(int i = 0 ; i < READS_LEN * READS_SIZE; i++){
        int idx = i >> 1;
        int offset = (i % 2) << 2;
        char c = m[reads_arr[i]] & 0xf;
        reads_arr_buffer[idx] = (reads_arr_buffer[idx] & ~(0xf << offset)) | (c << offset);
        //printf("reads_arr_buffer[%d]: %hhX \n", idx, (char)reads_arr_buffer[idx]);
    }

    std::vector<int,aligned_allocator<int>> con_len_buffer     (con_len, con_len + CON_SIZE);
    std::vector<int,aligned_allocator<int>> reads_len_buffer     (reads_len, reads_len + READS_SIZE);

    std::vector<int,aligned_allocator<int>> whd_buffer(con_size * reads_size << 1);
    //std::vector<int,aligned_allocator<int>> new_ref_idx(reads_size);
    int * new_ref_idx = new int [reads_size];

    for(int i = 0 ; i < reads_size; i++){
        new_ref_idx[i] = 0;
    }
    printf("\n");
    finish = std::chrono::high_resolution_clock::now();

    duration = std::chrono::duration_cast<std::chrono::milliseconds>(finish - start);
    std::cout << "Preprocess time is : " << duration.count() << " ms\n";
    //printf("Preprocess time is: %0.3f ms\n", duration.count());



    //OPENCL HOST CODE AREA START
    std::vector<cl::Device> devices = xcl::get_xil_devices();
    cl::Device device = devices[0];

    cl::Context context(device);
    cl::CommandQueue q(context, device, CL_QUEUE_PROFILING_ENABLE);
    std::string device_name = device.getInfo<CL_DEVICE_NAME>(); 

    //Create Program and Kernel
    std::string binaryFile = xcl::find_binary_file(device_name, "Indel_Accel");
    cl::Program::Binaries bins = xcl::import_binary_file(binaryFile);
    devices.resize(1);
    cl::Program program(context, devices, bins);
    cl::Kernel krnl_indel(program,"Indel_Accel");

    cl_mem_ext_ptr_t con_arr_buffer_ptr, reads_arr_buffer_ptr, weights_arr_buffer_ptr, con_len_buffer_ptr, reads_len_buffer_ptr, whd_buffer_ptr,  new_ref_idx_ptr; 
    con_arr_buffer_ptr.flags  = XCL_MEM_DDR_BANK0; 
    con_len_buffer_ptr.flags  = XCL_MEM_DDR_BANK0; 
    reads_arr_buffer_ptr.flags  = XCL_MEM_DDR_BANK1; 
    reads_len_buffer_ptr.flags  = XCL_MEM_DDR_BANK1; 
    weights_arr_buffer_ptr.flags  = XCL_MEM_DDR_BANK2; 
    //new_ref_idx_ptr.flags  = XCL_MEM_DDR_BANK2; 
    whd_buffer_ptr.flags = XCL_MEM_DDR_BANK3;
 

    // Setting input and output objects
    con_arr_buffer_ptr.obj = con_arr_buffer.data();
    con_len_buffer_ptr.obj = con_len_buffer.data();
    reads_arr_buffer_ptr.obj = reads_arr_buffer.data();
    reads_len_buffer_ptr.obj = reads_len_buffer.data();
    weights_arr_buffer_ptr.obj = weights_arr_buffer.data();
    whd_buffer_ptr.obj = whd_buffer.data();
    //new_ref_idx_ptr.obj = new_ref_idx.data();

    // Setting param to zero 
    con_arr_buffer_ptr.param = 0; con_len_buffer_ptr.param = 0;
    reads_arr_buffer_ptr.param = 0; reads_len_buffer_ptr.param = 0; 
    weights_arr_buffer_ptr.param = 0; new_ref_idx_ptr.param = 0;
    whd_buffer_ptr.param = 0;

    //Allocate Buffer in Global Memory

    std::vector<cl::Memory> inBufVec, outBufVec;
    cl::Buffer con_arr_input(context,CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY | CL_MEM_EXT_PTR_XILINX ,\ 
            CON_SIZE * CON_LEN, &con_arr_buffer_ptr);
    cl::Buffer reads_arr_input(context,CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY | CL_MEM_EXT_PTR_XILINX ,\ 
            READS_SIZE * READS_LEN, &reads_arr_buffer_ptr);
    cl::Buffer weights_arr_input(context,CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY | CL_MEM_EXT_PTR_XILINX ,\ 
            READS_SIZE * READS_LEN, &weights_arr_buffer_ptr);

    std::vector<cl::Memory> con_len_vec, reads_len_vec;
    cl::Buffer con_len_input(context,CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY | CL_MEM_EXT_PTR_XILINX , 
            CON_SIZE * sizeof(int), &con_len_buffer_ptr);
    cl::Buffer reads_len_input(context,CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY | CL_MEM_EXT_PTR_XILINX , 
            READS_SIZE * sizeof(int), &reads_len_buffer_ptr);

    cl::Buffer whd_output(context, CL_MEM_USE_HOST_PTR | CL_MEM_WRITE_ONLY| CL_MEM_EXT_PTR_XILINX ,
            reads_size * con_size * 2 * sizeof(int), &whd_buffer_ptr);
    //cl::Buffer new_ref_idx_output(context, CL_MEM_USE_HOST_PTR | CL_MEM_WRITE_ONLY | CL_MEM_EXT_PTR_XILINX ,
    //        READS_SIZE * sizeof(int), &new_ref_idx_ptr);

    inBufVec.push_back(con_len_input);
    inBufVec.push_back(reads_len_input);
    inBufVec.push_back(con_arr_input);
    inBufVec.push_back(reads_arr_input);
    inBufVec.push_back(weights_arr_input);
    outBufVec.push_back(whd_output);
    //outBufVec.push_back(new_ref_idx_output);

    //Copy input data to device global memory
    q.enqueueMigrateMemObjects(inBufVec, 0/* 0 means from host*/);

    //Set the Kernel Arguments
    int narg=0;
    krnl_indel.setArg(narg++, con_arr_input);
    krnl_indel.setArg(narg++, con_size);
    krnl_indel.setArg(narg++, con_len_input);
    krnl_indel.setArg(narg++, reads_arr_input);
    krnl_indel.setArg(narg++, reads_size);
    krnl_indel.setArg(narg++, reads_len_input);
    krnl_indel.setArg(narg++, weights_arr_input);
    krnl_indel.setArg(narg++, whd_output);

    //krnl_indel.setArg(narg++, new_ref_idx_output);

    //Launch the Kernel
    //cl::Event event;

    start = std::chrono::high_resolution_clock::now(); 
    int work_group = WORK_GROUP;

    cl::Event events[work_group];

    for(int i = 0; i < work_group; i++) {
        krnl_indel.setArg(narg+0, i);
        krnl_indel.setArg(narg+1, work_group);

        q.enqueueTask(krnl_indel, NULL, &events[i]);
    }

    q.finish();

    //event.wait();
    //Copy Result from Device Global Memory to Host Local Memory
    q.enqueueMigrateMemObjects(outBufVec, CL_MIGRATE_MEM_OBJECT_HOST);
    q.finish();

    //OPENCL HOST CODE AREA END
    
    int * whd_buffer_arr = &whd_buffer[0];
    Indel_Rank(con_size, reads_size, whd_buffer_arr, new_ref_idx); 
    finish = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(finish - start);

    std::cout << "OpenCl Execution time is: " << duration.count() << " ms\n";  

    // Compare the results of the Device to the simulation
    int match = 0;
    for (int i = 0 ; i < reads_size; i++){
        if (new_ref_idx[i] != new_ref_idx_ref[i]){
            std::cout << "Error: Result mismatch" << std::endl;
            std::cout << "i = " << i << " CPU result = " << new_ref_idx_ref[i]
                << " Device result = " << new_ref_idx[i] << std::endl;
            match = 1;
            break;
        }
    }

    delete[] new_ref_idx;
    std::cout << "TEST " << (match ? "FAILED" : "PASSED") << std::endl; 
    return (match ? EXIT_FAILURE :  EXIT_SUCCESS);
}
