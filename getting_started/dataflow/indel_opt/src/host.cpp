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
// JENNY TODO
// Put input and output onto different memory banks 
// https://github.com/Xilinx/SDAccel_Examples/blob/master/getting_started/kernel_to_gmem/
//
//typedef std::map<char, ap_uint<4>> BasePairMap;
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

int* parse_schedule (const char* file, int* num_tests) {

    FILE *fp;
    fp=fopen(file, "r");
    if (fp == NULL){
        printf("File not found: %s!\n", file);
        exit(EXIT_FAILURE);
    }
    char separators[] = " ";
    char *line = NULL;
    size_t line_len = 0;
    ssize_t read;

    int line_num = 0;
    int base = 0;
    int* file_index = NULL;
    int size = 0;
    while ((read = getline(&line, &line_len, fp)) != -1) {
        //printf("Retrieved line of length %zu :\n", read);
        //printf("%s", line);
        printf("String length%d\n", line_len);
        char* tmp = (char*) malloc(line_len * sizeof(char));
        strcpy(tmp, line);
 
        char *p = strtok(line, separators); 
        assert (p != NULL);
        int i = 0;
        while (p != NULL){
            p = strtok(NULL, separators);
            i++;
        }
        size = i;
        printf("Total Number of Target Files: %d\n", size);
        file_index = (int*) malloc(size * sizeof(int));

        p = strtok(tmp, separators); 
        assert (p != NULL);
        i = 0;
        while (p != NULL){
            file_index[i] = atoi(p);
            p = strtok(NULL, separators);
            i++;
        }
        
        printf("Last index: %d", file_index[size - 1]);
        printf("\n");
        break;
    }
  * num_tests =  size;
  return file_index;
}

void parse(const char* file_prefix, int col_num, char* con_arr, int** con_len_arr, int* con_size, int*total_len) {

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
    *total_len = base;
    *con_size = num_lines;
    free(line);
    fclose(fp); 
}

int main(int argc, char** argv)
{

    printf("Parse scheduele file\n");
    int num_tests = 0;
    int* test_indices = parse_schedule("../indel_tests/ir_toy-schedule.txt", &num_tests);

    if (argc < 2){
        return 1;
    }
    
    m['A'] = 0;
    m['T'] = 1;
    m['C'] = 2;
    m['G'] = 3;
    m['U'] = 4;
 
    //OPENCL HOST CODE AREA START
    std::vector<cl::Device> devices = xcl::get_xil_devices();
    cl::Device device = devices[0];

    cl::Context context(device);
    //cl::CommandQueue q(context, device, CL_QUEUE_PROFILING_ENABLE | CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE);
    cl::CommandQueue q(context, device, CL_QUEUE_PROFILING_ENABLE);
    std::string device_name = device.getInfo<CL_DEVICE_NAME>(); 

    //Create Program and Kernel
    std::string binaryFile = xcl::find_binary_file(device_name, "Indel_Accel");
    cl::Program::Binaries bins = xcl::import_binary_file(binaryFile);
    devices.resize(1);
    cl::Program program(context, devices, bins);
    //cl::Kernel krnl_indel(program,"Indel_Accel");

    std::vector<cl::Kernel> krnl_indels(1,  cl::Kernel(program,"Indel_Accel"));
    std::chrono::high_resolution_clock::time_point start, finish;
    std::chrono::milliseconds duration;

  for (int test_idx = 0; test_idx< num_tests; test_idx+= PARALLEL_UNITS) {

    std::chrono::milliseconds parse_time[PARALLEL_UNITS];
    //int kernel_idx = test_idx % 2;
    cl::Kernel krnl_indel = krnl_indels[0];
    //char* test_num = argv[1];
    std::vector<cl::Memory> inBufVec, outBufVec;

    int * min_whd, * min_whd_idx, * new_ref, * new_ref_idx_ref;
    int * new_ref_idx_ref_0;
    int * new_ref_idx_ref_1;
    int * new_ref_idx_ref_2;
    int * new_ref_idx_ref_3;
    int narg=0;

    Run_Unit(0, test_idx, context, q, krnl_indel, inBufVec, outBufVec, new_ref_idx_ref_0, parse_time, narg, num_tests, test_indices);
    Run_Unit(1, test_idx, context, q, krnl_indel, inBufVec, outBufVec, new_ref_idx_ref_1, parse_time, narg, num_tests, test_indices);
    Run_Unit(2, test_idx, context, q, krnl_indel, inBufVec, outBufVec, new_ref_idx_ref_2, parse_time, narg, num_tests, test_indices);
    Run_Unit(3, test_idx, context, q, krnl_indel, inBufVec, outBufVec, new_ref_idx_ref_3, parse_time, narg, num_tests, test_indices);
 
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

    q.finish();

    //event.wait();
    //Copy Result from Device Global Memory to Host Local Memory
    q.enqueueMigrateMemObjects(outBufVec, CL_MIGRATE_MEM_OBJECT_HOST);
    q.finish();

    //OPENCL HOST CODE AREA END   
    //int * whd_buffer_arr = &whd_buffer[0];
    //Indel_Rank(con_size, reads_size, whd_buffer_arr, new_ref_idx); 
    finish = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(finish - start);

    std::cout << "OpenCl Execution time is: " << duration.count() << " ms\n";  

    //Compare_Results(0);
    //Compare_Results(1);
    //Compare_Results(2);
    //Compare_Results(3);

  }
  return EXIT_SUCCESS;
}
