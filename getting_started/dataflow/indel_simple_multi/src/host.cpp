//OpenCL utility layer include
#include "xcl2.hpp"

#include "oclHelper.h"
#include <vector>
#include <assert.h>
#include <map>

#include <pthread.h>
#include <chrono>
#include "ap_int.h"

#include "Indel_Accel.h"
#include "Indel_Accel_SW.cpp"

//#define WORK_GROUP 4 
//#define WORK_ITEM_PER_GROUP 1
//#define PARALLEL_UNITS 1

// JENNY TODO
// Put input and output onto different memory banks 
// https://github.com/Xilinx/SDAccel_Examples/blob/master/getting_started/kernel_to_gmem/
//
//typedef std::map<char, ap_uint<8>> BasePairMap;
typedef std::map<char, char> BasePairMap;
BasePairMap m;

#include "Indel_Accel_Host.cpp"
// Wrap any OpenCL API calls that return error code(cl_int) with the below macro
// to quickly check for an error
#define OCL_CHECK(call)                                                        \
  do {                                                                         \
    cl_int err = call;                                                         \
    if (err != CL_SUCCESS) {                                                   \
      printf(__FILE__ ":%d: [ERROR] " #call " returned %s\n", __LINE__,        \
             oclErrorCode(err));                                               \
      exit(EXIT_FAILURE);                                                      \
    }                                                                          \
  } while (0);

// Checks OpenCL error codes
void check(cl_int err_code) {
  if (err_code != CL_SUCCESS) {
    printf("ERROR: %d\n", err_code);
    exit(EXIT_FAILURE);
  }
}

// An event callback function that prints the operations performed by the OpenCL
// runtime.
void event_cb(cl_event event, cl_int cmd_status, void *data) {
  cl_command_type command;
  clGetEventInfo(event, CL_EVENT_COMMAND_TYPE, sizeof(cl_command_type),
                 &command, nullptr);
  cl_int status;
  clGetEventInfo(event, CL_EVENT_COMMAND_EXECUTION_STATUS, sizeof(cl_int),
                 &status, nullptr);
  const char *command_str;
  const char *status_str;
  switch (command) {
  case CL_COMMAND_READ_BUFFER:
    command_str = "buffer read";
    break;
  case CL_COMMAND_WRITE_BUFFER:
    command_str = "buffer write";
    break;
  case CL_COMMAND_NDRANGE_KERNEL:
    command_str = "kernel";
    break;
  case CL_COMMAND_MAP_BUFFER:
    command_str = "kernel";
    break;
  case CL_COMMAND_COPY_BUFFER:
    command_str = "kernel";
    break;
  case CL_COMMAND_MIGRATE_MEM_OBJECTS:
        command_str = "buffer migrate";
      break;
  default:
    command_str = "unknown";
  }
  switch (status) {
  case CL_QUEUED:
    status_str = "Queued";
    break;
  case CL_SUBMITTED:
    status_str = "Submitted";
    break;
  case CL_RUNNING:
    status_str = "Executing";
    break;
  case CL_COMPLETE:
    status_str = "Completed";
    break;
  }
  printf("[%s]: %s %s\n", reinterpret_cast<char *>(data), status_str,
         command_str);
  fflush(stdout);
}

// Sets the callback for a particular event
void set_callback(cl_event event, const char *queue_name) {
  OCL_CHECK(
      clSetEventCallback(event, CL_COMPLETE, event_cb, (void *)queue_name));
}


//#define BATCH_SIZE 32
#define NUM_KERNELS 2

void get_filename (const char* file_prefix, char* con, char* reads, int number){

    char test_num[5];
    snprintf(test_num, sizeof(test_num), "%d", number);
    con[0] = '\0'; 
    strcat(con, file_prefix);
    strcat(con, test_num);
    strcat(con, ".SEQ");

    reads[0] = '\0';
    strcat(reads, file_prefix);
    strcat(reads, test_num);
    strcat(reads, ".READS");

    printf("TARGET %s\n", test_num);
}



void pack_ap_uint4(std::vector<char,aligned_allocator<char>>* inout){

    int size = inout->size();
    for (int i = 0; i < size; i++){
        int dst_idx = i >> 1;
        int offset = (i % 2) << 2;
	char cur = (*inout)[i];
	char prev = (*inout)[dst_idx];
        char c = cur & 0xf;
        c = (prev & ~(0xf << offset)) | (c << offset);
	(*inout)[dst_idx] = c;
	if (i < 8) printf("inout_prev %hhX inout_cur %hhX char %hhX %d %d\n", prev, cur, c, i, dst_idx);
       //printf("con_arr_buffer_0[%d + %d]: %hhX \n", idx, offset, con_arr_buffer_0[idx]);
    }
    inout->resize((size >> 1) + 1);
}


template <typename T> 
void parse_file(const char* file_prefix, int col_num, T* input_arr, std::vector<int,aligned_allocator<int>>* len_base, std::vector<int,aligned_allocator<int>>* size_base, int option) {
//void parse(const char* file_prefix, int col_num, char* con_arr, int** con_len_arr, int* con_size, int*total_len) {

    FILE *fp;
    char con[256]="";
    strcat(con, file_prefix);
    strcat(con, ".tbl");

    fp=fopen(con, "r");
    printf("File: %s\n", con);

    if (fp == NULL) {
        printf("File not found: %s!\n", con);
        exit(EXIT_FAILURE);
    }
    // Length of each segment 
    char separators[] = "|";
    char *line = NULL;
    size_t line_len = 0;
    ssize_t read;

    int prev_len_base = 0;
    int prev_size_base = 0;
    if (option < 3) {
        prev_len_base = len_base -> back();
        prev_size_base = size_base -> back(); 
    }
        
    int line_num = 0;
    int base = 0;
    while ((read = getline(&line, &line_len, fp)) != -1) {
        //printf("Retrieved line of length %zu :\n", read);
        //printf("%s", line);

        char *p = strtok(line, separators); 

        assert (p != NULL);

        int i = 1; 
        int input_len = 0;
        while (p != NULL){
            //printf("%s\n", p);
            //target[i++] = atoi(p); 
            p = strtok(NULL, separators);

            if(i == col_num){ 
                //if(line_num == 0){
                //  len = strlen(p); 
                //  input_arr = (char*) aligned_alloc(sizeof(char), len * num_lines * sizeof(char));
                //}
                input_len = strlen(p);
                //fprintf(stderr, "line_num %d strlen %d\n", line_num, input_len[line_num]);
                int k; 
                for (k = 0; k < input_len; k++){
                    //printf("%c", p[k]);
                    //strncpy(&input_arr[base + k], &p[k], 1);
                    if (option == 0){ // this is for elementwise 
   // for(int i = 0 ; i < con_total_len; i++){
   //     int idx = i >> 1;
   //     int offset = (i % 2) << 2;
   //     char c = m[con_arr[i]] & 0xf;
   //     con_arr_buffer_0[idx] = (con_arr_buffer_0[idx] & ~(0xf << offset)) | (c << offset);
   //     //printf("con_arr_buffer_0[%d + %d]: %hhX \n", idx, offset, con_arr_buffer_0[idx]);
                        ;
                    } else if (option == 1) {
                        input_arr -> push_back(m[p[k]]);
                    } else { // this is for weights 
                        input_arr -> push_back(p[k]);
                    }
                    //printf("%d\n", line_num* len + k);
                    //fprintf(stderr, "%c", input_arr[base + k]);
                }
               	//fprintf(stderr, "\n");
            }
            i++;
        }
        base += input_len;
        if(option < 3) len_base->push_back(prev_len_base + base);
        line_num++;
    }

    printf("what\n");
    if(option < 3) {
        size_base->push_back(prev_size_base + line_num); 
        printf("%d\t", size_base->back());
        printf("%d\t", prev_size_base + line_num);
    }
    printf("\n");
    //int arr_size = len * num_lines;
    //printf("%d\n", arr_size);
    //for(int i= 0; i < arr_size; i++ ){
    //  printf("%d ", input_arr[i]);
    //}
    free(line);
    fclose(fp); 
}



int main(int argc, char** argv ){
    printf("Parse scheduele file\n");
    //Total number of target files 
    int num_tests = 0;
    int* test_indices = parse_schedule("../indel_tests/ch22-schedule.txt", &num_tests);
    //int* test_indices = parse_schedule("../indel_tests/ir_toy-schedule.txt", &num_tests);

    num_tests  = 32;
    printf("num_tests: %d \n", num_tests);
    // Set up the mapping for 4bit representation
    m['A'] = 0;
    m['T'] = 1;
    m['C'] = 2;
    m['G'] = 3;
    m['U'] = 4;
    
    // Set up the reads size arr
    //int ** new_ref_idx_ref_arr = (int**) malloc(num_tests * sizeof(int*));
    //std::vector<std::vector<int,aligned_allocator<int>> *> new_ref_idx_arr(num_tests);

    const char* file_prefix = "../indel_tests/ch22-ir/";
    //const char* file_prefix = "../indel_tests/ir_toy/";
    char con[256] = "";
    char reads[256]="";

    // Set up timer
    std::chrono::high_resolution_clock::time_point start, finish;
    std::chrono::milliseconds duration;

//////OPENCL HOST CODE SETUP START////////
    std::vector<cl::Device> devices = xcl::get_xil_devices();
    cl::Device device = devices[0];

    cl::Context context(device);
    //cl::CommandQueue q(context, device, CL_QUEUE_PROFILING_ENABLE | CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE);
    //std::vector<cl::CommandQueue> qs( NUM_KERNELS, cl::CommandQueue(context, device, CL_QUEUE_PROFILING_ENABLE | CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE));
    //std::vector<cl::CommandQueue> qs( NUM_KERNELS, cl::CommandQueue(context, device, CL_QUEUE_PROFILING_ENABLE));
    cl::CommandQueue qs(context, device, CL_QUEUE_PROFILING_ENABLE| CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE);
    //cl::Event events[3];
    std::vector<std::vector<cl::Event>*> events(3);  
    events[0] = new std::vector<cl::Event>(1);
    //events[1] = new std::vector<cl::Event>(1);
    events[2] = new std::vector<cl::Event>(1);

    std::string device_name = device.getInfo<CL_DEVICE_NAME>(); 

    //Create Program and Kernel
    std::string binaryFile = xcl::find_binary_file(device_name, "Indel_Accel");
    cl::Program::Binaries bins = xcl::import_binary_file(binaryFile);
    devices.resize(1);
    cl::Program program(context, devices, bins);
    //cl::Kernel krnl_indel(program,"Indel_Accel");

    std::vector<cl::Kernel> krnl_indels(NUM_KERNELS,  cl::Kernel(program,"Indel_Accel"));
    int kernel_idx = 0;

    // Set up banks mapping
    unsigned ddr_banks[4] = {XCL_MEM_DDR_BANK0, XCL_MEM_DDR_BANK1, XCL_MEM_DDR_BANK2, XCL_MEM_DDR_BANK3};
    unsigned ddr_bank = XCL_MEM_DDR_BANK0;
 
    cl_mem_ext_ptr_t con_arr_buffer_ptr[NUM_KERNELS], reads_arr_buffer_ptr[NUM_KERNELS], weights_arr_buffer_ptr[NUM_KERNELS], con_size_buffer_ptr[NUM_KERNELS], reads_size_buffer_ptr[NUM_KERNELS], con_base_buffer_ptr[NUM_KERNELS], reads_base_buffer_ptr[NUM_KERNELS], new_ref_idx_ptr[NUM_KERNELS]; //whd_buffer_ptr[NUM_KERNELS],  
//////OPENCL HOST CODE SETUP END////////

    std::vector<std::vector<cl::Memory> > inBufVec_arr, outBufVec_arr;
    int BATCH_SIZE = num_tests; 
    for (int test_idx = 0; test_idx< num_tests; test_idx+=BATCH_SIZE ) {


        // Allocate Buffers
        std::vector<char,aligned_allocator<char>> * con_arr_buffer = new std::vector<char,aligned_allocator<char>>();
        std::vector<char,aligned_allocator<char>> * reads_arr_buffer = new std::vector<char,aligned_allocator<char>>();
        std::vector<char,aligned_allocator<char>> * weights_arr_buffer = new std::vector<char,aligned_allocator<char>>();

        std::vector<int,aligned_allocator<int>> * con_size_buffer = new std::vector<int,aligned_allocator<int>>() ;
        std::vector<int,aligned_allocator<int>> * reads_size_buffer = new std::vector<int,aligned_allocator<int>>() ;
        con_size_buffer->push_back(0);
        reads_size_buffer->push_back(0);
        
        std::vector<int,aligned_allocator<int>> * con_base_buffer = new std::vector<int,aligned_allocator<int>>() ;
        std::vector<int,aligned_allocator<int>> * reads_base_buffer = new std::vector<int,aligned_allocator<int>>() ;
        con_base_buffer->push_back(0);
        reads_base_buffer->push_back(0);

        int batch_size = ((num_tests - test_idx) <  BATCH_SIZE) ? (num_tests - test_idx) : BATCH_SIZE; 
        for (int i = 0; i < batch_size; i++){
            get_filename(file_prefix, con, reads, test_indices[test_idx + i]);
            parse_file<std::vector<char,aligned_allocator<char>>>(con, 1, con_arr_buffer, con_base_buffer, con_size_buffer, 1);
            parse_file<std::vector<char,aligned_allocator<char>>>(reads, 4, reads_arr_buffer, reads_base_buffer, reads_size_buffer, 1);
            parse_file<std::vector<char,aligned_allocator<char>> >(reads, 5, weights_arr_buffer, NULL, NULL, 3);
        }


        std::vector<char,aligned_allocator<char>> * con_arr_buffer_copy = new std::vector<char,aligned_allocator<char>>();
        std::vector<char,aligned_allocator<char>> * reads_arr_buffer_copy = new std::vector<char,aligned_allocator<char>>();
	con_arr_buffer_copy->insert(con_arr_buffer_copy->end(), con_arr_buffer->begin(), con_arr_buffer->end());
	reads_arr_buffer_copy->insert(reads_arr_buffer_copy->end(), reads_arr_buffer->begin(), reads_arr_buffer->end());

	pack_ap_uint4(con_arr_buffer);
	pack_ap_uint4(reads_arr_buffer);
        int reads_total_size = reads_size_buffer->back();

        printf("reads_total_size: %d\n", reads_total_size);
        std::vector<int,aligned_allocator<int>> * new_ref_idx_buffer = new std::vector<int,aligned_allocator<int>> (reads_total_size); 
	for (int i = 0; i < reads_total_size; i++) {
		(*new_ref_idx_buffer)[i] = 0;
	}
        std::vector<int,aligned_allocator<int>> * new_ref_idx_ref_buffer = new std::vector<int,aligned_allocator<int>> (); 

        con_arr_buffer_ptr[kernel_idx].flags  = ddr_bank; 
        con_size_buffer_ptr[kernel_idx].flags  = ddr_bank; 
        con_base_buffer_ptr[kernel_idx].flags  = ddr_bank; 
        reads_arr_buffer_ptr[kernel_idx].flags  = ddr_bank; 
        reads_size_buffer_ptr[kernel_idx].flags  = ddr_bank; 
        reads_base_buffer_ptr[kernel_idx].flags  = ddr_bank; 
        weights_arr_buffer_ptr[kernel_idx].flags  = ddr_bank; 
        new_ref_idx_ptr[kernel_idx].flags  = ddr_bank; 
        //whd_buffer_ptr[kernel_idx].flags = ddr_bank;
     
        // Setting input and output objects
        con_arr_buffer_ptr[kernel_idx].obj = con_arr_buffer->data();
	for(int i =0; i<4; i++){
		printf("con_buffer: %hhX\n", con_arr_buffer->data()[i]);
	}
        con_size_buffer_ptr[kernel_idx].obj = con_size_buffer->data();
        con_base_buffer_ptr[kernel_idx].obj = con_base_buffer->data();
        reads_arr_buffer_ptr[kernel_idx].obj = reads_arr_buffer->data();
        reads_size_buffer_ptr[kernel_idx].obj = reads_size_buffer->data();
        reads_base_buffer_ptr[kernel_idx].obj = reads_base_buffer->data();
        weights_arr_buffer_ptr[kernel_idx].obj = weights_arr_buffer->data();
        //whd_buffer_ptr[kernel_idx].obj = whd_buffer_0.data();
        new_ref_idx_ptr[kernel_idx].obj = new_ref_idx_buffer->data();

        std::cout << "Preprocess time is : " << duration.count() << " ms\n";
        // Setting param to zero 
        con_arr_buffer_ptr[kernel_idx].param = 0; con_size_buffer_ptr[kernel_idx].param = 0; con_base_buffer_ptr[kernel_idx].param = 0;
        reads_arr_buffer_ptr[kernel_idx].param = 0; reads_size_buffer_ptr[kernel_idx].param = 0; reads_base_buffer_ptr[kernel_idx].param = 0; 
        weights_arr_buffer_ptr[kernel_idx].param = 0; new_ref_idx_ptr[kernel_idx].param = 0;
        //whd_buffer_ptr[kernel_idx].param = 0;


        //Allocate Buffer in Global Memory
        cl::Buffer con_arr_input(context,CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY | CL_MEM_EXT_PTR_XILINX ,\
                con_arr_buffer->size(), &con_arr_buffer_ptr[kernel_idx]);
        cl::Buffer reads_arr_input(context,CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY | CL_MEM_EXT_PTR_XILINX ,\
                reads_arr_buffer->size(), &reads_arr_buffer_ptr[kernel_idx]);
        cl::Buffer weights_arr_input(context,CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY | CL_MEM_EXT_PTR_XILINX ,\
                weights_arr_buffer->size(), &weights_arr_buffer_ptr[kernel_idx]);
        cl::Buffer con_size_input(context,CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY | CL_MEM_EXT_PTR_XILINX ,\
                con_size_buffer->size() * sizeof(int), &con_size_buffer_ptr[kernel_idx]);
        cl::Buffer reads_size_input(context,CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY | CL_MEM_EXT_PTR_XILINX ,\ 
                reads_size_buffer->size() * sizeof(int), &reads_size_buffer_ptr[kernel_idx]);
        cl::Buffer con_base_input(context,CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY | CL_MEM_EXT_PTR_XILINX ,\
                con_base_buffer->size() * sizeof(int), &con_base_buffer_ptr[kernel_idx]);
        cl::Buffer reads_base_input(context,CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY | CL_MEM_EXT_PTR_XILINX ,\
                reads_base_buffer->size() * sizeof(int), &reads_base_buffer_ptr[kernel_idx]);

        //cl::Buffer whd_output(context, CL_MEM_USE_HOST_PTR | CL_MEM_WRITE_ONLY| CL_MEM_EXT_PTR_XILINX ,
        //        reads_size * con_size * 2 * sizeof(int), &whd_buffer_ptr[kernel_idx]);
        cl::Buffer new_ref_idx_output(context, CL_MEM_USE_HOST_PTR | CL_MEM_WRITE_ONLY | CL_MEM_EXT_PTR_XILINX ,
                reads_size_buffer->back() * sizeof(int), &new_ref_idx_ptr[kernel_idx]);


        printf("Finish creating buffers\n");
        std::vector<cl::Memory> inBufVec, outBufVec;

        inBufVec_arr.push_back(inBufVec);
        outBufVec_arr.push_back(outBufVec);
        inBufVec_arr.back().push_back(con_size_input);
        inBufVec_arr.back().push_back(reads_size_input);
        inBufVec_arr.back().push_back(con_base_input);
        inBufVec_arr.back().push_back(reads_base_input);
        inBufVec_arr.back().push_back(con_arr_input);
        inBufVec_arr.back().push_back(reads_arr_input);
        inBufVec_arr.back().push_back(weights_arr_input);
        //outBufVec.push_back(whd_output);
        outBufVec_arr.back().push_back(new_ref_idx_output);

        //Set the Kernel Arguments
        int narg=0;
        krnl_indels[kernel_idx].setArg(narg++, con_arr_input);
        krnl_indels[kernel_idx].setArg(narg++, con_size_input);
        krnl_indels[kernel_idx].setArg(narg++, con_base_input);
        krnl_indels[kernel_idx].setArg(narg++, reads_arr_input);
        krnl_indels[kernel_idx].setArg(narg++, reads_size_input);
        krnl_indels[kernel_idx].setArg(narg++, reads_base_input);
        krnl_indels[kernel_idx].setArg(narg++, weights_arr_input);
        //krnl_indels[kernel_idx].setArg(narg++, whd_output);
        krnl_indels[kernel_idx].setArg(narg++, new_ref_idx_output);


       
        //OCL_CHECK(qs[kernel_idx].enqueueMigrateMemObjects(inBufVec_arr.back(), 0/* 0 means from host*/, NULL, &((*events[0])[0])));
        OCL_CHECK(qs.enqueueMigrateMemObjects(inBufVec_arr.back(), 0/* 0 means from host*/, NULL, &((*events[0])[0])));
 	//set_callback((*events[0])[0], "mem_cpy");
        qs.finish();

        //qs[kernel_idx].finish();
     
        //Launch the Kernel
        //cl::Event event;
        start = std::chrono::high_resolution_clock::now(); 

        int work_group = batch_size >> 2;
        //events[1] = new std::vector<cl::Event>(work_group);
        for(int i = 0; i < work_group; i++) {
            krnl_indels[kernel_idx].setArg(narg+0, i);
            krnl_indels[kernel_idx].setArg(narg+1, 4);

            //q.enqueueTask(krnl_indels[kernel_idx], NULL, &events[0]);
            //qs[kernel_idx].enqueueTask(krnl_indels[kernel_idx], NULL, &events[kernel_idx]);
            //OCL_CHECK(qs[kernel_idx].enqueueTask(krnl_indels[kernel_idx], events[0], &((*events[1])[i])));
            //OCL_CHECK(qs[kernel_idx].enqueueTask(krnl_indels[kernel_idx], events[0], NULL ));
            //OCL_CHECK(qs.enqueueTask(krnl_indels[kernel_idx], events[0], NULL ));
            size_t global = 1;
            size_t local = 1;
            OCL_CHECK(qs.enqueueNDRangeKernel(krnl_indels[kernel_idx], 0, global , local, events[0], NULL ));

        	//qs.finish();
        //qs[kernel_idx].finish();
            //q.finish();
        }


        //event.wait();
        //Copy Result from Device Global Memory to Host Local Memory
        //q.enqueueMigrateMemObjects(outBufVec, CL_MIGRATE_MEM_OBJECT_HOST);
        
        
      	//qs[kernel_idx].finish();
        qs.finish();
        OCL_CHECK(qs.enqueueMigrateMemObjects(outBufVec_arr.back(), CL_MIGRATE_MEM_OBJECT_HOST, events[1], &((*events[2])[0])));

 	//set_callback((*events[2])[0], "mem_cpy_from_fpga");
        //OCL_CHECK(qs[kernel_idx].enqueueMigrateMemObjects(outBufVec_arr.back(), CL_MIGRATE_MEM_OBJECT_HOST, events[1], &((*events[2])[0])));
        //qs[kernel_idx].finish();
        qs.finish();
        //OPENCL HOST CODE AREA END   
        //int * whd_buffer_arr = &whd_buffer[0];
        //Indel_Rank(con_size, reads_size, whd_buffer_arr, new_ref_idx); 
        finish = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::milliseconds>(finish - start);
        std::cout << "OpenCl Execution time is: " << duration.count() << " ms\n";  

        //generate_ref(std::vector<ap_ustd::vector<int><4>>* consensus, std::vector<int>* con_size, std::vector<int>* con_base, \
        //std::vector<ap_ustd::vector<int><4>>* reads, std::vector<int>* reads_size, std::vector<int>* reads_base, char* qs, std::vector<int>* new_ref_idx) {
        //generate_ref(con_arr_buffer_copy, con_size_buffer, con_base_buffer, reads_arr_buffer_copy, reads_size_buffer, reads_base_buffer, weights_arr_buffer, new_ref_idx_ref_buffer);
        //compare_results(new_ref_idx_buffer, new_ref_idx_ref_buffer, reads_size_buffer);
        //print_results(new_ref_idx_ref_buffer, reads_size_buffer);
        print_results(new_ref_idx_buffer, reads_size_buffer);

	delete (new_ref_idx_buffer); 
    }

}



