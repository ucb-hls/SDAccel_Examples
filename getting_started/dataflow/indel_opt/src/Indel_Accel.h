#ifndef INDEL_ACCEL_H
#define INDEL_ACCEL_H

//#define CON_SIZE 2
//#define READS_SIZE 8 
#define CON_SIZE 32
#define READS_SIZE 256
#define CON_LEN 2048
#define READS_LEN 256


#define CREATE_INPUT(X) {
    //cl::CommandQueue q_ ## X (context, device, CL_QUEUE_PROFILING_ENABLE);
    int test_idx_## X= test_idx + X;

    if (test_idx ## X < num_tests){
        snprintf(test_num, sizeof(test_num), "%d", test_indices[test_idx ## X]);
    
        con[0] = "\0";
        strcat(con, file_prefix);
        strcat(con, test_num);
        strcat(con, ".SEQ");

        reads[0]= "\0";
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
        parse( reads, 4, reads_arr, &reads_len, &reads_size, &reads_total_len);
        parse( reads, 5, weights_arr, &reads_len, &reads_size, &weights_total_len);

        finish = std::chrono::high_resolution_clock::now();
        //std::chrono::duration<double> elapsed = finish - start;   
        duration = std::chrono::duration_cast<std::chrono::milliseconds>(finish - start);
        parse_time[test_local_idx] = duration;
        //std::cout << "Parsing time is :" << duration.count() << " ms\n";
        //printf("Parse time is: %0.3f ms\n", duration.count());

        min_whd = (int*) malloc(con_size * reads_size * sizeof(int));
        min_whd_idx = (int*) malloc(con_size * reads_size * sizeof(int));
        new_ref = (int*) malloc(reads_size * sizeof(int));
        //int* new_ref_idx = (int*) malloc(reads_size * sizeof(int));
        int* new_ref_idx_ref = (int*) malloc(reads_size * sizeof(int));

        whd(con_arr, con_size, con_len, reads_arr, reads_size, reads_len, weights_arr, min_whd, min_whd_idx);
        score_whd (min_whd,  min_whd_idx, con_size, reads_size, new_ref, new_ref_idx_ref);

        printf("Software Success!\n");
    } else {
        con_size=reads_size=0;
        con_total_len=reads_total_len=0;
        con_arr=reads_arr=weights_arr=NULL;
        con_len=reads_len=NULL;
        min_whd=min_whd_idx=new_ref=new_ref_idx_ref;   
    }

    //return 0;
    //Allocate Memory in Host Memory
    //std::vector<char,aligned_allocator<char>> con_arr_buffer     ( con_arr, con_arr + CON_SIZE * CON_LEN);
    //std::vector<char,aligned_allocator<char>> reads_arr_buffer     (reads_arr, reads_arr + READS_SIZE * READS_LEN);
    //std::vector<char,aligned_allocator<char>> weights_arr_buffer     (weights_arr, weights_arr + READS_SIZE * READS_LEN);

    start = std::chrono::high_resolution_clock::now();
    //std::vector<ap_uint<4>,aligned_allocator<ap_uint<4>>> con_arr_buffer     ( CON_SIZE * CON_LEN);
    std::vector<char,aligned_allocator<char>> con_arr_buffer_ ## X      ( con_size * con_total_len >> 1);
    //std::vector<ap_uint<4>,aligned_allocator<ap_uint<4>>> reads_arr_buffer     ( READS_SIZE * READS_LEN);
    std::vector<char,aligned_allocator<char>> reads_arr_buffer_ ## X    ( reads_size * reads_total_len >> 1);
    std::vector<char,aligned_allocator<char>> weights_arr_buffer_ ## X  (weights_arr, weights_arr + reads_size* reads_total_len);

  //  for(int i = 0 ; i < CON_SIZE * CON_LEN; i++){
  //      con_arr_buffer[i] = m[con_arr[i]];
  //  }

    //for(int i = 0 ; i < CON_SIZE * CON_LEN; i++){
    for(int i = 0 ; i < con_size*con_total_len; i++){
        int idx = i >> 1;
        int offset = (i % 2) << 2;
        char c = m[con_arr[i]] & 0xf;
        con_arr_buffer_ ## X[idx] = (con_arr_buffer_ ## X[idx] & ~(0xf << offset)) | (c << offset);
        //printf("con_arr_buffer_ ## X[%d + %d]: %hhX \n", idx, offset, con_arr_buffer_ ## X[idx]);
    }

 //   printf("Read Buffer:");     
 //   for(int i = 0 ; i < READS_LEN * READS_SIZE; i++){
 //       reads_arr_buffer[i] = m[reads_arr[i]];
 //       //unsigned char print_var = reads_arr_buffer[i];
 //       //printf("%x", print_var);     
 //   }
    for(int i = 0 ; i < reads_size*reads_total_len; i++){
        int idx = i >> 1;
        int offset = (i % 2) << 2;
        char c = m[reads_arr[i]] & 0xf;
        reads_arr_buffer_ ## X[idx] = (reads_arr_buffer_ ## X[idx] & ~(0xf << offset)) | (c << offset);
        //printf("reads_arr_buffer_ ## X[%d]: %hhX \n", idx, (char)reads_arr_buffer_ ## X[idx]);
    }

    std::vector<int,aligned_allocator<int>> con_len_buffer_ ## X     (con_len, con_len + con_size);
    std::vector<int,aligned_allocator<int>> reads_len_buffer_ ## X     (reads_len, reads_len + reads_size);

    //std::vector<int,aligned_allocator<int>> whd_buffer(con_size * reads_size << 1);
    std::vector<int,aligned_allocator<int>> new_ref_idx_ ## X(reads_size);
    //int * new_ref_idx = new int [reads_size];

    for(int i = 0 ; i < reads_size; i++){
        new_ref_idx_ ## X[i] = 0;
    }
    printf("\n");
    finish = std::chrono::high_resolution_clock::now();

    duration = std::chrono::duration_cast<std::chrono::milliseconds>(finish - start);
    std::cout << "Preprocess time is : " << duration.count() << " ms\n";
    //printf("Preprocess time is: %0.3f ms\n", duration.count());


    cl_mem_ext_ptr_ ## X_t con_arr_buffer_ptr_ ## X, reads_arr_buffer_ptr_ ## X, weights_arr_buffer_ptr_ ## X, con_len_buffer_ptr_ ## X, reads_len_buffer_ptr_ ## X, whd_buffer_ptr_ ## X,  new_ref_idx_ptr_ ## X; 
    con_arr_buffer_ptr_ ## X.flags  = XCL_MEM_DDR_BANK ## X; 
    con_len_buffer_ptr_ ## X.flags  = XCL_MEM_DDR_BANK ## X; 
    reads_arr_buffer_ptr_ ## X.flags  = XCL_MEM_DDR_BANK ## X; 
    reads_len_buffer_ptr_ ## X.flags  = XCL_MEM_DDR_BANK ## X; 
    weights_arr_buffer_ptr_ ## X.flags  = XCL_MEM_DDR_BANK ## X; 
    new_ref_idx_ptr_ ## X.flags  = XCL_MEM_DDR_BANK ## X; 
    //whd_buffer_ptr_ ## X.flags = XCL_MEM_DDR_BANK3;
 

    // Setting input and output objects
    con_arr_buffer_ptr_ ## X.obj = con_arr_buffer_ ## X.data();
    con_len_buffer_ptr_ ## X.obj = con_len_buffer_ ## X.data();
    reads_arr_buffer_ptr_ ## X.obj = reads_arr_buffer_ ## X.data();
    reads_len_buffer_ptr_ ## X.obj = reads_len_buffer_ ## X.data();
    weights_arr_buffer_ptr_ ## X.obj = weights_arr_buffer_ ## X.data();
    //whd_buffer_ptr_ ## X.obj = whd_buffer_ ## X.data();
    new_ref_idx_ptr_ ## X.obj = new_ref_idx_ ## X.data();

    // Setting param to zero 
    con_arr_buffer_ptr_ ## X.param = 0; con_len_buffer_ptr_ ## X.param = 0;
    reads_arr_buffer_ptr_ ## X.param = 0; reads_len_buffer_ptr_ ## X.param = 0; 
    weights_arr_buffer_ptr_ ## X.param = 0; new_ref_idx_ptr_ ## X.param = 0;
    whd_buffer_ptr_ ## X.param = 0;

    //Allocate Buffer in Global Memory

    cl::Buffer con_arr_input_ ## X(context,CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY | CL_MEM_EXT_PTR_XILINX ,\ 
            CON_SIZE * CON_LEN, &con_arr_buffer_ptr_ ## X);
    cl::Buffer reads_arr_input_ ## X(context,CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY | CL_MEM_EXT_PTR_XILINX ,\ 
            READS_SIZE * READS_LEN, &reads_arr_buffer_ptr_ ## X);
    cl::Buffer weights_arr_input_ ## X(context,CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY | CL_MEM_EXT_PTR_XILINX ,\ 
            READS_SIZE * READS_LEN, &weights_arr_buffer_ptr_ ## X);

    std::vector<cl::Memory> con_len_vec, reads_len_vec;
    cl::Buffer con_len_input_ ## X(context,CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY | CL_MEM_EXT_PTR_XILINX , 
            CON_SIZE * sizeof(int), &con_len_buffer_ptr_ ## X);
    cl::Buffer reads_len_input_ ## X(context,CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY | CL_MEM_EXT_PTR_XILINX , 
            READS_SIZE * sizeof(int), &reads_len_buffer_ptr_ ## X);

    //cl::Buffer whd_output_ ## X(context, CL_MEM_USE_HOST_PTR | CL_MEM_WRITE_ONLY| CL_MEM_EXT_PTR_XILINX ,
    //        reads_size * con_size * 2 * sizeof(int), &whd_buffer_ptr_ ## X);
    cl::Buffer new_ref_idx_output_ ## X(context, CL_MEM_USE_HOST_PTR | CL_MEM_WRITE_ONLY | CL_MEM_EXT_PTR_XILINX ,
            READS_SIZE * sizeof(int), &new_ref_idx_ptr_ ## X);

    inBufVec.push_back(con_len_input_ ## X);
    inBufVec.push_back(reads_len_input_ ## X);
    inBufVec.push_back(con_arr_input_ ## X);
    inBufVec.push_back(reads_arr_input_ ## X);
    inBufVec.push_back(weights_arr_input_ ## X);
    //outBufVec.push_back(whd_output_ ## X);
    outBufVec.push_back(new_ref_idx_output_ ## X);

    //Copy input data to device global memory
    q.enqueueMigrateMemObjects(inBufVec, 0/* 0 means from host*/);

    //Set the Kernel Arguments
    krnl_indel.setArg(narg++, con_arr_input_ ## X);
    krnl_indel.setArg(narg++, con_size);
    krnl_indel.setArg(narg++, con_len_input_ ## X);
    krnl_indel.setArg(narg++, reads_arr_input_ ## X);
    krnl_indel.setArg(narg++, reads_size);
    krnl_indel.setArg(narg++, reads_len_input_ ## X);
    krnl_indel.setArg(narg++, weights_arr_input_ ## X);
    //krnl_indel.setArg(narg++, whd_output_ ## X);
    krnl_indel.setArg(narg++, new_ref_idx_output_ ## X);

}


#define COMPARE_RESULTS(X) {



}


#endif
