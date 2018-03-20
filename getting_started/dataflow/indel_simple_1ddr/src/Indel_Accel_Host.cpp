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
        char* tmp = (char*) aligned_alloc(sizeof(char), line_len * sizeof(char));
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
        file_index = (int*) aligned_alloc(sizeof(int), size * sizeof(int));

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
    int * con_len = (int*) aligned_alloc(sizeof(int), sizeof(int) * num_lines);
    assert(con_len != NULL);
    //* con_len_arr = con_len;
    // Location of the base of each segment plus the end 
    int* pos = (int*) aligned_alloc(sizeof(int), sizeof(int) * (num_lines + 1));
    * con_len_arr = pos;
    assert(pos != NULL);

    char separators[] = "|";
    char *line = NULL;
    size_t line_len = 0;
    ssize_t read;

    int line_num = 0;
    int base = 0;
    pos[0] = 0;
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
                //  con_arr = (char*) aligned_alloc(sizeof(char), len * num_lines * sizeof(char));
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
        pos[line_num+1] = base; 
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

//#define CREATE_INPUT(X) {
void Run_Unit(int X, int test_idx, cl::Context &context, cl::CommandQueue &q, cl::Kernel &krnl_indel,\ 
     std::vector<cl::Memory> &inBufVec, std::vector<cl::Memory> &outBufVec, int * &new_ref_idx_X, \ 
     int * &new_ref_idx_ref, std::chrono::milliseconds* parse_time, int &narg, int num_tests, int* test_indices, int &reads_size_X) {

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
 
    int test_idx_X= test_idx + X;

    printf("test_idx_X: %d\n", test_idx_X);
    //int reads_size_X;

    if (test_idx_X < num_tests){
        snprintf(test_num, sizeof(test_num), "%d", test_indices[test_idx_X]);
    
        strcat(con, file_prefix);
        strcat(con, test_num);
        strcat(con, ".SEQ");

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
        parse( reads, 4, reads_arr, &reads_len, &reads_size_X , &reads_total_len);
        parse( reads, 5, weights_arr, &reads_len, &reads_size_X, &weights_total_len);

        finish = std::chrono::high_resolution_clock::now();
        //std::chrono::duration<double> elapsed = finish - start;   
        duration = std::chrono::duration_cast<std::chrono::milliseconds>(finish - start);
        parse_time[X] = duration;
        //std::cout << "Parsing time is :" << duration.count() << " ms\n";
        //printf("Parse time is: %0.3f ms\n", duration.count());

        min_whd = (int*) malloc(con_size * reads_size_X * sizeof(int));
        min_whd_idx = (int*) malloc(con_size * reads_size_X * sizeof(int));
        new_ref = (int*) malloc(reads_size_X * sizeof(int));
        //int* new_ref_idx = (int*) malloc(reads_size_X * sizeof(int));
        new_ref_idx_ref = (int*) malloc(reads_size_X * sizeof(int));

        whd(con_arr, con_size, con_len, reads_arr, reads_size_X, reads_len, weights_arr, min_whd, min_whd_idx);
        score_whd (min_whd,  min_whd_idx, con_size, reads_size_X, new_ref, new_ref_idx_ref);

        printf("Software Success!\n");
    } else {
        con_size=reads_size_X=0;
        con_total_len=reads_total_len=0;
        con_arr=reads_arr=weights_arr=NULL;
        con_len=reads_len=NULL;
        min_whd=min_whd_idx=new_ref=new_ref_idx_ref=NULL;   
    }

    //return 0;
    //Allocate Memory in Host Memory
    //std::vector<char,aligned_allocator<char>> con_arr_buffer     ( con_arr, con_arr + CON_SIZE * CON_LEN);
    //std::vector<char,aligned_allocator<char>> reads_arr_buffer     (reads_arr, reads_arr + READS_SIZE * READS_LEN);
    //std::vector<char,aligned_allocator<char>> weights_arr_buffer     (weights_arr, weights_arr + READS_SIZE * READS_LEN);

    start = std::chrono::high_resolution_clock::now();
    std::vector<ap_uint<4>,aligned_allocator<ap_uint<4>>> con_arr_buffer_X     ( CON_SIZE * CON_LEN >> 1);
    //std::vector<char,aligned_allocator<char>> con_arr_buffer_X      ( con_size * con_total_len >> 1);
    std::vector<ap_uint<4>,aligned_allocator<ap_uint<4>>> reads_arr_buffer_X     ( READS_SIZE * READS_LEN >> 1);
    //std::vector<char,aligned_allocator<char>> reads_arr_buffer_X    ( reads_size_X * reads_total_len >> 1);
    std::vector<char,aligned_allocator<char>> weights_arr_buffer_X  (weights_arr, weights_arr + READS_SIZE * READS_LEN);
    //std::vector<char,aligned_allocator<char>> weights_arr_buffer_X  (weights_arr, weights_arr + reads_size_X* reads_total_len);

  //  for(int i = 0 ; i < CON_SIZE * CON_LEN; i++){
  //      con_arr_buffer[i] = m[con_arr[i]];
  //  }

    //for(int i = 0 ; i < CON_SIZE * CON_LEN; i++){
    for(int i = 0 ; i < con_size*con_total_len; i++){
        int idx = i >> 1;
        int offset = (i % 2) << 2;
        char c = m[con_arr[i]] & 0xf;
        con_arr_buffer_X[idx] = (con_arr_buffer_X[idx] & ~(0xf << offset)) | (c << offset);
        //printf("con_arr_buffer_X[%d + %d]: %hhX \n", idx, offset, con_arr_buffer_X[idx]);
    }

 //   printf("Read Buffer:");     
 //   for(int i = 0 ; i < READS_LEN * READS_SIZE; i++){
 //       reads_arr_buffer[i] = m[reads_arr[i]];
 //       //unsigned char print_var = reads_arr_buffer[i];
 //       //printf("%x", print_var);     
 //   }
    for(int i = 0 ; i < reads_size_X*reads_total_len; i++){
        int idx = i >> 1;
        int offset = (i % 2) << 2;
        char c = m[reads_arr[i]] & 0xf;
        reads_arr_buffer_X[idx] = (reads_arr_buffer_X[idx] & ~(0xf << offset)) | (c << offset);
        //printf("reads_arr_buffer_X[%d]: %hhX \n", idx, (char)reads_arr_buffer_X[idx]);
    }

    //std::vector<int,aligned_allocator<int>> con_len_buffer_X     (con_len, con_len + con_size);
    std::vector<int,aligned_allocator<int>> con_len_buffer_X     (con_len, con_len + CON_SIZE);
    //std::vector<int,aligned_allocator<int>> reads_len_buffer_X     (reads_len, reads_len + reads_size_X);
    std::vector<int,aligned_allocator<int>> reads_len_buffer_X     (reads_len, reads_len + READS_SIZE);

    printf("narg: %d\n", narg);
    printf("reads_size_X: %d\n", reads_size_X);
    printf("con_size: %d\n", con_size );
    printf("con_len: \n");
    for(int i =0; i < con_size; i++) {
      printf("%d\t", con_len_buffer_X[i]);
    }
    printf("\n");

    //std::vector<int,aligned_allocator<int>> whd_buffer(con_size * reads_size_X << 1);
    new_ref_idx_X = (int*) malloc(sizeof(int) * reads_size_X);
    for(int i = 0 ; i < reads_size_X; i++){
        new_ref_idx_X[i] = 0;
    }
    printf("\n");

    //std::vector<int,aligned_allocator<int>> new_ref_idx(new_ref_idx_X, new_ref_idx_X + reads_size_X);
    std::vector<int,aligned_allocator<int>> new_ref_idx(new_ref_idx_X, new_ref_idx_X + READS_SIZE);
    //int * new_ref_idx = new int [reads_size_X];

    finish = std::chrono::high_resolution_clock::now();

    duration = std::chrono::duration_cast<std::chrono::milliseconds>(finish - start);
    std::cout << "Preprocess time is : " << duration.count() << " s\n";
    //printf("Preprocess time is: %0.3f ms\n", duration.count());


    cl_mem_ext_ptr_t con_arr_buffer_ptr_X, reads_arr_buffer_ptr_X, weights_arr_buffer_ptr_X, con_len_buffer_ptr_X, reads_len_buffer_ptr_X, whd_buffer_ptr_X,  new_ref_idx_ptr_X; 
    con_arr_buffer_ptr_X.flags  = XCL_MEM_DDR_BANK0; 
    con_len_buffer_ptr_X.flags  = XCL_MEM_DDR_BANK0; 
    reads_arr_buffer_ptr_X.flags  = XCL_MEM_DDR_BANK0; 
    reads_len_buffer_ptr_X.flags  = XCL_MEM_DDR_BANK0; 
    weights_arr_buffer_ptr_X.flags  = XCL_MEM_DDR_BANK0; 
    new_ref_idx_ptr_X.flags  = XCL_MEM_DDR_BANK0; 
    //whd_buffer_ptr_X.flags = XCL_MEM_DDR_BANK3;
 
    std::cout << "Preprocess time is : " << duration.count() << " ms\n";

    // Setting input and output objects
    con_arr_buffer_ptr_X.obj = con_arr_buffer_X.data();
    //con_len_buffer_ptr_X.obj = con_len;  // this would cause extra memcpy since it is not aligned
    con_len_buffer_ptr_X.obj = con_len_buffer_X.data();
    reads_arr_buffer_ptr_X.obj = reads_arr_buffer_X.data();
    //reads_len_buffer_ptr_X.obj = reads_len;
    reads_len_buffer_ptr_X.obj = reads_len_buffer_X.data();
    weights_arr_buffer_ptr_X.obj = weights_arr_buffer_X.data();
    //whd_buffer_ptr_X.obj = whd_buffer_X.data();
    new_ref_idx_ptr_X.obj = new_ref_idx.data();

    std::cout << "Preprocess time is : " << duration.count() << " ms\n";
    // Setting param to zero 
    con_arr_buffer_ptr_X.param = 0; con_len_buffer_ptr_X.param = 0;
    reads_arr_buffer_ptr_X.param = 0; reads_len_buffer_ptr_X.param = 0; 
    weights_arr_buffer_ptr_X.param = 0; new_ref_idx_ptr_X.param = 0;
    whd_buffer_ptr_X.param = 0;

    printf("Finish creating buffers\n");
    //Allocate Buffer in Global Memory

    cl::Buffer con_arr_input_X(context,CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY | CL_MEM_EXT_PTR_XILINX ,\ 
            CON_SIZE * CON_LEN, &con_arr_buffer_ptr_X);
    cl::Buffer reads_arr_input_X(context,CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY | CL_MEM_EXT_PTR_XILINX ,\ 
            READS_SIZE * READS_LEN, &reads_arr_buffer_ptr_X);
    cl::Buffer weights_arr_input_X(context,CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY | CL_MEM_EXT_PTR_XILINX ,\ 
            READS_SIZE * READS_LEN, &weights_arr_buffer_ptr_X);

    std::vector<cl::Memory> con_len_vec, reads_len_vec;
    cl::Buffer con_len_input_X(context,CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY | CL_MEM_EXT_PTR_XILINX , 
            CON_SIZE * sizeof(int), &con_len_buffer_ptr_X);
    cl::Buffer reads_len_input_X(context,CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY | CL_MEM_EXT_PTR_XILINX , 
            READS_SIZE * sizeof(int), &reads_len_buffer_ptr_X);

    //cl::Buffer whd_output_X(context, CL_MEM_USE_HOST_PTR | CL_MEM_WRITE_ONLY| CL_MEM_EXT_PTR_XILINX ,
    //        reads_size_X * con_size * 2 * sizeof(int), &whd_buffer_ptr_X);
    cl::Buffer new_ref_idx_output_X(context, CL_MEM_USE_HOST_PTR | CL_MEM_WRITE_ONLY | CL_MEM_EXT_PTR_XILINX ,
            READS_SIZE * sizeof(int), &new_ref_idx_ptr_X);

    inBufVec.push_back(con_len_input_X);
    inBufVec.push_back(reads_len_input_X);
    inBufVec.push_back(con_arr_input_X);
    inBufVec.push_back(reads_arr_input_X);
    inBufVec.push_back(weights_arr_input_X);
    //outBufVec.push_back(whd_output_X);
    outBufVec.push_back(new_ref_idx_output_X);

    //Set the Kernel Arguments
    krnl_indel.setArg(narg++, con_arr_input_X);
    krnl_indel.setArg(narg++, con_size);
    krnl_indel.setArg(narg++, con_len_input_X);
    krnl_indel.setArg(narg++, reads_arr_input_X);
    krnl_indel.setArg(narg++, reads_size_X);
    krnl_indel.setArg(narg++, reads_len_input_X);
    krnl_indel.setArg(narg++, weights_arr_input_X);
    //krnl_indel.setArg(narg++, whd_output_X);
    krnl_indel.setArg(narg++, new_ref_idx_output_X);
   
}


//#define Compare_Results(X) {
void Compare_Results(int X, int reads_size, std::vector<int, aligned_allocator<int> > new_ref_idx_X, int* & new_ref_idx_ref_X){
    // Compare the results of the Device to the simulation
    int match = 0;
    for (int i = 0 ; i < reads_size; i++){
        if (new_ref_idx_X [i] != new_ref_idx_ref_X [i]){
            std::cout << "Error: Result mismatch" << std::endl;
            std::cout << "i = " << i << " CPU result = " << new_ref_idx_ref_X[i]
                << " Device result = " << new_ref_idx_X[i] << std::endl;
            match = 1;
            //break;
        }
    }

    std::cout << "TEST " << (match ? "FAILED" : "PASSED") << std::endl; 
    //if (!match)
    //  return EXIT_FAILURE;
    //return (match ? EXIT_FAILURE :  EXIT_SUCCESS);

}

