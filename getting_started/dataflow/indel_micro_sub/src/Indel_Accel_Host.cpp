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

    //int line_num = 0;
    //int base = 0;
    int* file_index = NULL;
    int size = 0;
    while ((read = getline(&line, &line_len, fp)) != -1) {
        //printf("Retrieved line of length %zu :\n", read);
        //printf("%s", line);
        printf("String length: %d\n", line_len);
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

