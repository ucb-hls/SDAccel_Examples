void whd (char* consensus, const int consensus_size, int* consensus_length, \
    char* reads, const int reads_size, int* reads_length, char* qs, \
    int* min_whd, int* min_whd_idx) {
    int i, j, k, l;
    for (i = 0; i < consensus_size; i++) {
        int consensus_base =  consensus_length[i];
        int local_consensus_length =  consensus_length[i+1] - consensus_length[i];
        for (j = 0; j < reads_size; j++) {
            int reads_base = reads_length[j];
            int local_reads_length = reads_length[j+1]-reads_length[j];
            fprintf(stderr, "consensus size %d i %d consensus length %d, read size %d j %d reads length %d\n", \
                consensus_size, i, local_consensus_length, reads_size,  j, local_reads_length);
            int min = 0x7fffffff; 
            int min_idx = local_consensus_length - local_reads_length + 1;
            for (k = 0; k <= local_consensus_length - local_reads_length; k++) {

                // whd 
                int whd = 0;
                // Optimization tree based reduction
                for (l = 0; l < local_reads_length; l++) {

                    //printf("%c", consensus[consensus_base + k + l]);
                    //printf("%c", reads[reads_base + k + l]);
                    if (consensus[consensus_base + k + l] != reads[reads_base + l]){
                        whd += qs[reads_base + l];
                        //if(k == 8 & j == 1){
                        //    printf("whd: %d\t", whd);
                        //}
                    }                        
                    //printf("whd: %d\t", whd);
                    //printf("\t");

                }

                //printf("\n");
                if (whd < min) {
                    min =  whd; 
                    min_idx = k; 
               }

            }
            fprintf(stderr, "min_idx %d, min %d\n", min_idx, min);
            assert(min_idx <= local_consensus_length - local_reads_length);
            
            min_whd[i * reads_size + j] = min;
            min_whd_idx[i * reads_size + j] = min_idx;
            //reads_base += local_reads_length;
        }
        //consensus_base += local_consensus_length;
    }
}



void whd_vec (ap_uint<4>* consensus, const int consensus_size, int* consensus_length, \
    ap_uint<4>* reads, const int reads_size, int* reads_length, char* qs, \
    int* min_whd, int* min_whd_idx) {
    int i, j, k, l;
    for (i = 0; i < consensus_size; i++) {
        int consensus_base =  consensus_length[i];
        int local_consensus_length =  consensus_length[i+1] - consensus_length[i];
        for (j = 0; j < reads_size; j++) {
            int reads_base = reads_length[j];
            int local_reads_length = reads_length[j+1]-reads_length[j];
            fprintf(stderr, "consensus size %d i %d consensus length %d, read size %d j %d reads length %d\n", \
                consensus_size, i, local_consensus_length, reads_size,  j, local_reads_length);
            int min = 0x7fffffff; 
            int min_idx = local_consensus_length - local_reads_length + 1;
            for (k = 0; k <= local_consensus_length - local_reads_length; k++) {

                // whd 
                int whd = 0;
                // Optimization tree based reduction
                for (l = 0; l < local_reads_length; l++) {

                    //printf("%c", consensus[consensus_base + k + l]);
                    //printf("%c", reads[reads_base + k + l]);
                    if (consensus[consensus_base + k + l] != reads[reads_base + l]){
                        whd += qs[reads_base + l];
                        //if(k == 8 & j == 1){
                        //    printf("whd: %d\t", whd);
                        //}
                    }                        
                    //printf("\t");

                }
                //printf("whd: %d\t", whd);
                //printf("\n");
                if (whd < min) {
                    min =  whd; 
                    min_idx = k; 
               }

            }
            fprintf(stderr, "min_idx %d, min %d\n", min_idx, min);
            assert(min_idx <= local_consensus_length - local_reads_length);

            min_whd[i * reads_size + j] = min;
            min_whd_idx[i * reads_size + j] = min_idx;
            //reads_base += local_reads_length;
        }
        //consensus_base += local_consensus_length;
    }
}

void print_whd (int* min_whd, int consensus_size, int reads_size) {
    int i, j;
    fprintf(stderr, "whd: \n");
    for (i = 0; i < reads_size; i ++) {
        for (j = 0; j < consensus_size; j ++) {
            int tmp = min_whd[j * reads_size + i];
            fprintf(stderr, "%d\t", tmp);
        }
        fprintf(stderr, "\n");
    }
}

void score_whd (int* min_whd, int* min_whd_idx, int consensus_size, int reads_size, int* new_ref, int* new_ref_idx) {
    //print_whd(min_whd, consensus_size, reads_size);
    // might need to reduce the bits used in here 
    assert(new_ref_idx != NULL); 
    assert(min_whd != NULL); 
    assert(min_whd_idx != NULL); 
    int min_score = 0x7fffffff;
    int min_idx = consensus_size + 1;
    int i, j;
    for (i = 1; i < consensus_size; i++) {
        int score = 0;
        for (j = 0; j < reads_size; j++) {
            int tmp = min_whd[i * reads_size + j] - min_whd[j];
            score += (tmp > 0) ? tmp: -tmp;
        }
        min_idx = (score < min_score) ? i : min_idx;
        //scores[i] = score;
    }
    fprintf(stderr, "min_idx: %d\n", min_idx);
    assert(min_idx < consensus_size);
    
    for (j = 0; j < reads_size; j++) {
        //if ( min_whd[ min_idx * reads_size + j] < min_whd[j]){
            new_ref[j] = min_whd[min_idx * reads_size + j];
            new_ref_idx[j] = min_whd_idx[ min_idx * reads_size +j];
        //} else{
        //    new_ref[j] = min_whd [j];
        //    new_ref_idx[j] = min_whd_idx[j]; 
        //}
    }

    for (j = 0; j < reads_size; j++) {
        printf("Read %2d whd %4d  index %2d\n", j, new_ref[j], new_ref_idx[j]);
    }
}


void whd_ref (ap_uint<4>* consensus, const int consensus_size, int* consensus_length, \
    ap_uint<4>* reads, const int reads_size, int* reads_length, char* qs, \
    int* new_ref, int* new_ref_idx) {

    int* min_whd = (int*) malloc(consensus_size* reads_size * sizeof(int));
    int* min_whd_idx = (int*) malloc(consensus_size* reads_size * sizeof(int));

    whd_vec(consensus, consensus_size, consensus_length, reads, reads_size, reads_length, qs, min_whd, min_whd_idx); 
    score_whd(min_whd, min_whd_idx, consensus_size, reads_size, new_ref , new_ref_idx);
    free(min_whd);
    free(min_whd_idx);
}


void generate_ref(std::vector<ap_uint<4>, aligned_allocator<ap_uint<4>>>* consensus, std::vector<int, aligned_allocator<int>>* con_size, std::vector<int, aligned_allocator<int>>* con_base, \
    std::vector<ap_uint<4>, aligned_allocator<ap_uint<4>>>* reads, std::vector<int, aligned_allocator<int>>* reads_size, std::vector<int, aligned_allocator<int>>* reads_base, std::vector<char, aligned_allocator<char>>* qs, std::vector<int, aligned_allocator<int>>* new_ref_idx) {

    // Since size is actual size + 1 [0, base]
    for (size_t i = 0; i < reads_size->size() - 1; i++) {
        int con_size_base = (*con_size)[i];
        int con_size_local = (*con_size)[i + 1] - con_size_base;

        int reads_size_base = (*reads_size)[i];
        int rs_size_local = (*reads_size)[i + 1] - reads_size_base;

        printf("TEST %d base: %d\n", i, reads_size_base);
        int* new_ref = (int*) malloc(rs_size_local* sizeof(int));
        int* new_ref_idx_ref = (int*) malloc(rs_size_local* sizeof(int));
        whd_ref(consensus->data(), con_size_local, &(con_base->data())[con_size_base], reads->data(), rs_size_local, &(reads_base->data()[reads_size_base]), qs->data(), new_ref, new_ref_idx_ref);

        std::vector<int> new_ref_vec(new_ref, new_ref + rs_size_local);
        std::vector<int> new_ref_idx_vec(new_ref_idx_ref, new_ref_idx_ref + rs_size_local);

        new_ref_idx->insert(new_ref_idx->end(), new_ref_idx_vec.begin(), new_ref_idx_vec.end());
        free(new_ref);
        free(new_ref_idx_ref);
     }
}

void compare_results(std::vector<int, aligned_allocator<int>>*new_ref_idx_buffer, std::vector<int, aligned_allocator<int>>* new_ref_idx_ref_buffer, std::vector<int, aligned_allocator<int>>* reads_size){

    for (size_t i = 0; i < reads_size->size() - 1; i++) {
        int match = 0;
        std::cout << "TEST " << i << std::endl; 
        for (int j = (*reads_size)[i]; j < (*reads_size)[i+1] ; j++) {
            if ((*new_ref_idx_ref_buffer)[j] != (*new_ref_idx_buffer)[j]){
                std::cout << "Error: Result mismatch" << std::endl;
                std::cout << "i = " << j << " CPU result = " << (*new_ref_idx_ref_buffer)[j]
                    << " Device result = " << (*new_ref_idx_buffer)[j] << std::endl;
                match = 1;
            }
        }
        std::cout << "TEST " << (match ? "FAILED" : "PASSED") << std::endl;
    } 
}
