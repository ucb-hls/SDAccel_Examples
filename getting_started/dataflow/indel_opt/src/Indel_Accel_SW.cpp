
void whd (char* consensus, const int consensus_size, int* consensus_length, \
    char* reads, const int reads_size, int* reads_length, char* qs, \
    int* min_whd, int* min_whd_idx) {
    int consensus_base = 0; 
    int i, j, k, l;
    for (i = 0; i < consensus_size; i++) {
        int local_consensus_length =  consensus_length[i];
        int reads_base = 0;
        for (j = 0; j < reads_size; j++) {
            int local_reads_length = reads_length[j];
            fprintf(stderr, "consensus size %d i %d consensus length %d, read size %d j %d reads length %d\n", \
                consensus_size, i, local_consensus_length, reads_size,  j, local_reads_length);
            int min = 0x7fffffff; 
            int min_idx = local_consensus_length - local_reads_length + 1;
            for (k = 0; k <= local_consensus_length - local_reads_length; k++) {

                // whd 
                int whd = 0;
                // Optimization tree based reduction
                for (l = 0; l < reads_length[j]; l++) {

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
            fprintf(stderr, "min_idx %d\n", min_idx);
            assert(min_idx <= local_consensus_length - local_reads_length);
            
            min_whd[i * reads_size + j] = min;
            min_whd_idx[i * reads_size + j] = min_idx;
            reads_base += local_reads_length;
        }
        consensus_base += local_consensus_length;
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


