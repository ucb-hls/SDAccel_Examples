#include <hls_stream.h>
#include <stdio.h>
#include <string.h>
#include <ap_int.h>

#include "Indel_Accel.h"
//#define CON_SIZE 2
//#define READS_SIZE 8

extern "C" {
//whd(con_arr, con_size, con_len, reads_arr, reads_size, reads_len, weights_arr, min_whd, min_whd_idx);
void Indel_Accel (char* consensus, const int consensus_size, int* consensus_length, \
    char* reads, const int reads_size, int* reads_length, char* qs, int* new_ref_idx) {
#pragma HLS INTERFACE m_axi port=consensus offset=slave bundle=gmem
#pragma HLS INTERFACE m_axi port=consensus_length offset=slave bundle=gmem
#pragma HLS INTERFACE m_axi port=reads offset=slave bundle=gmem
#pragma HLS INTERFACE m_axi port=reads_length offset=slave bundle=gmem
#pragma HLS INTERFACE m_axi port=qs offset=slave bundle=gmem
#pragma HLS INTERFACE m_axi port=new_ref_idx offset=slave bundle=gmem

#pragma HLS INTERFACE s_axilite port=consensus bundle=control
#pragma HLS INTERFACE s_axilite port=consensus_size bundle=control
#pragma HLS INTERFACE s_axilite port=consensus_length bundle=control
#pragma HLS INTERFACE s_axilite port=reads bundle=control
#pragma HLS INTERFACE s_axilite port=reads_size bundle=control
#pragma HLS INTERFACE s_axilite port=reads_length bundle=control
#pragma HLS INTERFACE s_axilite port=qs bundle=control
#pragma HLS INTERFACE s_axilite port=new_ref_idx bundle=control
#pragma HLS INTERFACE s_axilite port=return bundle=control


    //array of stream declaration
    /*hls::stream<char> condStreams[CON_SIZE];
    //#pragma HLS STREAM variable=condStreams depth=32

    hls::stream<char> readStreams[READS_SIZE];
    //#pragma HLS STREAM variable=readStreams depth=32


    #pragma HLS dataflow
    for (int i = 0 ; i < CON_SIZE; i ++){
        read_input(consensus, condStreams[i], consensus_length[i]);
    }
   
    #pragma HLS dataflow
    for (int i = 0 ; i < READS_SIZE; i ++){
        read_input(reads, readStreams[i], reads_length[i]);
    }*/
    int min_whd[CON_SIZE * READS_SIZE];
    int min_whd_idx[CON_SIZE * READS_SIZE];
    int new_ref[READS_SIZE];
    int consensus_base = 0; 
    int i, j, k, l;
    for (i = 0; i < consensus_size; i++) {
        int local_consensus_length =  consensus_length[i];
        int reads_base = 0;
        for (j = 0; j < reads_size; j++) {
            int local_reads_length = reads_length[j];
            //printf( "consensus size %d i %d consensus length %d, read size %d j %d reads length %d\n", \
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
            //printf( "min_idx %d\n", min_idx);
            //assert(min_idx <= local_consensus_length - local_reads_length);
            
            //min_whd[i][j] = min;
            //min_whd_idx[i][j] = min_idx;
            min_whd[i * reads_size + j] = min;
            min_whd_idx[i * reads_size + j] = min_idx;
 
            reads_base += local_reads_length;
        }
        consensus_base += local_consensus_length;
    }
    
    int min_score = 0x7fffffff;
    int min_idx = consensus_size + 1;
    for (i = 1; i < consensus_size; i++) {
        int score = 0;
        for (j = 0; j < reads_size; j++) {
            int tmp = min_whd[i * reads_size + j] - min_whd[j];
            score += (tmp > 0) ? tmp: -tmp;
        }
        min_idx = (score < min_score) ? i : min_idx;
        //scores[i] = score;
    }
    //printf( "min_idx: %d\n", min_idx);
    //assert(min_idx < consensus_size);
    
    for (j = 0; j < reads_size; j++) {
        //if ( min_whd[ min_idx * reads_size + j] < min_whd[j]){
            //new_ref[j] = min_whd[min_idx][j];
            //new_ref_idx[j] = min_whd_idx[min_idx][j];
            new_ref[j] = min_whd[min_idx * reads_size + j];
            new_ref_idx[j] = min_whd_idx[ min_idx * reads_size +j];
 
        //} else{
        //    new_ref[j] = min_whd [j];
        //    new_ref_idx[j] = min_whd_idx[j]; 
        //}
    }
    for (j = 0; j < reads_size; j++) {
     printf("Read %2d whd %2d index %2d\n", j, new_ref[j], new_ref_idx[j]);
    }

}

}

