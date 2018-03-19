#include <hls_stream.h>
#include <stdio.h>
#include <string.h>
#include <ap_int.h>

#include "Indel_Accel.h"
//#define CON_SIZE 2
//#define READS_SIZE 8

// JENNY TODO
// ap_uint<3> for ATGCU 
// ap_uint<512> to round it to power of 2 
// ap_int<512> is the most efficient for global memory 
// Assign reads and consensus to different bundles for memory banking optimizations



void Indel_Rank (const int consensus_size, const int reads_size, int*  min_whd, int* __restrict new_ref_idx) {
    int new_ref[READS_SIZE];

    int min_score = 0x7fffffff;
    int min_idx = consensus_size + 1;
    int i, j;
    score: for (i = 1; i < consensus_size; i++) {
        int score = 0;
        for (j = 0; j < reads_size; j++) {
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

extern "C" {
//whd(con_arr, con_size, con_len, reads_arr, reads_size, reads_len, weights_arr, min_whd, min_whd_idx);
void Indel_Accel_Krnl (ap_uint<8>* consensus, const int consensus_size, int* consensus_length, \
    ap_uint<8>* reads, const int reads_size, int* reads_length, char* qs, int* new_ref_idx){
#pragma HLS INLINE
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


void Indel_Accel (ap_uint<8>* consensus_0, const int consensus_size_0, int* consensus_length_0, \
    ap_uint<8>* reads_0, const int reads_size_0, int* reads_length_0, char* qs_0, int* new_ref_idx_0, \
    ap_uint<8>* consensus_1, const int consensus_size_1, int* consensus_length_1, \
    ap_uint<8>* reads_1, const int reads_size_1, int* reads_length_1, char* qs_1, int* new_ref_idx_1, \
    ap_uint<8>* consensus_2, const int consensus_size_2, int* consensus_length_2, \
    ap_uint<8>* reads_2, const int reads_size_2, int* reads_length_2, char* qs_2, int* new_ref_idx_2, \ 
    ap_uint<8>* consensus_3, const int consensus_size_3, int* consensus_length_3, \
    ap_uint<8>* reads_3, const int reads_size_3, int* reads_length_3, char* qs_3, int* new_ref_idx_3, int global_id, int global_threads) {
 
 
    //ap_uint<8>* reads, const int reads_size, int* reads_length, char* qs, int* new_ref_idx, int* new_ref_idx, int global_id, int global_threads) {
#pragma HLS INTERFACE m_axi port=consensus_0 offset=slave bundle=gmem0
#pragma HLS INTERFACE m_axi port=consensus_length_0 offset=slave bundle=gmem0
#pragma HLS INTERFACE m_axi port=reads_0 offset=slave bundle=gmem0
#pragma HLS INTERFACE m_axi port=reads_length_0 offset=slave bundle=gmem0
#pragma HLS INTERFACE m_axi port=qs_0 offset=slave bundle=gmem0
#pragma HLS INTERFACE m_axi port=new_ref_idx_0 offset=slave bundle=gmem0
#pragma HLS INTERFACE s_axilite port=consensus_0 bundle=control
#pragma HLS INTERFACE s_axilite port=consensus_size_0 bundle=control
#pragma HLS INTERFACE s_axilite port=consensus_length_0 bundle=control
#pragma HLS INTERFACE s_axilite port=reads_0 bundle=control
#pragma HLS INTERFACE s_axilite port=reads_size_0 bundle=control
#pragma HLS INTERFACE s_axilite port=reads_length_0 bundle=control
#pragma HLS INTERFACE s_axilite port=qs_0 bundle=control
#pragma HLS INTERFACE s_axilite port=new_ref_idx_0 bundle=control

#pragma HLS INTERFACE m_axi port=consensus_1 offset=slave bundle=gmem0
#pragma HLS INTERFACE m_axi port=consensus_length_1 offset=slave bundle=gmem0
#pragma HLS INTERFACE m_axi port=reads_1 offset=slave bundle=gmem0
#pragma HLS INTERFACE m_axi port=reads_length_1 offset=slave bundle=gmem0
#pragma HLS INTERFACE m_axi port=qs_1 offset=slave bundle=gmem0
#pragma HLS INTERFACE m_axi port=new_ref_idx_1 offset=slave bundle=gmem0
#pragma HLS INTERFACE s_axilite port=consensus_1 bundle=control
#pragma HLS INTERFACE s_axilite port=consensus_size_1 bundle=control
#pragma HLS INTERFACE s_axilite port=consensus_length_1 bundle=control
#pragma HLS INTERFACE s_axilite port=reads_1 bundle=control
#pragma HLS INTERFACE s_axilite port=reads_size_1 bundle=control
#pragma HLS INTERFACE s_axilite port=reads_length_1 bundle=control
#pragma HLS INTERFACE s_axilite port=qs_1 bundle=control
#pragma HLS INTERFACE s_axilite port=new_ref_idx_1 bundle=control

#pragma HLS INTERFACE m_axi port=consensus_2 offset=slave bundle=gmem0
#pragma HLS INTERFACE m_axi port=consensus_length_2 offset=slave bundle=gmem0
#pragma HLS INTERFACE m_axi port=reads_2 offset=slave bundle=gmem0
#pragma HLS INTERFACE m_axi port=reads_length_2 offset=slave bundle=gmem0
#pragma HLS INTERFACE m_axi port=qs_2 offset=slave bundle=gmem0
#pragma HLS INTERFACE m_axi port=new_ref_idx_2 offset=slave bundle=gmem0
#pragma HLS INTERFACE s_axilite port=consensus_2 bundle=control
#pragma HLS INTERFACE s_axilite port=consensus_size_2 bundle=control
#pragma HLS INTERFACE s_axilite port=consensus_length_2 bundle=control
#pragma HLS INTERFACE s_axilite port=reads_2 bundle=control
#pragma HLS INTERFACE s_axilite port=reads_size_2 bundle=control
#pragma HLS INTERFACE s_axilite port=reads_length_2 bundle=control
#pragma HLS INTERFACE s_axilite port=qs_2 bundle=control
#pragma HLS INTERFACE s_axilite port=new_ref_idx_2 bundle=control

#pragma HLS INTERFACE m_axi port=consensus_3 offset=slave bundle=gmem0
#pragma HLS INTERFACE m_axi port=consensus_length_3 offset=slave bundle=gmem0
#pragma HLS INTERFACE m_axi port=reads_3 offset=slave bundle=gmem0
#pragma HLS INTERFACE m_axi port=reads_length_3 offset=slave bundle=gmem0
#pragma HLS INTERFACE m_axi port=qs_3 offset=slave bundle=gmem0
#pragma HLS INTERFACE m_axi port=new_ref_idx_3 offset=slave bundle=gmem0
#pragma HLS INTERFACE s_axilite port=consensus_3 bundle=control
#pragma HLS INTERFACE s_axilite port=consensus_size_3 bundle=control
#pragma HLS INTERFACE s_axilite port=consensus_length_3 bundle=control
#pragma HLS INTERFACE s_axilite port=reads_3 bundle=control
#pragma HLS INTERFACE s_axilite port=reads_size_3 bundle=control
#pragma HLS INTERFACE s_axilite port=reads_length_3 bundle=control
#pragma HLS INTERFACE s_axilite port=qs_3 bundle=control
#pragma HLS INTERFACE s_axilite port=new_ref_idx_3 bundle=control
#pragma HLS INTERFACE s_axilite port=global_id bundle=control
#pragma HLS INTERFACE s_axilite port=global_threads bundle=control
#pragma HLS INTERFACE s_axilite port=return bundle=control
//#pragma HLS INTERFACE m_axi port=new_ref_idx offset=slave bundle=gmem2
//#pragma HLS INTERFACE s_axilite port=new_ref_idx bundle=control
 

printf("DEBUG:\n");
#pragma HLS expression_balance
Indel_Accel_Krnl(consensus_0, consensus_size_0, consensus_length_0, reads_0, reads_size_0, reads_length_0, qs_0, new_ref_idx_0);
Indel_Accel_Krnl(consensus_1, consensus_size_1, consensus_length_1, reads_1, reads_size_1, reads_length_1, qs_1, new_ref_idx_1);
Indel_Accel_Krnl(consensus_2, consensus_size_2, consensus_length_2, reads_2, reads_size_2, reads_length_2, qs_2, new_ref_idx_2);
Indel_Accel_Krnl(consensus_3, consensus_size_3, consensus_length_3, reads_3, reads_size_3, reads_length_3, qs_3, new_ref_idx_3);
//ap_uint<8>* consensus;
//const int consensus_size;
//int* consensus_length;
//ap_uint<8>* reads;
//const int reads_size;
//int* reads_length;
//char* qs;
//int* new_ref_idx;

printf("FINISH:\n");
}
}
