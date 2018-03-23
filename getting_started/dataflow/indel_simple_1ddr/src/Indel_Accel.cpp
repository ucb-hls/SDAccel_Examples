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

extern "C" {
//whd(con_arr, con_size, con_len, reads_arr, reads_size, reads_len, weights_arr, min_whd, min_whd_idx);
void Indel_Accel_Krnl (ap_uint<4>* consensus, const int consensus_size, int* consensus_length, \
    ap_uint<4>* reads, const int reads_size, int* reads_length, char* qs, int* new_ref_idx, int new_ref_idx_base){

printf("DEBUG:\n");
//#pragma HLS INLINE
    //#pragma HLS DATAFLOW
    #pragma HLS expression_balance
    int min_whd[CON_SIZE * READS_SIZE];
    int min_whd_idx[CON_SIZE * READS_SIZE];
    int new_ref[READS_SIZE];
    //int i, j, k, l;
    for (int i = 0; i < consensus_size; i++) {

        int consensus_base = consensus_length[i];
        int local_consensus_length =  consensus_length[i+1] - consensus_length[i];
        for (int j = 0; j < reads_size; j++) {
        #pragma HLS unroll factor=4
            int reads_base = reads_length[j];
            int local_reads_length = reads_length[j+1] - reads_length[j];
            printf( "consensus size %d i %d consensus length %d, read size %d j %d reads length %d\n", \
                consensus_size, i, local_consensus_length, reads_size,  j, local_reads_length);
            int min = 0x7fffffff; 
            int min_idx = 0x7fffffff;

            printf("Length Diff: %d\n", local_consensus_length - local_reads_length);
            for (int k = 0; k <= local_consensus_length - local_reads_length; k++) {
            #pragma HLS loop_tripcount min=0 max=2048


                // whd 
                int whd = 0;
                // Optimization tree based reduction
                for (int l = 0; l < local_reads_length; l++) {
                #pragma HLS loop_tripcount min=0 max=256
                    char con_char = consensus[consensus_base + k + l];
                    con_char = con_char & 0xf;
                    char reads_char = reads[reads_base + k + l];
                    reads_char = reads_char & 0xf;
                    printf("reads_char: %x -- con_char: %x; ", reads_char, con_char);
                    printf("qs %d %d con %d %d,", reads_base, l,  consensus_base, k + l);
            
                    if (consensus[consensus_base + k + l] != reads[reads_base + l]){
                        whd += qs[reads_base + l];
                        //if(k == 8 & j == 1){
                        //    printf("whd: %d\t", whd);
                        //}
                    }                        
                    printf("whd: %d\t", whd);
                    printf("\t");

                }

                printf("\n");
                if (whd < min) {
                    min =  whd; 
                    min_idx = k; 
               }

            }
            //printf( "min_idx %d\n", min_idx);
            //assert(min_idx <= local_consensus_length - local_reads_length);
            printf( "[%d, %d]-Kernel:min_idx %d min %d\n", i, j, min_idx, min);
            
            //min_whd[i][j] = min;
            //min_whd_idx[i][j] = min_idx;
            min_whd[i * reads_size + j] = min;
            min_whd_idx[i * reads_size + j] = min_idx;
 
            //reads_base += local_reads_length;
        }
        //consensus_base += local_consensus_length;
    }
    
    int min_score = 0x7fffffff;
    int min_idx = consensus_size + 1;
    for (int i = 1; i < consensus_size; i++) {
        int score = 0;
        for (int j = 0; j < reads_size; j++) {
        #pragma HLS unroll factor=4
            int tmp = min_whd[i * reads_size + j] - min_whd[j];
            score += (tmp > 0) ? tmp: -tmp;
        }
        min_idx = (score < min_score) ? i : min_idx;
        //scores[i] = score;
    }
    //printf( "min_idx: %d\n", min_idx);
    //assert(min_idx < consensus_size);
    
    for (int j = 0; j < reads_size; j++) {
    #pragma HLS unroll factor=4
        //if ( min_whd[ min_idx * reads_size + j] < min_whd[j]){
            //new_ref[j] = min_whd[min_idx][j];
            //new_ref_idx[j] = min_whd_idx[min_idx][j];
            new_ref[new_ref_idx_base + j] = min_whd[min_idx * reads_size + j];
            new_ref_idx[new_ref_idx_base + j] = min_whd_idx[ min_idx * reads_size +j];
 
        //} else{
        //    new_ref[j] = min_whd [j];
        //    new_ref_idx[j] = min_whd_idx[j]; 
        //}
    }
    for (int j = 0; j < reads_size; j++) {
     printf("Read %2d whd %2d index %2d\n", j, new_ref[new_ref_idx_base + j], new_ref_idx[new_ref_idx_base + j]);
    }
}


void Indel_Accel (ap_uint<4>* consensus, int* consensus_size, int* consensus_length, \
    ap_uint<4>* reads, int* reads_size, int* reads_length, char* qs, int* new_ref_idx, \
    int global_id, int global_threads) {
 
    //ap_uint<4>* reads, const int reads_size, int* reads_length, char* qs, int* new_ref_idx, int* new_ref_idx, int global_id, int global_threads) {
#pragma HLS INTERFACE m_axi port=consensus offset=slave bundle=gmem0
#pragma HLS INTERFACE m_axi port=consensus_size offset=slave bundle=gmem0
#pragma HLS INTERFACE m_axi port=consensus_length offset=slave bundle=gmem0
#pragma HLS INTERFACE m_axi port=reads offset=slave bundle=gmem0
#pragma HLS INTERFACE m_axi port=reads_size offset=slave bundle=gmem0
#pragma HLS INTERFACE m_axi port=reads_length offset=slave bundle=gmem0
#pragma HLS INTERFACE m_axi port=qs offset=slave bundle=gmem0
#pragma HLS INTERFACE m_axi port=new_ref_idx offset=slave bundle=gmem0
#pragma HLS INTERFACE s_axilite port=consensus bundle=control
#pragma HLS INTERFACE s_axilite port=consensus_size bundle=control
#pragma HLS INTERFACE s_axilite port=consensus_length bundle=control
#pragma HLS INTERFACE s_axilite port=reads bundle=control
#pragma HLS INTERFACE s_axilite port=reads_size bundle=control
#pragma HLS INTERFACE s_axilite port=reads_length bundle=control
#pragma HLS INTERFACE s_axilite port=qs bundle=control
#pragma HLS INTERFACE s_axilite port=new_ref_idx bundle=control

#pragma HLS INTERFACE s_axilite port=global_id bundle=control
#pragma HLS INTERFACE s_axilite port=global_threads bundle=control
#pragma HLS INTERFACE s_axilite port=return bundle=control
//#pragma HLS INTERFACE m_axi port=new_ref_idx offset=slave bundle=gmem2
//#pragma HLS INTERFACE s_axilite port=new_ref_idx bundle=control
 

printf("DEBUG:\n");
#pragma HLS expression_balance
//#pragma HLS DATAFLOW
//Copy_Con(consensus, consensus_size, consensus_length);
// 0-2 2-4
// Number of lengths 
int con_size_base = consensus_size[global_id];
int con_size = consensus_size[global_id + 1] - con_size_base;

int reads_size_base = reads_size[global_id];
int rs_size = reads_size[global_id + 1] - reads_size_base;
    
printf("DEBUG:\n");
Indel_Accel_Krnl(consensus, con_size, &consensus_length[con_size_base], reads, rs_size, &reads_length[reads_size_base], qs, new_ref_idx, reads_size_base);
//Indel_Accel_Krnl(consensus_1, consensus_size_1, consensus_length_1, reads_1, reads_size_1, reads_length_1, qs_1, new_ref_idx_1);
//Indel_Accel_Krnl(consensus_2, consensus_size_2, consensus_length_2, reads_2, reads_size_2, reads_length_2, qs_2, new_ref_idx_2);
//Indel_Accel_Krnl(consensus_3, consensus_size_3, consensus_length_3, reads_3, reads_size_3, reads_length_3, qs_3, new_ref_idx_3);

printf("FINISH:\n");
}
}
