#include <hls_stream.h>
#include <stdio.h>
#include <string.h>
#include <ap_int.h>

#include "Indel_Accel.h"
//#define CON_SIZE 2
//#define READS_SIZE 8
#define BLOCK_SIZE 128

// JENNY TODO
// ap_uint<3> for ATGCU 
// ap_uint<512> to round it to power of 2 
// ap_int<512> is the most efficient for global memory 
// Assign reads and consensus to different bundles for memory banking optimizations

extern "C" {

void Indel_Rank(int * min_whd, int * min_whd_idx, /*int* con_size_base, int* con_size_base_1, int* reads_size_base, int* reads_size_base_1, */int* new_ref_idx, int global_id, int con_size_local_base, int consensus_size, int reads_size_local_base, int reads_size){

    //int con_size_local_base = con_size_base[global_id];
    //int consensus_size = con_size_base_1[global_id] - con_size_local_base;
    //int reads_size_local_base = reads_size_base[global_id];
    //int reads_size = reads_size_base_1[global_id] - reads_size_local_base;

    int min_score = 0x7fffffff;
    int min_idx = -0x7fffffff;
    int iterations = (consensus_size - 1) * reads_size; 
    int score = 0;
   // for (int itr = 0, i = 1, j = 0; itr < iterations; itr ++){
   // #pragma HLS PIPELINE
   // #pragma HLS unroll factor=4

   //     int tmp = min_whd[i * reads_size + j] - min_whd[j];
   //     score += (tmp > 0) ? tmp: -tmp;
   //     if (j = reads_size - 1){
   //         min_idx = (score < min_score) ? i : min_idx;
   //         j = 0; 
   //         score = 0;
   //         i++;
   //     } else{
   //         j++; 
   //     }
   // }

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
    printf( "min_idx: %d\n", min_idx);
    //assert(min_idx < consensus_size);
    
    for (int j = 0; j < reads_size; j++) {
    #pragma HLS unroll factor=4
        //if ( min_whd[ min_idx * reads_size + j] < min_whd[j]){
            //new_ref[j] = min_whd[min_idx][j];
            //new_ref_idx[j] = min_whd_idx[min_idx][j];
            //new_ref[new_ref_idx_base + j] = min_whd[min_idx * reads_size + j];
            new_ref_idx[reads_size_local_base + j] = min_whd_idx[ min_idx * reads_size +j];
 
        //} else{
        //    new_ref[j] = min_whd [j];
        //    new_ref_idx[j] = min_whd_idx[j]; 
        //}
    }
    //for (int j = 0; j < reads_size; j++) {
    // printf("Read %2d whd %2d index %2d\n", j, new_ref[new_ref_idx_base + j], new_ref_idx[new_ref_idx_base + j]);
    //}
}
//whd(con_arr, con_size, con_len, reads_arr, reads_size, reads_len, weights_arr, min_whd, min_whd_idx);

void Indel_Accel_Krnl (ap_uint<4>* consensus, int* con_size_base, int* con_size_base_1, int* abs_consensus_length, \
    ap_uint<4>* reads, int* reads_size_base, int* reads_size_base_1, int* abs_reads_length, char* qs, int* new_ref_idx, int global_id, int global_threads, int* min_whd, int* min_whd_idx, int* con_size_local_base_ptr, int* consensus_size_ptr, int* reads_size_local_base_ptr, int* reads_size_ptr){

int con_size_local_base = con_size_base[global_id];
*con_size_local_base_ptr = con_size_local_base;
int consensus_size = con_size_base_1[global_id] - con_size_local_base;
*consensus_size_ptr = consensus_size;
int reads_size_local_base = reads_size_base[global_id];
*reads_size_local_base_ptr = reads_size_local_base;
int reads_size = reads_size_base_1[global_id] - reads_size_local_base;
*reads_size_ptr = reads_size;
 
 
int * consensus_length = &abs_consensus_length[con_size_local_base];
int * reads_length = &abs_reads_length[reads_size_local_base];
//void Indel_Accel_Krnl (ap_uint<4>* consensus, const int consensus_size, int* consensus_length, \
//    ap_uint<4>* reads, const int reads_size, int* reads_length, char* qs, int* new_ref_idx, int* min_whd, int* min_whd_idx, int new_ref_idx_base){

//#pragma HLS INLINE
    //#pragma HLS DATAFLOW
    #pragma HLS expression_balance
   //int i, j, k, l;
    for (int i = 0; i < consensus_size; i++) {
        int consensus_base = consensus_length[i];

        printf("con_base: %d\n", consensus_base);
        int local_consensus_length =  consensus_length[i+1] - consensus_length[i];
        for (int j = 0; j < reads_size; j++) {
        #pragma HLS unroll factor=16
            int reads_base = reads_length[j];
            int local_reads_length = reads_length[j+1] - reads_length[j];
            printf( "consensus size %d i %d consensus length %d, read size %d j %d reads length %d\n", \
                consensus_size, i, local_consensus_length, reads_size,  j, local_reads_length);
            int min = 0x7fffffff; 
            int min_idx = 0x7fffffff;

            printf("Length Diff: %d\n", local_consensus_length - local_reads_length);
            for (int k = 0; k <= local_consensus_length - local_reads_length; k++) {
            #pragma HLS loop_tripcount min=0 max=2048

            int iterations = local_reads_length * (local_consensus_length - local_reads_length) ;
               // whd 
                int whd = 0;
                // Optimization tree based reduction
                for (int l = 0; l < local_reads_length; l+=BLOCK_SIZE) {
                //#pragma HLS unroll factor=4
                        
                    ap_uint<4> reads_buffer[BLOCK_SIZE];
                    #pragma HLS array_partition variable=reads_buffer cyclic factor=4
                    ap_uint<4> con_buffer[BLOCK_SIZE];
                    #pragma HLS array_partition variable=con_buffer cyclic factor=4
                    char qs_buffer[BLOCK_SIZE]; 
                    #pragma HLS array_partition variable=qs_buffer cyclic factor=4
                    char whd_buffer[BLOCK_SIZE];
                    #pragma HLS array_partition variable=whd_buffer cyclic factor=4
                        
                    int block_size = (local_reads_length - l) > BLOCK_SIZE ? BLOCK_SIZE : (local_reads_length  - l);
                    for (int ll = 0; ll < block_size; ll++) {
                    //#pragma HLS DATAFLOW
                    //#pragma HLS unroll factor=4
                    #pragma HLS loop_tripcount min=0 max=128
                        reads_buffer[ll] = reads[reads_base + l + ll]; 
                        qs_buffer[ll] = qs[reads_base + l + ll]; 
                        con_buffer[ll] = consensus[consensus_base + k + l + ll];
                    }

                   
                    for (int ll = 0; ll < block_size; ll++) {
                    #pragma HLS loop_tripcount min=0 max=128
                    #pragma HLS unroll factor=4
                    //yy#pragma HLS unroll factor=8
                        char con_char = con_buffer[ll];
                        con_char = con_char & 0xf;
                        char reads_char = reads_buffer[ll];
                        reads_char = reads_char & 0xf;
                        printf("reads_char: %x -- con_char: %x; ", reads_char, con_char);
                        printf("qs %d %d con %d %d,", reads_base, l + ll,  consensus_base, k + l + ll);
                        whd_buffer[ll] = (con_buffer[ll] != reads_buffer[ll]) ? qs_buffer[ll] : 0;
                    }

                    for (int ll = 0; ll < block_size; ll++) {
                    #pragma HLS loop_tripcount min=0 max=128
                    #pragma HLS unroll factor=4
                    #pragma HLS expression_balance
                        whd += whd_buffer[ll];   
                    }
                }
                //int whd = whd_ptr[0];
                printf("whd: %d\t", whd);
                printf("\t");



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
 

#pragma HLS expression_balance
//Copy_Con(consensus, consensus_size, consensus_length);
// 0-2 2-4
// Number of lengths 
//int con_size_base = consensus_size[global_id];
//int con_size = consensus_size[global_id + 1] - con_size_base;
//printf("con_size_base[%d]: %d\n", global_id, con_size);
//
//int reads_size_base = reads_size[global_id];
//int rs_size = reads_size[global_id + 1] - reads_size_base;
//printf("reads_size_base: %d\n", rs_size);
//
int min_whd[CON_SIZE * READS_SIZE];
//#pragma HLS array_partition variable=min_whd cyclic factor=16
int min_whd_idx[CON_SIZE * READS_SIZE];
//#pragma HLS array_partition variable=min_whd_idx cyclic factor=16

int con_size_local_base, consensus_size_ptr, reads_size_local_base, reads_size_ptr;

printf("DEBUG:\n");
#pragma HLS DATAFLOW

Indel_Accel_Krnl(consensus, consensus_size, &consensus_size[1], consensus_length, reads, reads_size, &reads_size[1], reads_length, qs, new_ref_idx, global_id, global_threads, min_whd, min_whd_idx, & con_size_local_base, & consensus_size_ptr, & reads_size_local_base, &reads_size_ptr);
//Indel_Accel_Krnl(consensus, con_size, &consensus_length[con_size_base], reads, rs_size, &reads_length[reads_size_base], qs, new_ref_idx, min_whd, min_whd_idx, reads_size_base);
Indel_Rank( min_whd,  min_whd_idx, /*consensus_size, &consensus_size[1], reads_size, &reads_size[1],*/  new_ref_idx, global_id, con_size_local_base, consensus_size_ptr, reads_size_local_base, reads_size_ptr);

//Indel_Accel_Krnl(consensus_1, consensus_size_1, consensus_length_1, reads_1, reads_size_1, reads_length_1, qs_1, new_ref_idx_1);
//Indel_Accel_Krnl(consensus_2, consensus_size_2, consensus_length_2, reads_2, reads_size_2, reads_length_2, qs_2, new_ref_idx_2);
//Indel_Accel_Krnl(consensus_3, consensus_size_3, consensus_length_3, reads_3, reads_size_3, reads_length_3, qs_3, new_ref_idx_3);

printf("FINISH:\n");
}
}
