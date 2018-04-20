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
//whd(con_arr, con_size, con_len, reads_arr, reads_size, reads_len, weights_arr, min_whd, min_whd_idx);
void Indel_Accel_Krnl (ap_uint<4>* consensus, const int consensus_size, int* consensus_length, \
    ap_uint<4>* reads, const int reads_size, int* reads_length, char* qs, int* new_ref_idx, int new_ref_idx_base){
#pragma HLS INLINE
//#pragma HLS expression_balance
    ap_uint<4> con_buffer_0 [CON_SIZE * CON_LEN / 4];
    ap_uint<4> con_buffer_1 [CON_SIZE * CON_LEN / 4];
    ap_uint<4> con_buffer_2 [CON_SIZE * CON_LEN / 4];
    ap_uint<4> con_buffer_3 [CON_SIZE * CON_LEN / 4];
    //#pragma HLS array_partition variable=con_buffer block factor=8
    //#pragma HLS array_partition variable=con_buffer cyclic factor=4
    ap_uint<4> reads_buffer_0 [READS_SIZE * READS_LEN / 4];
    ap_uint<4> reads_buffer_1 [READS_SIZE * READS_LEN / 4];
    ap_uint<4> reads_buffer_2 [READS_SIZE * READS_LEN / 4];
    ap_uint<4> reads_buffer_3 [READS_SIZE * READS_LEN / 4];
    //#pragma HLS array_partition variable=reads_buffer block factor=8
    //#pragma HLS array_partition variable=reads_buffer cyclic factor=4
    char weights_buffer_0[READS_SIZE * READS_LEN / 4];
    char weights_buffer_1[READS_SIZE * READS_LEN / 4];
    char weights_buffer_2[READS_SIZE * READS_LEN / 4];
    char weights_buffer_3[READS_SIZE * READS_LEN / 4];
    //#pragma HLS array_partition variable=weights_buffer block factor=8
    //#pragma HLS array_partition variable=weights_buffer cyclic factor=8

    int consensus_length_buffer[CON_SIZE];
    #pragma HLS array_partition variable=consensus_length_buffer cyclic factor=2
    int reads_length_buffer[READS_SIZE];
    #pragma HLS array_partition variable=reads_length_buffer cyclic factor=2

    for (int i = 0; i < consensus_size + 1; i++) {
        consensus_length_buffer[i] = consensus_length[i];
    }

    for (int i = 0; i < reads_size + 1; i++) {
        reads_length_buffer[i] = reads_length[i];
    }

    int con_base = consensus_length_buffer[0];
    int con_end = consensus_length_buffer[consensus_size];
    
    int con_diff = con_end - con_base;
    int con_diff_bound = con_diff >> 2 + 1;
    for (int i = 0; i < con_diff_bound; i+=1) {
        //#pragma HLS unroll skip_exit_check factor=2
        int i4 = i << 2;
        con_buffer_0[i] = consensus[con_base + i4];
        con_buffer_1[i] = consensus[con_base + i4 + 1];
        con_buffer_2[i] = consensus[con_base + i4 + 2];
        con_buffer_3[i] = consensus[con_base + i4 + 3];
    }

    int rs_base = reads_length_buffer[0];
    int rs_end = reads_length_buffer[reads_size];
    int rs_diff = rs_end - rs_base; 
    int rs_diff_bound = rs_diff >> 2 + 1;
    for (int i = 0; i < rs_diff_bound; i++) {
        //#pragma HLS unroll skip_exit_check factor=2
        int i4 = i << 2;
        reads_buffer_0[i] = reads[rs_base + i4 + 0];
        reads_buffer_1[i] = reads[rs_base + i4 + 1];
        reads_buffer_2[i] = reads[rs_base + i4 + 2];
        reads_buffer_3[i] = reads[rs_base + i4 + 3];
        weights_buffer_0[i] = qs[rs_base + i4 + 0];
        weights_buffer_1[i] = qs[rs_base + i4 + 1];
        weights_buffer_2[i] = qs[rs_base + i4 + 2];
        weights_buffer_3[i] = qs[rs_base + i4 + 3];
    }

    //#pragma HLS DATAFLOW
    int min_whd[CON_SIZE * READS_SIZE];
    #pragma HLS array_partition variable=min_whd cyclic factor=4
    int min_whd_idx[CON_SIZE * READS_SIZE];
    #pragma HLS array_partition variable=min_whd_idx cyclic factor=4
    int new_ref[READS_SIZE];

    
    //int i, j, k, l;
    for (int i = 0; i < consensus_size; i++) {
        //int consensus_base = consensus_length[i];
        int consensus_base = consensus_length_buffer[i] - con_base;

        printf("con_base: %d\n", consensus_base);
        int local_consensus_length =  consensus_length_buffer[i+1] - consensus_length_buffer[i];
        for (int j = 0; j < reads_size; j++) {
        //#pragma HLS unroll factor=4
            //int reads_base = reads_length[j];
            int reads_base = reads_length_buffer[j] - rs_base;
            int local_reads_length = reads_length_buffer[j+1] - reads_length_buffer[j];
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
                //for (int l = 0; l < local_reads_length; l+=BLOCK_SIZE) {
//loop_conbuffer: for (int l = 0; l < local_reads_length; l++) {
              int bound = (local_reads_length >> 2) << 2; 
loop_conbuffer: for (int l = 0; l < bound; l+=4) {
                //#pragma HLS unroll factor=4
                //#pragma HLS unroll skip_exit_check factor=4
                              //printf("%c", consensus[consensus_base + k + l]);
                    //printf("%c", reads[reads_base + k + l]);
                    int con_start_addr =  consensus_base + k + l; 
                    int con_start_buf = con_start_addr & 0x3;
                    int con_start_addr_idx = con_start_addr >> 2;
                    int con_0, con_1, con_2, con_3;
                    con_0 = (con_start_buf > 0) ? con_start_addr_idx + 1 : 0;
                    con_1 = (con_start_buf > 1) ? con_start_addr_idx + 1 : 0;
                    con_2 = (con_start_buf > 2) ? con_start_addr_idx + 1 : 0;
                    con_3 = (con_start_buf > 3) ? con_start_addr_idx + 1 : 0;
                    
                    int r_start_addr =  reads_base + l; 
                    int r_start_buf = r_start_addr & 0x3;
                    int r_start_addr_idx = r_start_addr >> 2;
                    int r_0, r_1, r_2, r_3;
                    r_0 = (r_start_buf > 0) ? r_start_addr_idx + 1 : 0;
                    r_1 = (r_start_buf > 1) ? r_start_addr_idx + 1 : 0;
                    r_2 = (r_start_buf > 2) ? r_start_addr_idx + 1 : 0;
                    r_3 = (r_start_buf > 3) ? r_start_addr_idx + 1 : 0;
 

                    int whd_0 = (con_buffer_0[con_0] != reads_buffer_0[r_0]) ? weights_buffer_0[r_0] : 0;
                    int whd_1 = (con_buffer_1[con_1] != reads_buffer_1[r_1]) ? weights_buffer_1[r_1] : 0;
                    int whd_2 = (con_buffer_2[con_2] != reads_buffer_2[r_2]) ? weights_buffer_2[r_2] : 0;
                    int whd_3 = (con_buffer_3[con_3] != reads_buffer_3[r_3]) ? weights_buffer_3[r_3] : 0;
                    int whd_01 = whd_0 + whd_1;
                    int whd_23 = whd_2 + whd_3;
                    int whd_all = whd_01 + whd_23;
                    whd += whd_all;
                    //ap_uint<4> con_char = con_buffer[consensus_base + k + l];
                    //ap_uint<4> reads_char = reads_buffer[reads_base + l];
                    //if (con_buffer[consensus_base + k + l] != reads_buffer[reads_base + l]){
                    //if (con_char != reads_char) {
                    //    whd += weights_buffer[reads_base + l];
                        //if(k == 8 & j == 1){
                        //    printf("whd: %d\t", whd);
                        //}
                    //}                        

                    //printf("whd: %d\t", whd);
                    //printf("\t");

              
                   // ap_uint<4> reads_buffer[BLOCK_SIZE];
                   // #pragma HLS array_partition variable=reads_buffer cyclic factor=2
                   // ap_uint<4> con_buffer[BLOCK_SIZE];
                   // #pragma HLS array_partition variable=con_buffer cyclic factor=2
                   // char qs_buffer[BLOCK_SIZE]; 
                   // #pragma HLS array_partition variable=qs_buffer cyclic factor=2
                   // char whd_buffer[BLOCK_SIZE];
                   // #pragma HLS array_partition variable=whd_buffer cyclic factor=2
                   //     
                   // int block_size = (local_reads_length - l) > BLOCK_SIZE ? BLOCK_SIZE : (local_reads_length  - l);
                   // for (int ll = 0; ll < block_size; ll++) {
                   // //#pragma HLS unroll factor=2
                   // #pragma HLS loop_tripcount min=0 max=128
                   //     reads_buffer[ll] = reads[reads_base + l + ll]; 
                   //     qs_buffer[ll] = qs[reads_base + l + ll]; 
                   //     con_buffer[ll] = consensus[consensus_base + k + l + ll];
                   // }

                   //
                   // for (int ll = 0; ll < block_size; ll++) {
                   // #pragma HLS loop_tripcount min=0 max=128
                   // #pragma HLS unroll factor=2
                   // //#pragma HLS unroll factor=8
                   //     char con_char = con_buffer[ll];
                   //     con_char = con_char & 0xf;
                   //     char reads_char = reads_buffer[ll];
                   //     reads_char = reads_char & 0xf;
                   //     printf("reads_char: %x -- con_char: %x; ", reads_char, con_char);
                   //     printf("qs %d %d con %d %d,", reads_base, l + ll,  consensus_base, k + l + ll);
                   //     whd_buffer[ll] = (con_buffer[ll] != reads_buffer[ll]) ? qs_buffer[ll] : 0;
                   // }

                   // for (int ll = 0; ll < block_size; ll++) {
                   // #pragma HLS loop_tripcount min=0 max=128
                   // #pragma HLS unroll factor=2
                   // #pragma HLS expression_balance
                   //     whd += whd_buffer[ll];   
                   // }
                }
                //int whd = whd_ptr[0];
                //printf("whd: %d\t", whd);
                //printf("\t");


                //printf("\n");
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
   
    printf("Finish WHD\n");
    
    int min_score = 0x7fffffff;
    int min_idx = consensus_size + 1;
    for (int i = 1; i < consensus_size; i++) {
        int score = 0;
        for (int j = 0; j < reads_size; j++) {
        #pragma HLS unroll factor=2
            int tmp = min_whd[i * reads_size + j] - min_whd[j];
            score += (tmp > 0) ? tmp: -tmp;
        }
        min_idx = (score < min_score) ? i : min_idx;
        //scores[i] = score;
    }
    //printf( "min_idx: %d\n", min_idx);
    //assert(min_idx < consensus_size);
    
    for (int j = 0; j < reads_size; j++) {
    #pragma HLS unroll factor=2
        //if ( min_whd[ min_idx * reads_size + j] < min_whd[j]){
            //new_ref[j] = min_whd[min_idx][j];
            //new_ref_idx[j] = min_whd_idx[min_idx][j];
            //new_ref[new_ref_idx_base + j] = min_whd[min_idx * reads_size + j];
            new_ref[j] = min_whd[min_idx * reads_size + j];
            new_ref_idx[new_ref_idx_base + j] = min_whd_idx[ min_idx * reads_size +j];
 
        //} else{
        //    new_ref[j] = min_whd [j];
        //    new_ref_idx[j] = min_whd_idx[j]; 
        //}
    }
    for (int j = 0; j < reads_size; j++) {
     printf("Read %2d whd %2d index %2d\n", j, new_ref[j], new_ref_idx[new_ref_idx_base + j]);
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
for (int itr = 0; itr < global_threads; itr ++){
#pragma HLS unroll factor=4
int index = (global_id << 2) + itr;
int con_size_base = consensus_size[index];
int con_size = consensus_size[index + 1] - con_size_base;
printf("itr: %d index: %d, global_id %d \n", itr, index, global_id);
printf("con_size_base: %d \t con_size: %d \n", con_size_base, con_size);

int reads_size_base = reads_size[index];
int rs_size = reads_size[index + 1] - reads_size_base;
printf("reads_size_base: %d \t reads_size: %d \n", reads_size_base, rs_size);
    
printf("DEBUG:\n");
printf("\n");
Indel_Accel_Krnl(consensus, con_size, &consensus_length[con_size_base], reads, rs_size, &reads_length[reads_size_base], qs, new_ref_idx, reads_size_base);
//Indel_Accel_Krnl(consensus_1, consensus_size_1, consensus_length_1, reads_1, reads_size_1, reads_length_1, qs_1, new_ref_idx_1);
//Indel_Accel_Krnl(consensus_2, consensus_size_2, consensus_length_2, reads_2, reads_size_2, reads_length_2, qs_2, new_ref_idx_2);
//Indel_Accel_Krnl(consensus_3, consensus_size_3, consensus_length_3, reads_3, reads_size_3, reads_length_3, qs_3, new_ref_idx_3);

printf("FINISH:\n");
}
}
}
