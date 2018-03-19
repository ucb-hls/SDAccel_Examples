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
    score: for (int i = 1; i < consensus_size; i++) {
        int score = 0;
        for (int j = 0; j < reads_size; j++) {
            int tmp = min_whd[(i * reads_size + j) << 1] - min_whd[j << 1];
            score += (tmp > 0) ? tmp: -tmp;
        }
        min_idx = (score < min_score) ? i : min_idx;
        //scores[i] = score;
    }
    //printf( "min_idx: %d\n", min_idx);
    //assert(min_idx < consensus_size);
    
    rank: for (int j = 0; j < reads_size; j++) {
            new_ref[j] = min_whd[(min_idx * reads_size + j) << 1];
            new_ref_idx[j] = min_whd[((min_idx * reads_size +j) << 1) + 1];
            }

    print: for (int j = 0; j < reads_size; j++) {
     //printf("Read %2d whd %2d index %2d\n", j, new_ref[j], new_ref_idx[j]);
        printf("Kernel: Read %2d whd %4d  index %2d\n", j, new_ref[j], new_ref_idx[j]);
    }
}

extern "C" {
//whd(con_arr, con_size, con_len, reads_arr, reads_size, reads_len, weights_arr, min_whd, min_whd_idx);
void Indel_Accel_Krnl (ap_uint<4>* consensus, const int consensus_size, int* consensus_length, \
    ap_uint<4>* reads, const int reads_size, int* reads_length, char* qs, int* new_ref_idx){
#pragma HLS INLINE
    
         
    //Set buffer to 512 bit -> 128 char 
    // Tile the output with the consensus 
    ap_uint<4> con_buffer_0[128];
    ap_uint<4> con_buffer_1[128];
    char qs_buffer_0[128];
    char qs_buffer_1[128];
    ap_uint<4> reads_buffer_0[128];
    ap_uint<4> reads_buffer_1[128];


    // Buffer the reads_length 
    int consensus_length_buffer[CON_SIZE];
    int reads_length_buffer[READS_SIZE];

    int consensus_offset_buffer[CON_SIZE];
    int reads_offset_buffer[READS_SIZE];


    int count = 0; 
    for (int i = 0; i < consensus_size; i++) {
        consensus_offset_buffer[i] = count;               
        int local_consensus_length = consensus_length[i];
        consensus_length_buffer[i] = local_consensus_length;
        count += local_consensus_length;
    }
   
    count = 0;
    for (int i = 0; i < reads_size; i++) {
        reads_offset_buffer[i] = count;               
        int local_reads_length = reads_length[i];
        reads_length_buffer[i] = local_reads_length;
        count += local_reads_length;
    }
 
    //array of stream declaration
    /*hls::stream<ap_uint<512>> condStreams[CON_SIZE];
    //#pragma HLS STREAM variable=condStreams depth=32

    hls::stream<ap_uint<512>> readStreams[READS_SIZE];
    //#pragma HLS STREAM variable=readStreams depth=32


    #pragma HLS dataflow
    for (int i = 0 ; i < CON_SIZE; i ++){
        read_input(consensus, condStreams[i], consensus_length[i]);
    }
   
    #pragma HLS dataflow
    for (int i = 0 ; i < READS_SIZE; i ++){
        read_input(reads, readStreams[i], reads_length[i]);
    }*/
    int min_whd[CON_SIZE * READS_SIZE << 1];
    //int min_whd_idx[CON_SIZE * READS_SIZE];
    int new_ref[READS_SIZE];
    //int consensus_base = 0; 

    // 512 / 4 -> 128 vectors 
    //int consensus_size_in512 = (consensus_size-1)/128 + 1;
    //int reads_size_in512 = ((reads_size-1)/128 + 1) * 128;
 
    //int trunck_size = (reads_size % global_threads == 0) ? (reads_size / global_threads) : ((reads_size  + global_threads) / global_threads);
    //int thread_work_start = global_id * trunck_size;
    //int thread_work_end = (global_id + 1) * trunck_size;
    //thread_work_end = (thread_work_end < reads_size) ? thread_work_end: reads_size;
   
    //reads_size: for (j = thread_work_start; j < thread_work_end; j++ ) {
    reads_size: for (int j = 0; j < reads_size; j++ ) {
    #pragma HLS unroll factor=8
    #pragma HLS LOOP_TRIPCOUNT min=1 max=256

        int reads_base = reads_offset_buffer[j];
        int local_reads_length = reads_length_buffer[j];

        consensus_size: for (int i = 0; i < consensus_size; i++) {
        #pragma HLS LOOP_TRIPCOUNT min=1 max=32

            int consensus_base = consensus_offset_buffer[i]; 
            int local_consensus_length =  consensus_length_buffer[i];
     

            //printf( "consensus size %d i %d consensus length %d, read size %d j %d reads length %d\n", \
                consensus_size, i, local_consensus_length, reads_size,  j, local_reads_length);
            //int reads_base = 0
            // Prefetech cons 


            int min = 0x7fffffff; 
            int min_idx = local_consensus_length - local_reads_length + 1;

            // Assume not tile accross the whd 
            //for (int kk = 0; kk < 128; k++){
            //    con_buffer_0[kk] = consensus[consensus_base + kk];
            //} 
            // kth element location in the con_buffer
            int con_start = 0;

            printf("Local con_len: %d, reads_len: %d\n", local_consensus_length, local_reads_length);
            printf("Length Diff: %d\n", local_consensus_length - local_reads_length);
           //con_reads_diff: for (k = 0; k <= local_consensus_length - local_reads_length; k++) {
            con_reads_diff: for (int k = 0; k <= local_consensus_length - local_reads_length; k+=1) {
            #pragma HLS unroll factor=2
            #pragma HLS LOOP_TRIPCOUNT min=1 max=1792

                printf("k=%d\n",k);
                // whd 
                int whd = 0;

                int vec_begin, vec_end;
                int vec_begin_offset, vec_end_offset; 
                vec_begin = reads_base >> 7;
                vec_begin_offset = reads_base % 128;
                vec_end = ((local_reads_length + reads_base) >> 7) + 1;
                vec_end_offset = (local_reads_length + reads_base) % 128;
                /*reads_buffer[0] = reads[vec_begin]; 
                if (vec_begin_offset != 0) {
                    for (int v = vec_begin_offset; v < 128; v++){
                        if (consensus[consensus_base + k + v - vec_begin_offset] !=reads_buffer[v]){
                            whd += qs[reads_base + v];
                        }
                    }
                    vec_begin++;
                }*/
                
               int abs_reads_base= vec_begin << 7;
               int rlt_con_base = 0 - vec_begin_offset;
                vec_reads: for (int ll_r = 0; ll_r < 128; ll_r++){
                    if (vec_begin % 2 == 0){
                        reads_buffer_0[ll_r] = reads[abs_reads_base + ll_r];
                        qs_buffer_0[ll_r] =qs[abs_reads_base + ll_r];
                        //if (rlt_con_base + ll_r >= 0)
                            con_buffer_0[ll_r] = consensus[consensus_base + k + rlt_con_base + ll_r];
                    } else {
                        reads_buffer_1[ll_r] = reads[abs_reads_base + ll_r];
                        qs_buffer_1[ll_r] =qs[abs_reads_base + ll_r];
                        //if (rlt_con_base + ll_r >= 0)
                            con_buffer_1[ll_r] = consensus[consensus_base + k + rlt_con_base + ll_r];
                    }
                }

                int rlt_index = 0;
                vec_run: for(int ll = vec_begin; ll < vec_end; ll+=1){

                    abs_reads_base= ll << 7;
                    rlt_con_base = ((ll-vec_begin) << 7) - vec_begin_offset;

                    // Fill the read buffer, hope to enable burst
                    for (int ll_r = 0; ll_r < 128; ll_r++){
                        if (ll % 2 == 0){
                            reads_buffer_1[ll_r] = reads[abs_reads_base + ll_r + 128];
                            qs_buffer_1[ll_r] = qs[abs_reads_base + ll_r + 128];
                            con_buffer_1[ll_r] = consensus[consensus_base + k + rlt_con_base + ll_r + 128];
                        } else {
                            reads_buffer_0[ll_r] = reads[abs_reads_base + ll_r + 128];
                            qs_buffer_0[ll_r] =qs[abs_reads_base + ll_r + 128];
                            con_buffer_0[ll_r] = consensus[consensus_base + k + rlt_con_base + ll_r + 128];
                        }
                    }

                    int chunk_begin = (ll == vec_begin)? vec_begin_offset : 0;
                    int chunk_end = (ll == vec_end - 1)? vec_end_offset: 128;
                    

                    chunk_run: for (int v = chunk_begin; v < chunk_end; v++){
                    #pragma HLS UNROLL
                        //long long  print_var = reads_buffer[0].range(63, 0);
                        //printf("buffer: %x", reads_buffer[0].range(31,0));
                        //printf("buffer: %lx ", print_var);
                        //char reads_char = (reads_buffer[0] >> (v << 3)) & 0xff; 
                        char reads_char = (ll % 2 == 0) ? reads_buffer_0[v]: reads_buffer_1[v];
                        reads_char = reads_char & 0xf;
                        //char con_char =  consensus[consensus_base + k + rlt_index]; 
                        int con_idx = (k + rlt_index) % 128;
                        char con_char =  (ll % 2 == 0) ? con_buffer_0[v]: con_buffer_1[v]; 
                        con_char = con_char & 0xf;
                        printf("reads_char: %x -- con_char: %x; ", reads_char, con_char);
                        printf("qs %d %d con %d %d,", abs_reads_base, abs_reads_base+ v,  consensus_base, consensus_base + k + rlt_index);
                        if (con_char != reads_char){
                
                            //whd += qs[abs_reads_base + v];
                            whd += (ll % 2 == 0) ? qs_buffer_0[v] : qs_buffer_1[v];
                        }
                        rlt_index++;
                    }
    
                    printf("qs chunk end\n");

                }
                printf("whd %d\n\n", whd);
                
                // Optimization tree based reduction
                // ll starts from reads_base to reads_base + local_reads_length 
                /*for(int ll = 0; ll < local_reads_length; ll+=128){
                    int chunk_size = 128;
                    int remain_size = local_reads_length - ll; 
                    chunk_size = (size < chunk_size) ? size: chunk_size;    

                    // If this is the beginning of the 128 block
                    int reads_address = (ll + reads_base) >> 7;
                    int reads_offset = (ll + reads_base) % 128;
                    reads_buffer[0] = reads[reads_address];

                    // If this is not a multiple of 
                    reads_length: for (l = 0; l < chunk_size; l++) {
                    //reads_length: for (l = 0; l < local_reads_length; l++) {
                    #pragma HLS LOOP_TRIPCOUNT min=1 max=256
                    #pragma HLS UNROLL 


                    //printf("%c", consensus[consensus_base + k + l]);
                    //printf("%c", reads[reads_base + k + l]);
                    //if (consensus[consensus_base + k + l] !=reads[reads_base + l]){
                    int chunk_offset = (reads_base + l) % chunk_size;
                    if (consensus[consensus_base + k + l] !=reads_buffer[chunk_offset]){
                        whd += qs[reads_base + l];
                        //if(k == 8 & j == 1){
                        //    printf("whd: %d\t", whd);
                        //}
                    }                        
                    //printf("whd: %d\t", whd);
                    //printf("\t");

                }*/

                //printf("\n");
                if (whd < min) {
                    min =  whd; 
                    min_idx = k; 
               }

            }
            printf( "[%d, %d]-Kernel:min_idx %d min %d\n", i, j, min_idx, min);
            //assert(min_idx <= local_consensus_length - local_reads_length);
            
            //min_whd[i][j] = min;
            //min_whd_idx[i][j] = min_idx;
            //min_whd[i * reads_size + j] = min;
            //min_whd_idx[i * reads_size + j] = min_idx;

            min_whd[(i * reads_size + j) << 1] = min;
            min_whd[((i * reads_size + j) << 1) + 1] = min_idx;
 
 
            //reads_base += local_reads_length;
            //reads_base = reads_offset_buffer[j];
    
        }
        //consensus_base += local_consensus_length;
        //printf("%d %d\n", consensus_base, consensus_offset_buffer[i]);
        //consensus_base = consensus_offset_buffer[i];
    }
    
    int min_score = 0x7fffffff;
    int min_idx = consensus_size + 1;
    score: for (int i = 1; i < consensus_size; i++) {
        int score = 0;
        for (int j = 0; j < reads_size; j++) {
            int tmp = min_whd[(i * reads_size + j) << 1] - min_whd[j << 1];
            score += (tmp > 0) ? tmp: -tmp;
        }
        min_idx = (score < min_score) ? i : min_idx;
        //scores[i] = score;
    }
    //printf( "min_idx: %d\n", min_idx);
    //assert(min_idx < consensus_size);
    
    rank: for (int j = 0; j < reads_size; j++) {
        //if ( min_whd[ min_idx * reads_size + j] < min_whd[j]){
            //new_ref[j] = min_whd[min_idx][j];
            //new_ref_idx[j] = min_whd_idx[min_idx][j];
            new_ref[j] = min_whd[(min_idx * reads_size + j) << 1];
            new_ref_idx[j] = min_whd[ ((min_idx * reads_size +j) << 1) + 1 ];
 
        //} else{
        //    new_ref[j] = min_whd [j];
        //    new_ref_idx[j] = min_whd_idx[j]; 
        //}
    }
    print: for (int j = 0; j < reads_size; j++) {
     //printf("Read %2d whd %2d index %2d\n", j, new_ref[j], new_ref_idx[j]);
        printf("Kernel: Read %2d whd %4d  index %2d\n", j, new_ref[j], new_ref_idx[j]);
    }



}


void Indel_Accel (ap_uint<4>* consensus_0, const int consensus_size_0, int* consensus_length_0, \
    ap_uint<4>* reads_0, const int reads_size_0, int* reads_length_0, char* qs_0, int* new_ref_idx_0, \
    ap_uint<4>* consensus_1, const int consensus_size_1, int* consensus_length_1, \
    ap_uint<4>* reads_1, const int reads_size_1, int* reads_length_1, char* qs_1, int* new_ref_idx_1, \
    ap_uint<4>* consensus_2, const int consensus_size_2, int* consensus_length_2, \
    ap_uint<4>* reads_2, const int reads_size_2, int* reads_length_2, char* qs_2, int* new_ref_idx_2, \ 
    ap_uint<4>* consensus_3, const int consensus_size_3, int* consensus_length_3, \
    ap_uint<4>* reads_3, const int reads_size_3, int* reads_length_3, char* qs_3, int* new_ref_idx_3, int global_id, int global_threads) {
 
 
    //ap_uint<4>* reads, const int reads_size, int* reads_length, char* qs, int* new_ref_idx, int* new_ref_idx, int global_id, int global_threads) {
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
 

#pragma HLS expression_balance
Indel_Accel_Krnl(consensus_0, consensus_size_0, consensus_length_0, reads_0, reads_size_0, reads_length_0, qs_0, new_ref_idx_0);
Indel_Accel_Krnl(consensus_1, consensus_size_1, consensus_length_1, reads_1, reads_size_1, reads_length_1, qs_1, new_ref_idx_1);
Indel_Accel_Krnl(consensus_2, consensus_size_2, consensus_length_2, reads_2, reads_size_2, reads_length_2, qs_2, new_ref_idx_2);
Indel_Accel_Krnl(consensus_3, consensus_size_3, consensus_length_3, reads_3, reads_size_3, reads_length_3, qs_3, new_ref_idx_3);
//ap_uint<4>* consensus;
//const int consensus_size;
//int* consensus_length;
//ap_uint<4>* reads;
//const int reads_size;
//int* reads_length;
//char* qs;
//int* new_ref_idx;


}
}
