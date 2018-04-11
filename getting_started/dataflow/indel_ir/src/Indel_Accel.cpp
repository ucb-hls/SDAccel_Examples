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

typedef ap_uint<8> TYPE;
typedef ap_uint<1024> CON_SIZE_TYPE;
#define nConsensus 32

//extern "C" {
void Indel_Accel_Krnl (
                unsigned int consAddr, unsigned int readAddr, unsigned int qualAddr, \
                unsigned int consLen, unsigned int readLen, \
                unsigned int posAddr, unsigned int swapAddr, \
                unsigned int targetPos, CON_SIZE_TYPE conSize, \
                TYPE* cons, TYPE* reads, TYPE* qual, \
                unsigned int* pos, unsigned char* swap){
#pragma HLS INLINE
#pragma HLS expression_balance
        int min_whd[CON_SIZE][READS_SIZE];
        int min_whd_idx[CON_SIZE][READS_SIZE];
        int new_ref[READS_SIZE];
        int new_ref_idx[READS_SIZE];
        int scores[CON_SIZE];

        for (size_t i = 0; i < consLen; i++) {
                //int local_consensus_length = conSize[i];
                int local_consensus_length = conSize.range( i << 5 + 31,i << 5);
                for (size_t j = 0; j < readLen; j++) {
                //#pragma HLS pipeline II=1
                //#pragma HLS unroll factor=8 region
                        //int local_reads_length = reads_length[j];
                        int local_reads_length = 0;

                        //printf( "consensus size %d i %d consensus length %d, read size %d j %d reads length %d\n", \
                        consLen, i, local_consensus_length, readLen,  j, local_reads_length);
                        int min = 0x7fffffff; 
                        int min_idx = 0x7fffffff;;

                        int whd_0 = 0;
                        // get local_reads_length based on quality
                        for (size_t l = 0; l < READS_LEN; l++) {
                        
                                unsigned char qs = qual[readAddr + l];
                                if (qs == 0){ 
                                        local_reads_length = l; 
                                        break;
                                }

                                if (cons[consAddr + 0 + l] != reads[readAddr + l]){
                                        whd_0 += qs;
                                        //if(k == 8 & j == 1){
                                        //    printf("whd: %d\t", whd);
                                        //}
                                }                        
                                //printf("whd: %d\t", whd);
                                //printf("\t");

                        }
                        //printf("\n");
                        if (whd_0 < min) {
                                min =  whd_0; 
                                min_idx = 0; 
                        }



                        for (size_t k = 1; k <= local_consensus_length - local_reads_length; k++) {

                                // whd 
                                int whd = 0;
                                // Optimization tree based reduction
                                //for (int l = 0; l < reads_length[j]; l++) {
                                for (size_t l = 0; l < local_reads_length; l++) {
                                        //printf("%c", consensus[consAddr + k + l]);
                                        //printf("%c", reads[readAddr + k + l]);
                                        if (cons[consAddr + k + l] != reads[readAddr + l]){
                                                whd += qual[readAddr + l];
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

                        min_whd[i][j] = min;
                        min_whd_idx[i][j] = min_idx;
                        //min_whd[i * readLen + j] = min;
                        //min_whd_idx[i * readLen + j] = min_idx;
                        readAddr += local_reads_length;
                }

                consAddr += local_consensus_length;

        }

        int min_score = 0x7fffffff;
        int min_idx = consLen + 1;
        for (size_t i = 1; i < consLen; i++) {
                int score = 0;
                for (size_t j = 0; j < readLen; j++) {
                        //int tmp = min_whd[i * readLen + j] - min_whd[j];
                        int tmp = min_whd[i][j] - min_whd[0][j];
                        score += (tmp > 0) ? tmp: -tmp;
                }
                min_idx = (score < min_score) ? i : min_idx;
                //scores[i] = score;
        }
        //printf( "min_idx: %d\n", min_idx);
        //assert(min_idx < consLen);

        for (size_t j = 0; j < readLen; j++) {
                //if ( min_whd[ min_idx * readLen + j] < min_whd[j]){
                if ( min_whd [min_idx][j] < min_whd[0][j]){
                        swap[swapAddr + j] = 1;
                        //pos[posAddr + j] = min_whd_idx[ min_idx * readLen +j] + targetPos;
                        pos[posAddr + j] = min_whd_idx[min_idx][j] + targetPos;
                        //new_ref[j] = min_whd[min_idx * readLen + j];
                        new_ref[j] = min_whd[min_idx][j];
                        //new_ref_idx[j] = min_whd_idx[ min_idx * readLen +j];
                        new_ref_idx[j] = min_whd_idx[min_idx][j];
                } else {
                        swap[swapAddr + j] = 0;
                }
        }
        for (size_t j = 0; j < readLen; j++) {
                printf("Read %2zu whd %2d index %2d\n", j, new_ref[j], new_ref_idx[j]);
        }
}

// Len is the size, size is the len
void Indel_Accel(unsigned int id, unsigned int rd, \
                unsigned int consAddr, unsigned int readAddr, unsigned int qualAddr, \
                unsigned int consLen, unsigned int readLen, \
                unsigned int posAddr, unsigned int swapAddr, \
                unsigned int targetPos, CON_SIZE_TYPE conSize, \
                TYPE* cons, TYPE* reads, TYPE* qual, \
                unsigned int* pos, unsigned char* swap, \
                unsigned int ret_id, unsigned int ret_rd) {
#pragma HLS INTERFACE m_axi depth=64 port=cons offset=off bundle=cons
#pragma HLS INTERFACE m_axi depth=64 port=reads offset=off bundle=reads
#pragma HLS INTERFACE m_axi depth=64 port=qual offset=off bundle=qual
#pragma HLS INTERFACE m_axi depth=64 port=pos offset=off bundle=pos
#pragma HLS INTERFACE m_axi depth=64 port=swap offset=off bundle=swap
//#pragma HLS INTERFACE m_axi depth=64 port=conSize offset=off bundle=gmem

        Indel_Accel_Krnl(consAddr, readAddr, qualAddr, consLen, readLen, posAddr, swapAddr, targetPos, conSize, cons, reads, qual, pos, swap);


        ret_id = id;
        ret_rd = rd;
        return;
}
