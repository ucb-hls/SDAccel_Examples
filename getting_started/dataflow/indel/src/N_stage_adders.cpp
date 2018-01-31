/**********
Copyright (c) 2017, Xilinx, Inc.
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice,
this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors
may be used to endorse or promote products derived from this software
without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
**********/

/*******************************************************************************
Description: 
    HLS Dataflow Example using array of HLS Stream datatype
    This example is N Stages of Vector Addition to demonstrate Array of 
    Stream usage in HLS C Kernel. 

    1) read_input(): 
        This API reads the input vector from Global Memory and writes it into 
        HLS Stream inStream using blocking write command.

    2) Multiple instance of Adder(): 
        This API reads the input vector from inStream using blocking
        read command and increment the value by user specified increment. It writes 
        the results into outStream using blocking write command.

    3) write_result(): 
        This API reads the result vector from outStream using blocking read
        command and write the result into Global Memory Location.
        
    Four Stage Adder will be implemented as below:

                     _____________
                    |             |<----- Input Vector from Global Memory
                    |  read_input |       __
                    |_____________|----->|  |
                     _____________       |  | streamArray[0]
                    |             |<-----|__|
                    |   adder_0   |       __
                    |_____________|----->|  |
                     _____________       |  | streamArray[1]
                    |             |<-----|__|
                    |   adder_1   |       __
                    |_____________|----->|  |
                     _____________       |  | streamArray[2]
                    |             |<-----|__|
                    |   adder_2   |       __
                    |_____________|----->|  |
                     _____________       |  | streamArray[3]
                    |             |<-----|__|
                    |   adder_3   |       __
                    |_____________|----->|  |
                     ______________      |  | streamArray[4]
                    |              |<----|__|
                    | write_result |       
                    |______________|-----> Output result to Global Memory


*******************************************************************************/
#include <hls_stream.h>
#include <string.h>
#include <ap_int.h>

#define STAGES 4 //Number of Aadder Stages

//read_input(): Read Data from Global Memory and write into Stream inStream
static void read_input(int *input, hls::stream<char> &inStream , int size)
{
    mem_rd: for (int i = 0 ; i < size ; i++){
#pragma HLS pipeline
#pragma HLS LOOP_TRIPCOUNT min=4096 max=4096
        //Blocking write command to inStream 
        inStream << input[i];
    }
}

//adder(): Read Input data from inStream and write the result into outStream
//
static void adder(hls::stream<int> &inStream , hls::stream<int> &outStream, int incr, int size)
{
    execute: for (int i = 0 ; i < size ; i ++){
#pragma HLS pipeline
#pragma HLS LOOP_TRIPCOUNT min=4096 max=4096
        int inVar = inStream.read();
        int adderedVal = inVar + incr;
        //Blocking read command from inStream and Blocking write command 
        //to outStream 
        outStream << adderedVal;
    }
}

// obiviously this does not work 
static void indel_stream (hls::stream<char> conStreams[CON_SIZE], hls::stream<char> readsStreams[READS_SIZE]){
    char conShift[CON_SIZE][READS_SIZE];
    #pragma HLS ARRAY_PARTITION variable=conShift complete dim=0
    char readsShift[READS_SIZE][CON_SIZE];
    #pragma HLS ARRAY_PARTITION variable=readsShift complete dim=0
    //ap_shift_reg<char, CON_SIZE> conShift;
    //ap_shift_reg<char, READS_SIZE> readsShift;

    init_con_loop: 
    for (int i = 0; i < CON_SIZE; i++) {
        for (int j = 0; j < READS_SIZE; j++) {
            conShift[i][j] = conStreams[i].reads();
        }
    }

    init_read_loop: 
    for (int i = 0; i < READS_SIZE; i++) {
        for (int j = 0; j < CON_SIZE; j++) {
            readsShift[i][j] = readStreams[i].read();
        } 
    }
}

//write_result(): Read result from outStream and write the result to Global Memory
static void write_result(int *output, hls::stream<int> &outStream , int size)
{
    mem_wr: for (int i = 0 ; i < size ; i++){
#pragma HLS pipeline
#pragma HLS LOOP_TRIPCOUNT min=4096 max=4096
        //Blocking read command from OutStream  
        output[i] = outStream.read();
    }
}

//----------------------------------------- INDEL functions -------------------------------------------//

/*
int min_whd [NUM_CON][NUM_READ];

void whd (int min_whd[NUM_CON][NUM_READ]) {
        // JENNY is consensus_length read_length know at compile time? 
        for (int i = 0; i < NUM_CON; i++) {
            for (int j = 0; j < NUM_READ; j++) {
                //int whd_local [CON_LEN - READ_LEN];
                int min; 
                for (int k = 0; k <= consensus_length - read_length; k++) {
                    // this function? 
             
                    int ret = cal_whd(i, j, starting_index);
                    min = (ret < min) ? ret : min; 
                }
                min_whd[i][j] = min;
            }
        }
}

void score_whd (int min_whd[NUM_CON][NUM_READ], int scores[NUM_CON]) {

        // JENNY is consensus_length read_length know at compile time? 
        int min_idx;
        for (int i = 0; i < NUM_CON; i++) {
            for (int j = 0; j < NUM_READ; j++) {
                tmp = min_whd[i][j] - min_whd[0][j];
                // JENNY sum? 
                scores[i] += (tmp > 0) ? tmp: -tmp;
            }
        }
}*/

extern "C" {
//whd(con_arr, con_size, con_len, reads_arr, reads_size, reads_len, weights_arr, min_whd, min_whd_idx);
void Indel_Accel (char* consensus, const int consensus_size, int* consensus_length, \
    char* reads, const int reads_size, int* reads_length, char* qs, int* min_whd_idx) {
#pragma HLS INTERFACE m_axi port=consensus offset=slave bundle=gmem
#pragma HLS INTERFACE m_axi port=consensus_length offset=slave bundle=gmem
#pragma HLS INTERFACE m_axi port=reads offset=slave bundle=gmem
#pragma HLS INTERFACE m_axi port=reads_length offset=slave bundle=gmem
#pragma HLS INTERFACE m_axi port=qs offset=slave bundle=gmem
#pragma HLS INTERFACE m_axi port=min_whd_idx offset=slave bundle=gmem

#pragma HLS INTERFACE s_axilite port=consensus bundle=control
#pragma HLS INTERFACE s_axilite port=consensus_size bundle=control
#pragma HLS INTERFACE s_axilite port=consensus_length bundle=control
#pragma HLS INTERFACE s_axilite port=reads bundle=control
#pragma HLS INTERFACE s_axilite port=reads_size bundle=control
#pragma HLS INTERFACE s_axilite port=reads_length bundle=control
#pragma HLS INTERFACE s_axilite port=qs bundle=control
#pragma HLS INTERFACE s_axilite port=min_whd_idx bundle=control
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
    int min_whd[CON_SIZE][READS_SIZE];
    int min_whd_idx[CON_SIZE][READS_SIZE];
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
            //assert(min_idx <= local_consensus_length - local_reads_length);
            
            min_whd[i * reads_size + j] = min;
            min_whd_idx[i * reads_size + j] = min_idx;
            reads_base += local_reads_length;
        }
        consensus_base += local_consensus_length;
    }

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
}

/*extern "C" {
void N_stage_Adders(int *input, int *output, int incr, int size)
{
#pragma HLS INTERFACE m_axi port=input offset=slave bundle=gmem
#pragma HLS INTERFACE m_axi port=output offset=slave bundle=gmem
#pragma HLS INTERFACE s_axilite port=input bundle=control
#pragma HLS INTERFACE s_axilite port=output bundle=control
#pragma HLS INTERFACE s_axilite port=incr bundle=control
#pragma HLS INTERFACE s_axilite port=size bundle=control
#pragma HLS INTERFACE s_axilite port=return bundle=control

    //array of stream declaration
    hls::stream<int> streamArray[STAGES+1];

    #pragma HLS dataflow
    //one read input unit for data read
    read_input(input,streamArray[0],size);
    compute_loop: for (int i = 0 ; i < STAGES ; i++){
    #pragma HLS UNROLL
        // total 4 units of adder(). each is compute vector addition 
        // and sending result to immediate next unit using stream
        // datatype
        adder(streamArray[i],streamArray[i+1],incr,size);
    }

    //one write result unit to write result back to global Memory
    write_result(output,streamArray[STAGES],size);
}
}*/
