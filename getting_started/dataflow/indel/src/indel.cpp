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
static void read_input(int *input, hls::stream<int> &inStream , int size)
{
    mem_rd: for (int i = 0 ; i < size ; i++){
#pragma HLS pipeline
#pragma HLS LOOP_TRIPCOUNT min=4096 max=4096
        //Blocking write command to inStream 
        inStream << input[i];
    }
}

//adder(): Read Input data from inStream and write the result into outStream
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
// Target file representation 
// target number| chromo contig # |starting position | ending position
// offset is added to the starting position to produce the result
// 
// returns two values
// true,false
// if true update the position to the new k + starting position

int cal_whd(int i, int j, int starting_index){
    int whd = 0;
    int min = 0; 
    int k = 0; 
    // Is this correct?
    for (int i = 0; i <= consensus_length - read_length; i++) {
            
            // Can make it and? for four bit rep 
            for (int j = 0; j < read_length; j++) {
                    if (consensus[i] != read[j]){
                        whd += qs[j];
                    }
            }
            // kepp the min 
            
            
    }

    return whd;
}

#define NUM_CON 32
#define NUM_READ 256
// These are variables 
#define CON_LEN 2048
#define READ_LEN 256

// starts with 33 - 126
#define QS_LEN 256


int min_whd [NUM_CON][NUM_READ];

void whd (int consensus[NUM_CON][CON_LEN], int read[NUM_READ][READ_LEN], int qs[NUM_READ][READ_LEN], \
    int min_whd[NUM_CON][NUM_READ], int min_idx[NUM_CON][NUM_READ]) {
    for (int i = 0; i < NUM_CON; i++) {
        for (int j = 0; j < NUM_READ; j++) {
            //int whd_local [CON_LEN - READ_LEN];
            int min = 0xffff_ffff; 
            int min_idx = CON_LEN - READ_LEN + 1;
            //for (int k = 0; k <= consensus_length - read_length; k++) {
            for (int k = 0; k <= CON_LEN - READ_LEN; k++) {

                // whd 
                int whd = 0;
                for (int l = 0; l <= read_length; l++) {
                    if (consensus[k + l] != read[l]){
                        whd += qs[l];
                    }                        
                }
                if (whd < min) {
                    min =  whd; 
                    min_idx = k; 
               }

            }
            assert(min_idx <= CON_LEN - READ_LEN);
            min_whd[i][j] = min;
            min_idx[i][j] = min_idx;
        }
    }
}

void score_whd (int min_whd[NUM_CON][NUM_READ], int scores[NUM_CON]) {
    // might need to reduce the bits used in here 
    int min_score = 0xffff_ffff;
    int min_idx = NUM_CON + 1;
    for (int i = 1; i < NUM_CON; i++) {
        int score = 0;
        for (int j = 0; j < NUM_READ; j++) {
            tmp = min_whd[i][j] - min_whd[0][j];
            score += (tmp > 0) ? tmp: -tmp;
        }
        min_idx = (score < min_score) ? i : min_idx;
    }
    assert(min_idx < NUM_CON);
    return min_idx;
}

void indel (int consensus[NUM_CON][CON_LEN], int read[NUM_READ][READ_LEN], int qs[NUM_READ][READ_LEN] ){
    int min_whd[NUM_CON][[NUM_READ];

    


}


extern "C" {
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
}
