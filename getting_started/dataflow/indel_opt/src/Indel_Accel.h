#ifndef INDEL_ACCEL_H
#define INDEL_ACCEL_H

//#define CON_SIZE 2
//#define READS_SIZE 8 
#define CON_SIZE 32
#define READS_SIZE 256
#define CON_LEN 2048
#define READS_LEN 256

int* parse_schedule (const char* file, int* num_tests);
void parse(const char* file_prefix, int col_num, char* con_arr, int** con_len_arr, int* con_size, int*total_len);
#endif

