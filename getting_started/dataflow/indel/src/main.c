#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include "indel_ref.c"

int* read_target (){

  FILE *fp;
  char* con = "./1025.TARGET.tbl"; 
  fp=fopen(con, "r");

  if (fp == NULL)
    exit(EXIT_FAILURE);

  char line[2046];
  char separators[] = "|";

  int* target = malloc(sizeof(int)*4);

  fgets(line, sizeof(line), fp);
  char *p = strtok(line, separators); 

  int i = 0; 

  while (p != NULL){
    target[i++] = atoi(p); 
    //printf("%s\n", p);
    p = strtok(NULL, separators);
  } 
  printf("%d %d\n", i, target[3]);
  fclose(fp);

  int target_num = target[0];
  int contig_num = target[1];
  int starting_pos = target[2];
  int ending_pos = target[3];

  int length = ending_pos - starting_pos;
  printf("length: %d\n", length);

  return target;
}

int count_lines(char* filename){

  FILE *fp;
  fp = fopen(filename, "r");

  // Check if file exists
  if (fp == NULL)
  {
    printf("Could not open file %s\n", filename);
    return 0;
  }
  int count = 0;
  char c; 
  // Extract characters from file and store in character c
  for (c = getc(fp); c != EOF; c = getc(fp))
    if (c == '\n') // Increment count if this character is newline
      count = count + 1;

  // Close the file
  fclose(fp);
  return count+1;
}

void parse(const char* file_prefix, int col_num, char* con_arr, int** con_len_arr, int* con_size) {

  FILE *fp;
  char con[256]="";
  strcat(con, file_prefix);
  strcat(con, ".tbl");

  int num_lines = count_lines(con);
  //printf("num_lines %d\n", num_lines);
  fp=fopen(con, "r");

  if (fp == NULL)
    exit(EXIT_FAILURE);

  // Length of each segment 
  int * con_len = (int*) malloc(sizeof(int) * num_lines);
  assert(con_len != NULL);
  * con_len_arr = con_len;
  // Location of the end of each segment  
  int* pos = (int*) malloc(sizeof(int) * num_lines);
  assert(pos != NULL);

  char separators[] = "|";
  char *line = NULL;
  size_t line_len = 0;
  ssize_t read;

  int line_num = 0;
  int base = 0;
  while ((read = getline(&line, &line_len, fp)) != -1) {
    //printf("Retrieved line of length %zu :\n", read);
    //printf("%s", line);

    char *p = strtok(line, separators); 

    assert (p != NULL);

    int i = 1; 
    while (p != NULL){
      //printf("%s\n", p);
      //target[i++] = atoi(p); 
      p = strtok(NULL, separators);

      if(i == col_num){ 
        //if(line_num == 0){
        //  len = strlen(p); 
        //  con_arr = (char*) malloc(len * num_lines * sizeof(char));
        //}
        con_len[line_num] = strlen(p);
        int k; 
        for (k = 0; k < con_len[line_num]; k++){
          //printf("%c", p[k]);
          strncpy(&con_arr[base + k], &p[k], 1);
          //printf("%d\n", line_num* len + k);
          printf("%c", con_arr[base + k]);
        }
        printf("\n");
      }
      i++;
    }
    base += con_len[line_num];
    pos[line_num] = base; 
    line_num++;
  }
  //printf("%d\n", len);
  //int arr_size = len * num_lines;
  //printf("%d\n", arr_size);
  //for(int i= 0; i < arr_size; i++ ){
  //  printf("%d ", con_arr[i]);
  //}

  *con_size = num_lines;
  free(line);
  fclose(fp); 
}


#define NUM_CON 32
#define NUM_READ 256
#define CON_LEN 2048
#define READ_LEN 256

int main() {
  const char* file_prefix = "./ir_toy/00";
  char con[256]="";
  strcat(con, file_prefix);
  strcat(con, ".SEQ");

  char reads[256]="";
  strcat(reads, file_prefix);
  strcat(reads, ".READS");

  char* con_arr, *reads_arr, *weights_arr; 
  int *con_len, con_size, *reads_len, reads_size; 

  // Malloc the largest array 
  con_arr = (char*) malloc( NUM_CON * CON_LEN * sizeof(char));
  reads_arr = (char*) malloc( NUM_READ * READ_LEN * sizeof(char));
  weights_arr = (char*) malloc( NUM_READ * READ_LEN * sizeof(char));

  parse( con,1, con_arr, &con_len, &con_size);
  parse( reads, 4, reads_arr, &reads_len, &reads_size);
  parse( reads, 5, weights_arr, &reads_len, &reads_size);
   
  int* min_whd = (int*) malloc(con_size * reads_size * sizeof(int));
  int* min_whd_idx = (int*) malloc(con_size * reads_size * sizeof(int));
  int* new_ref= (int*) malloc(reads_size * sizeof(int));
  int* new_ref_idx= (int*) malloc(reads_size * sizeof(int));

  whd(con_arr, con_size, con_len, reads_arr, reads_size, reads_len, weights_arr, min_whd, min_whd_idx);
  score_whd (min_whd,  min_whd_idx, con_size, reads_size, new_ref, new_ref_idx);
 
  return 0;

}

