#ifndef OUTPUT_INCLUDE
#define OUTPUT_INCLUDE

#include <stdio.h>

#include "transmat.h"
#include "transform.h"

// transition matrix
void onebit(FILE* file, char* info,  double** w, int size){
    fprintf(file, "%s\n", info);
    for(int i = 0; i < size; i++){
        fprintf(file, "    %3d, %.6f, %.6f\n", i, w[i][0], w[i][1]);
    }
    fprintf(file, "\n");
}

void twobit(FILE* file, char* info, transmat* trmat){
    double** w = trmat->w;
    int size = trmat->size;
    fprintf(file, "%s\n", info);
    for(int i = 0; i < size; i++){
        fprintf(file, "    %3d, %.6f, %.6f, %.6f, %.6f\n",i,  w[i][0], w[i][1], w[i][2], w[i][3]);
    }
    fprintf(file, "\n");
}


void out_doubles(FILE* file, char* info, double* arr, int n){
    fprintf(file, "%s\n", info);
    for(int i = 0; i < n; i++){
        fprintf(file, "    %3d, %.6f\n",i, arr[i]);
    }
    fprintf(file, "\n");
}

void out_bits(FILE* file, char* info, int* arr, int n){
    fprintf(file, "%s ", info);
    for(int i = 0; i < n; i++){
        fprintf(file, "%d", arr[i]);
    }
    fprintf(file, "\n");
}

void out_permutatation(FILE* file, char* info, int** permutation, int m){
    fprintf(file, "%s\n", info);
    for(int i = 0; i < m; i++){
        if (permutation[i]){
            fprintf(file, "    layer %2d:  ", i);
            for(int j = 0; j < (1<<m); j++){
                fprintf(file, "%5d ", permutation[i][j]);
            }
            fprintf(file, "\n");
        }else{
            fprintf(file, "    layer %2d:  null\n", i);
        }
    }
}

void out_swap(FILE* file, char* info, int** swap, int m){
    fprintf(file, "%s\n", info);
    fprintf(file, "    layer %2d:   NULL\n", 0);
    int number;
    for(int i = 1; i < m; i++){
        number = (1<<(m-i))-1;
        fprintf(file, "    layer %2d:  ", i);
        for(int j = 0; j < number; j++){
            fprintf(file, "%2d ", swap[i][j]);
        }
        fprintf(file, "\n");
    }
}

#endif // #ifndef OUTPUT_INCLUDE