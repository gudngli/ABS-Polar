//==============================================================
//
// Interfaces of code construction of ABS Polar Codes.
//  
// Copyright 2022 and onwards Guodong Li
// 
// Licensed under the Apache License, Version 2.0.
// 
//==============================================================
#ifndef CONSTRUCTION_INCLUDE
#define CONSTRUCTION_INCLUDE

#include "channel.h"

// I is the infotmation bits mask, a bit vector of length n.
// I[i] = 1 iff i-th bit is a information bit.
// 
// permutation is a m x n matrix recording the swap information
// after every polar transform in layer 0 ... m-1.
// 
// swap is a dynamic 2-D array which will be used in decoding of
// ABS polar codes. For all 1 <= i < m, swap[i] is a bit vector of length
// 2^(m-i) - 1.
// swap[i][j] = 1 iff the order of decoding of two bits of
// double bits (DB) channel V_{2i+1}^{(n)} will be swaped,
// where n = 2^{i+1}.
// 
// Different from our paper, all of the subscripts used here start from 0.
// 
// u is the upper bound of the quantized output alphabet size in the code
// construction algorithm.
void construct_abs(int n, int k, int c, double snr, int* I, int** permutation, int** swap, int u);

//====================construction.c================================

#include "transmat.h"
#include <future>
#define file_name_length 300

void wdele(double**w, int size);

// Search for non adjacent subsequences (which may be empty) with the largest sum
int* max_sum_subseq(int n, double* seq);

// swap[cur_m] to permutation[cur_m-1];
int* swap_to_permutation(int number, int* swap, int n);

typedef struct index_value{ int index; double value; }index_value;
int  index_cmp(const void *_a, const void *_b);
void choose_information_bits(int n, int k, double* conditional_entropy, int* I);

double* construction(double** w, int size, int m, int n, int** permutation, int** swap, int u);

void construct_abs(int n, int k, int c, double snr, int* I, int** permutation, int** swap, int u){
    printf("Constructing (%d, %d, %d)-ABS-Polar Code on BPSK-AWGN %.2fdB, upper bound of output alphabet size: %d...\n", n, k, c, snr, u);
    int m = LOG2(n);
    int size = 256; // original output alphabet size of BPSK-AWGN channel.
    double rate = ((double)k)/((double)n);
    double ssl  = snr_sqrt_linear(snr, rate);
    double** w  = discretization(ssl, size);
    double* arr = construction(w, size, m, n, permutation, swap, u);

    choose_information_bits(n, k+c, arr, I);

    FREE(arr);
    wdele(w, size);
    printf("End of Construction.\n");
}


double* construction(double** w, int size, int m, int n, int** permutation, int** swap, int u){
    int d = (int)pow((double)u, 1.0/3.0)+1;
    transmat*** matrix = MALLOC(m, transmat**);
    
    int cur_m, number;
    for(cur_m = 0; cur_m < m; cur_m++){
        number = (1 << (m-cur_m))-1;
        matrix[cur_m] = MALLOC(number, transmat*);
        for(int branch = 0; branch < number; branch++)matrix[cur_m][branch] = NULL;
    }
    cur_m  = m-1;
    number =   1;
    merge_map* map = map_init(d);
    matrix[cur_m][0] = trmat_init(w, size, map);//INFO;
    permutation[cur_m] = NULL;
    
    // profit (increase of degree of polarization) of switch order of decoding of two channels
    // computing for profit of a combine channel see the definition of function profit(transmat* trmat);
    double* prof = MALLOC(n, double);
    transmat**   help = MALLOC(n, transmat*);
    transmat**  phelp = MALLOC(n, transmat*);
        for(int i = 0; i < n; i++){help[i] = phelp[i] = NULL;}

    while(cur_m){
        // old branches of layer cur_m-1
        for(int i = 0; i < number; i++){
            help[i] =   trans1(matrix[cur_m][i], map);//INFO;
           phelp[i] = p_trans1(matrix[cur_m][i], map);//INFO;
            prof[i] = profit(help[i], phelp[i]);
        }
        swap[cur_m] = max_sum_subseq(number, prof);
        permutation[cur_m-1] = swap_to_permutation(number, swap[cur_m], n);
        for(int i = 0; i < number; i++){
            if (swap[cur_m][i]){
                matrix[cur_m-1][(i<<1)+1] = phelp[i];
                trmat_dele( help[i]);
            }else{
                matrix[cur_m-1][(i<<1)+1] =  help[i];
                trmat_dele(phelp[i]);
            }
            help[i] = phelp[i] = NULL;
        }
        // even branches of layer cur_m-1
        for(int i = 0; i < number; i++){
            if (swap[cur_m][i]){
                matrix[cur_m-1][(i<<1)  ] = p_trans0(matrix[cur_m][i], map);//INFO;
                matrix[cur_m-1][(i<<1)+2] = p_trans2(matrix[cur_m][i], map);//INFO;
            }else{
                if (i==0 || swap[cur_m][i-1]==0){
                    matrix[cur_m-1][(i<<1)  ] = trans0(matrix[cur_m][i], map);//INFO;
                }
                if (i==number-1){
                    matrix[cur_m-1][(i<<1)+2] = trans2(matrix[cur_m][i], map);//INFO;
                }
            }
        }

        for(int i = 0; i < number; i++)trmat_dele(matrix[cur_m][i]);
        cur_m--;
        number = (1<<(m-cur_m))-1;
    }
    
    swap[0] = NULL;
    map_dele(map);

    double* conditional_entropy = MALLOC(n, double);
    for (int i = 0; i < n-1; i++){
        conditional_entropy[i] = minus_conditional_entropy(matrix[0][i]);
    }
    conditional_entropy[n-1] = plus_conditional_entropy(matrix[0][n-2]);

    for(int i = 0; i < number; i++)trmat_dele(matrix[0][i]);
    FREE( help);
    FREE(phelp);
    FREE(prof);
    for(int i = 0; i < m; i++){
        FREE(matrix[i]);
    }
    FREE(matrix);

    return conditional_entropy;
}



void wdele(double**w, int size){
    for (int i = 0; i < size; i++){
        FREE(w[i]);
    }
    FREE(w);
}

int* max_sum_subseq(int n, double* seq){
    double sum0, sum1;
    int *I0 = MALLOC(n, int), *I1 = MALLOC(n, int);
    
    //init
    sum0 = sum1 = 0.0;
    for (int i = 0; i < n; i++){
        I0[i] = I1[i] = 0;
    }
    
    // seq[0]
    if (GT(seq[0],0)){ I0[0] = 1; sum0 = seq[0]; }
    // seq[1]
    if (GT(seq[1], sum0)){ I1[0] = 0; I1[1] = 1; sum1 = seq[1]; }
    else{ I1[0] = I0[0]; sum1 = sum0; }

    for (int i = 2; i < n; i++){
        //seq[i]
        if (GT(seq[i]+sum0, sum1)){
            I0[i] = 1;
            sum0+=seq[i];
            int   *tempI =   I0;   I0 =   I1;   I1 = tempI;  // swap   I0 and   I1
            double temps = sum0; sum0 = sum1; sum1 = temps;  // swap sum0 and sum1
        }else{
            sum0 = sum1;
            for (int j = 0; j <= i; j++)I0[j] = I1[j];
        }
    }
    FREE(I0);
    return I1;
}

int* swap_to_permutation(int number, int* swap, int n){
    int flag = 0;
    for (int i = 0; i < number; i++){
        if (swap[i]==1){
            flag = 1;
            break;
        }
    }
    if (!flag)return NULL;

    int  cur_n = (number<<1) + 2;
    int* prepermu = MALLOC(cur_n, int);
    for (int i = 0; i < cur_n;  i++) prepermu[i] = i;
    for (int i = 0; i < number; i++){
        if (swap[i]==1){
            // swap prepermu[2i+1] and prepermu[2i+2]
            prepermu[(i<<1)+1] = (i<<1)+2;
            prepermu[(i<<1)+2] = (i<<1)+1;
        }
    }
    int* permutation = MALLOC(n, int);
    int factor = n / cur_n;
    for (int i = 0; i < cur_n; i++){
        for (int j = 0; j < factor; j++){
            permutation[factor*i+j] = factor * prepermu[i] + j;
        }
    }
    FREE(prepermu);
    return permutation;
}

int index_cmp(const void *_a, const void *_b){
    index_value* a = (index_value*) _a;
    index_value* b = (index_value*) _b;
    if (LT(a->value, b->value)) return -1;
    if (GT(a->value, b->value)) return  1;
    return 0;
}

void choose_information_bits(int n, int k, double* conditional_entropy, int* I){
    index_value* iv = MALLOC(n, index_value);
    for (int i = 0; i < n; i++){
        iv[i].index = i;
        iv[i].value = conditional_entropy[i];
    }
    qsort(iv, n, sizeof(index_value), index_cmp);
    for (int i = 0; i < n; i++)I[i] = 0;

    for (int i = 0; i < k; i++){
        I[iv[i].index] = 1;
    }
    FREE(iv);
}



#endif // #ifndef CONSTRUCTION_INCLUDE