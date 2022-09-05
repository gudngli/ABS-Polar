//==============================================================
//
// Interfaces of ABS(+) Polar Codes encoder.
//  
// Copyright 2022 and onwards Guodong Li
// 
// Licensed under the Apache License, Version 2.0.
// 
//==============================================================

#ifndef ENCODING_INCLUDE
#define ENCODING_INCLUDE

#include "common.h"

typedef struct encoder{
    int   m;
    int   n; //n = 2^m;
    int   k;
    int*  I; // There are k bits in I equal 1.
    int** transform;
}encoder;

encoder* enc_init(int m, int k, int *I, int ** transform);
void     enc_dele(encoder* enc);

// message vector to codeword
void         encode(int* message, int* codeword, encoder* enc);
// codeword to messagr vector
void inverse_encode(int* message, int* codeword, encoder* enc);

//===============================encoding.c========================

encoder* enc_init(int m, int k, int *I, int ** transform){
    encoder* enc = MALLOC(1, encoder);
    enc->m = m;
    enc->n = 1<<m;
    enc->k = k;
    enc->I = I;
    enc->transform = transform;
    return enc;
}

void enc_dele(encoder* enc){
    FREE(enc);
}



void encode(int* message, int* codeword, encoder* enc){
    for(int i = enc->n-1, j = enc->k-1; i>=0; i--){
        codeword[i] = (enc->I[i])?(message[j--]):(0);
    }
    for(int i = 0; i < enc->m; i++){
        // transform
        if (enc->transform[i]){
            for(int j = 0; j < enc->n; j++){
                int pj = enc->transform[i][j];
                if (pj < 0) {
                    // Arikan transform
                    codeword[j]^=codeword[-pj];
                }else if (pj > j){
                    // swapping transform
                    // swap codeword[j] and codeword[pj]
                    codeword[ j] ^= codeword[pj];
                    codeword[pj] ^= codeword[ j];
                    codeword[ j] ^= codeword[pj];
                }
            }
        }
        // polar transform
        for(int j = 0; j < enc->n; j++){
            if (((j>>i)&1)==0){
                codeword[j]^=codeword[j^(1<<i)];
            }
        }
    }
}


void inverse_encode(int* message, int* codeword, encoder* enc){
    for(int i = enc->m-1; i >= 0; i--){
        // polar transform
        for(int j = 0; j < enc->n; j++){
            if (((j>>i)&1)==0){
                codeword[j]^=codeword[j^(1<<i)];
            }
        }
        // transform
        if (enc->transform[i]){
            for(int j = 0; j < enc->n; j++){
                int pj = enc->transform[i][j];
                if (pj < 0) {
                    // Arikan transform
                    codeword[j]^=codeword[-pj];
                }else if (pj > j){
                    // swapping transform
                    // swap codeword[j] and codeword[pj]
                    codeword[ j] ^= codeword[pj];
                    codeword[pj] ^= codeword[ j];
                    codeword[ j] ^= codeword[pj];
                }
            }
        }
    }
    for(int i = 0, j = 0; i < enc->n; i++){
        if (enc->I[i]){
            message[j++] = codeword[i];
        }
    }
}

#endif // #ifndef ENCDOING_INCLUDE