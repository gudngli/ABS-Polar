//==============================================================
//
// Interfaces of simulations.
//  
// Copyright 2022 and onwards Guodong Li
// 
// Licensed under the Apache License, Version 2.0.
// 
//==============================================================

#ifndef SIMULATION_INCLUDE
#define SIMULATION_INCLUDE

#include <stdio.h>

#include "construction.h"
#include "channel.h"
#include "encoding.h"
#include "decoding.h"

#include "output.h"

/**
 * Only output the final result after every simulation, if define FINAL,
 * update the intermediate results in each simulation,  else.
 */
#define FINAL

typedef struct simulation_instance{
    int n;
    int k;
    int c;

    int m; // n = 2^m
    double rate;

    int  *I;
    int **transform;
    int **state;

    int L;

    double snr;
    int         rounds;
    int    fail_rounds;
    int ML_fail_rounds;

}simu_ins;

simu_ins* ins_init(int n, int k, int c, int L);
void      ins_dele(simu_ins* ins);
void      ins_simu(double snr, int rounds, simu_ins* ins);

//================simulation.c============================

simu_ins* ins_init(int n, int k, int c, int L){
    simu_ins* ins = MALLOC(1, simu_ins);
    ins->n = n;
    ins->k = k;
    ins->c = c;

    ins->m = LOG2(n);
    ins->rate = ((double)k)/((double)n);


    ins->I = MALLOC(n, int);
    ins->transform = MALLOC(ins->m, int*);
    ins->state = MALLOC(ins->m, int*);
    
    ins->L = L;

    return ins;
}

void ins_dele(simu_ins* ins){
    FREE(ins->I);
    for(int i = 0; i < ins->m; i++){
        FREE(ins->transform[i]);
        FREE(ins->state[i]);
    }
    FREE(ins->transform);
    FREE(ins->state);
}

void ins_simu(double snr, int rounds, simu_ins* ins){
    ins->snr = snr;
    ins->rounds = rounds;
    double ssl = snr_sqrt_linear(snr, ins->rate);
    
    CRC*     crc = crc_init(ins->k, ins->c);
    encoder* enc = enc_init(ins->m, ins->k+ins->c, ins->I, ins->transform);
    decoder* dec = dec_init(ins->m, ins->L, ins->I, ins->state);

    int* message_crc = MALLOC(ins->n, int);
    int* codeword    = MALLOC(ins->n, int);
    double* llr      = MALLOC(ins->n, double);
    int* deResult    = MALLOC(ins->n, int);

    int   failRounds = 0;
    int MLfailRounds = 0;
    int    err;
    int ML_err;
    
    printf("Simulation results of (%d, %d, %d)-ABS(+) Polar code, L = %d, snr = %.2fdB:\n", ins->n, ins->k, ins->c, ins->L, snr);
            printf("    Rounds     error      MLerr        FER       MLFER \n");

    for(int r = 0; r < rounds; r++){
        
        get_message(message_crc, ins->k);
        set_crc(message_crc, crc);
        encode(message_crc, codeword, enc);
        transmit(codeword, llr, ssl, ins->n);
        
        decode(llr, deResult, crc, enc, dec);

        error(codeword, deResult, llr, ins->n, &err, &ML_err);
          failRounds +=    err;
        MLfailRounds += ML_err;
        #ifndef FINAL
        if(r%1000==999){
            printf("    %-10d %-10d %-10d %.2e  %.2e.\n", r+1, failRounds, MLfailRounds, ((double)failRounds)/((double)(r+1)), ((double)MLfailRounds)/((double)(r+1)));
        }
        #endif
    }

#ifdef FINAL
    printf("    %-10d %-10d %-10d %.2e  %.2e.\n", rounds, failRounds, MLfailRounds, ((double)failRounds)/((double)(rounds)), ((double)MLfailRounds)/((double)(rounds)));
#endif

    
    crc_dele(crc);
    enc_dele(enc);
    dec_dele(dec);

    FREE(message_crc);
    FREE(codeword);
    FREE(llr);
    FREE(deResult);

    ins->fail_rounds = failRounds;
    ins->ML_fail_rounds = MLfailRounds;
    
    printf("Simulation End.\n\n");
}

#endif // #ifndef SIMULATION_INCLUDE