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
//#define FINAL

typedef struct abs_simu_ins{
    int n;
    int k;
    int c;

    int m; // n = 2^m
    double rate;
    
    int u;
    double cons_snr;

    int  *I;
    int **permutation;
    int **swap;

    int L;

    double snr;
    int         rounds;
    int    fail_rounds;
    int ML_fail_rounds;

}abs_simu_ins;

abs_simu_ins* asi_init(int n, int k, int c, int u, double cons_snr, int L);
void          asi_dele(abs_simu_ins* asi);
void          asi_simu(double snr, int rounds, abs_simu_ins* asi);

//================simulation.c============================

abs_simu_ins* asi_init(int n, int k, int c, int u, double cons_snr, int L){
    abs_simu_ins* asi = MALLOC(1, abs_simu_ins);
    asi->n = n;
    asi->k = k;
    asi->c = c;

    asi->m = LOG2(n);
    asi->rate = ((double)k)/((double)n);

    asi->u = u;
    asi->cons_snr = cons_snr;

    asi->I = MALLOC(n, int);
    asi->permutation = MALLOC(asi->m, int*);
    asi->swap = MALLOC(asi->m, int*);
    
    construct_abs(n, k, c, cons_snr, asi->I, asi->permutation, asi->swap, u);
    out_swap(stdout, (char*)"swap array: ", asi->swap, asi->m);
    out_bits(stdout, (char*)"info bits mask: ", asi->I, asi->n);
    asi->L = L;

    return asi;
}

void asi_dele(abs_simu_ins* asi){
    FREE(asi->I);
    for(int i = 0; i < asi->m; i++){
        FREE(asi->permutation[i]);
        FREE(asi->swap[i]);
    }
    FREE(asi->permutation);
    FREE(asi->swap);
}

void asi_simu(double snr, int rounds, abs_simu_ins* asi){
    asi->snr = snr;
    asi->rounds = rounds;
    double ssl = snr_sqrt_linear(snr, asi->rate);
    
    CRC*     crc = crc_init(asi->k, asi->c);
    encoder* enc = enc_init(asi->m, asi->k+asi->c, asi->I, asi->permutation);
    decoder* dec = dec_init(asi->m, asi->L, asi->I, asi->swap);

    int* message_crc = MALLOC(asi->n, int);
    int* codeword    = MALLOC(asi->n, int);
    double* llr      = MALLOC(asi->n, double);
    int* deResult    = MALLOC(asi->n, int);

    int   failRounds = 0;
    int MLfailRounds = 0;
    int    err;
    int ML_err;
    
    printf("Simulation results of (%d, %d, %d)-ABS-Polar code, L = %d, snr = %.2fdB:\n", asi->n, asi->k, asi->c, asi->L, snr);
            printf("    Rounds     error      MLerr        FER       MLFER \n");

    for(int r = 0; r < rounds; r++){
        
        get_message(message_crc, asi->k);
        set_crc(message_crc, crc);
        encode(message_crc, codeword, enc);
        transmit(codeword, llr, ssl, asi->n);
        
        decode(llr, deResult, crc, enc, dec);

        error(codeword, deResult, llr, asi->n, &err, &ML_err);
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

    asi->fail_rounds = failRounds;
    asi->ML_fail_rounds = MLfailRounds;
    
    printf("Simulation End.\n\n");
}

#endif // #ifndef SIMULATION_INCLUDE