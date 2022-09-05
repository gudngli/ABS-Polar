#include "simulation.h"

#define cons_from_file
#define ABS_Polar

int main(){
    int n = 64; // code length
    int k = 19; // code dimension
    int c =  0; // CRC length
    
    // List size in SCL decoder
    int L = 32;   

    // simulation round of every snr.
    int rounds = 10000;

    simu_ins* ins = ins_init(n, k, c, L);
    #ifdef cons_from_file
        #ifdef ABS
            char* consfile = (char*)"consfile/ABS_BPSK-AWGN_2.000dB_n64_R0.3_u250000.txt";
            construct_abs_from_file(n, k, c, consfile, ins->I, ins->transform, ins->state);
        #elif
            char* consfile = (char*)"consfile/ABS_Plus_BPSK-AWGN_2.000dB_n64_R0.3_u250000.txt";
            construct_abs_plus_from_file(n, k, c, consfile, ins->I, ins->transform, ins->state);
        #endif
    #else
        // the upper bound of the quantized output
        // alphabet size in the code construction algorithm.
        int u = 256;
        // Signal to Noise Ratio (SNR) [dB]
        // used in code construction
        double cons_snr = 2.00; 
        #ifdef ABS_Polar
            construct_abs(n, k, c, cons_snr, ins->I, ins->transform, ins->state, u);
        #else
            construct_abs_plus(n, k, c, cons_snr, ins->I, ins->transform, ins->state, u);
        #endif
        out_swap(stdout, (char*)"state array: ", ins->state, ins->m);
        out_bits(stdout, (char*)"info bits mask: ", ins->I, ins->n);
    #endif
    
    for(double snr = 1.00; snr < 3.10; snr+=0.25){
        ins_simu(snr, rounds, ins);
    }

    ins_dele(ins);
    return 0;
}  