#include "simulation.h"

int main(){
    int n = 1024; // code length
    int k =  307; // code dimension
    int c =    8; // CRC length
    
    // the upper bound of the quantized output
    // alphabet size in the code construction algorithm.
    int u = 256;
    // Signal to Noise Ratio (SNR) [dB]
    // used in code construction
    double cons_snr = 2.00; 
	
    // List size in SCL decoder
    int L = 32;   

    // simulation round of every snr.
    int rounds = 10000;

    abs_simu_ins* asi = asi_init(n, k, c, u, cons_snr, L);
    
    for(double snr = 1.00; snr < 3.10; snr+=0.25){
        asi_simu(snr, rounds, asi);
    }

    asi_dele(asi);
    return 0;
}  