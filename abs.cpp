#include "simulation.h"

int main(){
    int n = 1024;
    int k =  307;
    int c =    8;
    
    int u = 256;
    double cons_snr = 2.00; //dB

    int L = 32;

    int rounds = 10000;

    abs_simu_ins* asi = asi_init(n, k, c, u, cons_snr, L);
    
    for(double snr = 1.00; snr < 3.10; snr+=0.25){
        asi_simu(snr, rounds, asi);
    }

    asi_dele(asi);
    return 0;
}  