# ABS Polar Codes and ABS+ Polar Codes
This is an implementation of ***Adjacent-Bits-Swapped (ABS) Polar Codes*** ([arXiv](https://arxiv.org/abs/2202.04454)) and ***ABS+ Polar Codes: Exploiting More Linear Transforms on Adjacent Bits*** ([arXiv]()) proposed by Guodong Li, Min Ye and Sihuang Hu . This repo includes code construction, encoding and decoding (SCL and CRC-Aided SCL) of these two families of codes.

## Build

ABS+ Polar Codes: Exploiting More Linear Transforms on Adjacent Bits

Only the randomization part uses C++ library `<random>`.  All the other parts of this repo are implemented in C. We don't separate the implementations and definitions of  interfaces, so you can compile this program on the command line easily.

```
$ g++ main.cpp -o main
```

## Running

```
$ ./main
```

### Simulation parameters

+ `n` - code length.
+ `k` - code dimension
+ `c` - CRC length
+ `L` - List size in SCL decoder
+ `u` - the upper bound of the quantized output alphabet size in the code construction algorithm.
+ `cons_snr` - signal to noise Ratio (SNR) [dB] used in code construction. In our paper, we chose cons_snr to be 2.00 dB in the code construction for all choices of n and k.
+ `rounds` - simulation round of every snr.

All parameters  to use in simulation are written in `abs.cpp` . For example, 

```C++
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
```

## Simulation Results

We set the upper bound of the quantized output alphabet size to be $\mu=250000$ in the code construction algorithm to obtain the simulation results in our paper. However, running the code construction algorithm with $\mu=250000$ takes up to more than one month on a personal computer. In our simulations, we used a server with 128 threads to reduce the running time to several hours. We recommend setting the parameter $\mu$ to be 8000 when running the code construction algorithm on a personal computer. When $\mu=8000$, the construction only takes 3 to 4 hours, and the resulting codes have essentially the same performance as the ones constructed in our paper using $\mu=250000$. Below we show the comparison between the performance of ABS polar codes and standard polar codes over binary-input AWGN channels. 

In the figures below, "ST" refers to standard polar codes, and "ABS" refers to ABS polar codes.

+ n = 2048, k = 614

<img src="/fig/2048_614.png?raw=true" alt="2048_614" title="Performance comparison between standard polar codes and ABS polar codes" style="zoom:100%;" />



+ n = 2048, k = 1024

<img src="/fig/2048_1024.png?raw=true" alt="2048_1024" title="Performance comparison between standard polar codes and ABS polar codes" style="zoom:100%;" />



+ n = 2048, k = 1434

<img src="/fig/2048_1434.png?raw=true" alt="2048_1434" title="Performance comparison between standard polar codes and ABS polar codes" style="zoom:100%;" />

## Acknowledgement

In the implementation of our decoding algorithm, we have learned a lot from the GitHub project [ecclab](https://github.com/kshabunov/ecclab)  maintained by Kirill Shabunov. Shabunov’s GitHub project mainly presents the implementation of the Reed-Muller decoder proposed in [3]. Due to the similarity between (ABS) polar codes and Reed-Muller codes, some of the accelerating techniques for Reed-Muller decoders can also be used to speed up (ABS) polar decoders.

## References

1. E. Arıkan, “Channel polarization: A method for constructing capacity-achieving codes for symmetric binary-input memoryless channels,” IEEE Transactions on Information Theory, vol. 55, no. 7, pp. 3051–3073, 2009. [IEEE Xplore](https://ieeexplore.ieee.org/document/5075875)
2. G. Li, Y. Me and S. Hu, "Adjacent-Bits-Swapped Polar codes: A new code construction to speed up polarization,"  [arxiv](https://arxiv.org/abs/2202.04454)
3. I. Dumer and K. Shabunov, “Soft-decision decoding of Reed-Muller codes: Recursive lists,” IEEE Transactions on Information Theory, vol. 52, no. 3, pp. 1260–1266, 2006. [IEEE Xplore](https://ieeexplore.ieee.org/document/1603792)
4. I. Tal and A. Vardy, “How to construct polar codes,” IEEE Transactions on Information Theory, vol. 59, no. 10, pp. 6562–6582, 2013. [IEEE Xplore](https://ieeexplore.ieee.org/document/6557004)
5. I. Tal and A. Vardy, “List decoding of polar codes,” IEEE Transactions on Information Theory, vol. 61, no. 5, pp. 2213– 2226, 2015. [IEEE Xplore](https://ieeexplore.ieee.org/document/7055304)
