# ABS Polar Codes and ABS+ Polar Codes
This is an implementation of ***Adjacent-Bits-Swapped (ABS) Polar Codes*** ([arXiv](https://arxiv.org/abs/2202.04454)) and ***ABS+ Polar Codes: Exploiting More Linear Transforms on Adjacent Bits*** ([arXiv]()) proposed by Guodong Li, Min Ye and Sihuang Hu . This repo includes code construction, encoding and decoding (SCL and CRC-Aided SCL) of these two code families.

## Build

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
+ `rounds` - simulation round of every snr.
+ `u` - the upper bound of the quantized output alphabet size in the code construction algorithm.
+ `cons_snr` - signal to noise Ratio (SNR) [dB] used in code construction. In our paper, we chose cons_snr to be 2.00 dB in the code construction for all choices of n and k.


All parameters  to use in simulation are written in `main.cpp` . For example, 

```C++
#include "simulation.h"

#define cons_from_file
#define ABS_Polar

int main(){
    int n = 1024; // code length
    int k =  512; // code dimension
    int c =    0; // CRC length
    
    // List size in SCL decoder
    int L =   20;   

    // simulation round of every snr.
    int rounds = 10000;

    simu_ins* ins = ins_init(n, k, c, L);
    #ifdef cons_from_file
        #ifdef ABS_Polar
            char* consfile = (char*)"consfile/ABS_BPSK-AWGN_2.000dB_n1024_R0.5_u250000.txt";
            construct_abs_from_file(n, k, c, consfile, ins->I, ins->transform, ins->state);
        #elif
            char* consfile = (char*)"consfile/ABS_Plus_BPSK-AWGN_2.000dB_n1024_R0.5_u250000.txt";
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
### Macros
+ `cons_from_file` - If we define this macro, the main function will construct ABS(+) polar code from given construction file in the folder `consfile`.
                     Otherwise, If we conmment out this macro, then the main function will construct the code using the function `construct_abs` or the function `construct_abs_plus`.
+ `ABS_Polar` -  If we define this macro, the main function will construct ABS polar code. Otherwise, If we conmment out this macro, the main function will construct ABS+ polar code.

### The folder `consfile`

Since the running time of the code construction algorithm is very long, we provide the outputs of the code construction algorithm in the folder `consfile` for various choices of parameters. The file name in the folder `consfile` contains the parameters of the constructed codes.
For example, the file `consfile/ABS_Plus_BPSK-AWGN_2.00dB_n1024_R0.5_u250000.txt` includes the following parameters:

+ `ABS_Plus` - ABS+ polar code
+ `BPSK-AWGN`, `2.00dB`, `R0.5`-  mean that the code is constructed for the BPSK-AWGN channel with `SNR`=2.00dB and rate $R = 0.5$
+ `u250000` - the upper bound of the quantized output alphabet size is $\mu=250000$
+ `n1024` - the code length is $n = 1024$

### implementation

Since the ABS polar codes are special cases of the ABS+ polar codes, the decoding algorithm of ABS+ polar codes can also be used to decode the ABS polar codes.
In this version, we only call the decoding function of the ABS+ polar codes, and we put the original ABS decoder separately in file `abs_interfaces.h`.


## Simulation Results

We set the upper bound of the quantized output alphabet size to be $\mu=250000$ in the code construction algorithm to obtain the simulation results in our paper. However, running the code construction algorithm with $\mu=250000$ takes up to more than one month on a personal computer. In our simulations, we used a server with 128 threads to reduce the running time to several hours.
As mentioned above, we save the outputs of the code construction algorithm to the folder `consfile` so that you can reproduce our simulation results.
Below we show the comparison between the performance of ABS+ polar codes, ABS polar codes and standard polar codes over binary-input AWGN channels.

In the figures below, "ST" refers to standard polar codes, "ABS" refers to ABS polar codes, and "ABS+" refers to ABS+ polar codes.
The CRC length is chosen from the set ${4, 8, 12, 16, 20}$ to minimize the decoding error probability.

+ n = 1024, k = 307

<img src="/fig/1024_307.png?raw=true" alt="1024_307" title="Performance comparison between standard polar codes and ABS polar codes" style="zoom:100%;" />



+ n = 2048, k = 1024

<img src="/fig/2048_1024.png?raw=true" alt="2048_1024" title="Performance comparison between standard polar codes and ABS polar codes" style="zoom:100%;" />



+ n = 2048, k = 1434

<img src="/fig/2048_1434.png?raw=true" alt="2048_1434" title="Performance comparison between standard polar codes and ABS polar codes" style="zoom:100%;" />

## Acknowledgement

In the implementation of our decoding algorithm, we have learned a lot from the GitHub project [ecclab](https://github.com/kshabunov/ecclab)  maintained by Kirill Shabunov. Shabunov’s GitHub project mainly presents the implementation of the Reed-Muller decoder proposed in [3]. Due to the similarity between (ABS) polar codes and Reed-Muller codes, some of the accelerating techniques for Reed-Muller decoders can also be used to speed up (ABS) polar decoders.

## References

1. E. Arıkan, “Channel polarization: A method for constructing capacity-achieving codes for symmetric binary-input memoryless channels,” IEEE Transactions on Information Theory, vol. 55, no. 7, pp. 3051–3073, 2009. [IEEE Xplore](https://ieeexplore.ieee.org/document/5075875)
2. G. Li, Y. Me and S. Hu, "Adjacent-Bits-Swapped Polar codes: A new code construction to speed up polarization,"  [arXiv](https://arxiv.org/abs/2202.04454)
3. G. Li, Y. Me and S. Hu, "ABS+ Polar Codes: Exploiting More Linear Transforms on Adjacent Bits," [arXiv]()
4. I. Dumer and K. Shabunov, “Soft-decision decoding of Reed-Muller codes: Recursive lists,” IEEE Transactions on Information Theory, vol. 52, no. 3, pp. 1260–1266, 2006. [IEEE Xplore](https://ieeexplore.ieee.org/document/1603792)
5. I. Tal and A. Vardy, “How to construct polar codes,” IEEE Transactions on Information Theory, vol. 59, no. 10, pp. 6562–6582, 2013. [IEEE Xplore](https://ieeexplore.ieee.org/document/6557004)
6. I. Tal and A. Vardy, “List decoding of polar codes,” IEEE Transactions on Information Theory, vol. 61, no. 5, pp. 2213– 2226, 2015. [IEEE Xplore](https://ieeexplore.ieee.org/document/7055304)
