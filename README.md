# ABS-Polar
This is an implementation of ***Adjacent-Bits-Swapped (ABS) Polar Code*** Proposed by Guodong Li, Min Ye, Sihuang Hu ([arXiv address]()) at 2022. This repo include construction, encoding and decoding (SCL and CA-SCL) of the new code.

## Build

Except that the randomization part uses C+ +, the rest of this repo is implemented in C language. We don't separate the implementations and definitions of  interfaces, so you can compile this program on the command line easily.

```
$ g++ abs.cpp -o abs
```

## Running

```
$ ./abs
```

All parameters to use in simulation are written in `abs.cpp` . For example

```c++
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
```



## Simulation Result

It should be noted that in order to obtain data faster and more accurately, we use the parallel implementation of the program when obtaining data. In particular,  we set the upper bound of output alphabet size $\mu = 250000$, and  we used 128 threads to test. However, in order to show algorithms easily, the parallel implementation of the program is omitted in the repository. Even if we set $\mu = 8000$ , we can get the similar results as in the paper.

(256， 77) ABS-Polar Code

![256_77](C:\Users\PlumJ\Desktop\Task\ABS-Polar\fig\256_77.png)

 (256， 128) ABS-Polar Code

![256_128](C:\Users\PlumJ\Desktop\Task\ABS-Polar\fig\256_128.png)

 (256， 179) ABS-Polar Code

![256_179](C:\Users\PlumJ\Desktop\Task\ABS-Polar\fig\256_179.png)

 (512, 154) ABS-Polar Code

![512_154](C:\Users\PlumJ\Desktop\Task\ABS-Polar\fig\512_154.png)

 (512, 256) ABS-Polar Code

![512_256](C:\Users\PlumJ\Desktop\Task\ABS-Polar\fig\512_256.png)

 (512, 358) ABS-Polar Code

![512_358](C:\Users\PlumJ\Desktop\Task\ABS-Polar\fig\512_358.png)

 (1024, 307) ABS-Polar Code

![1024_307](C:\Users\PlumJ\Desktop\Task\ABS-Polar\fig\1024_307.png)

 (1024, 512) ABS-Polar Code

![1024_512](C:\Users\PlumJ\Desktop\Task\ABS-Polar\fig\1024_512.png)

 (1024, 717) ABS-Polar Code

![1024_717](C:\Users\PlumJ\Desktop\Task\ABS-Polar\fig\1024_717.png)

 (2048, 614) ABS-Polar Code

![2048_614](C:\Users\PlumJ\Desktop\Task\ABS-Polar\fig\2048_614.png)

 (2048, 1024) ABS-Polar Code

![2048_1024](C:\Users\PlumJ\Desktop\Task\ABS-Polar\fig\2048_1024.png)

 (2048, 1434) ABS-Polar Code

![2048_1434](C:\Users\PlumJ\Desktop\Task\ABS-Polar\fig\2048_1434.png)

## Acknowledgement

In the implementation of our decoding algorithm, we have learned a lot from the GitHub project https://github.com/kshabunov/ecclab  maintained by Kirill Shabunov. Shabunov’s GitHub project mainly presents the implementation of the Reed-Muller decoder proposed in "Soft-decision decoding of Reed-Muller codes: Recursive lists". Due to the similarity between (ABS) polar codes and Reed-Muller codes, some of the accelerating techniques for Reed-Muller decoders can also be used to speed up (ABS) polar decoders.


## References

* E. Arikan, “Channel Polarization: A Method for Constructing Capacity-Achieving Codes for Symmetric Binary-Input Memoryless Channels,”
 IEEE Transactions on Information Theory, vol. 55, no. 7, July 2009.  [IEEE Xplore](https://ieeexplore.ieee.org/document/5075875)
* I. Tal and A. Vardy, "List Decoding of Polar Codes",
 IEEE Transactions on Information Theory, vol. 61, no. 5, May 2015. [IEEE Xplore](https://ieeexplore.ieee.org/document/7055304)
* I. Dumer and K. Shabunov, “Soft-decision decoding of Reed-Muller codes: Recursive lists,” IEEE Transactions on Information Theory, vol. 52, no. 3, pp. 1260–1266, 2006. [IEEE Xplore](https://ieeexplore.ieee.org/document/1603792)
* I. Tal and A. Vardy, “How to construct polar codes,” IEEE Transactions on Information Theory, vol. 59, no. 10, pp. 6562–6582, 2013 [IEEE Xplore](https://ieeexplore.ieee.org/document/6557004)
