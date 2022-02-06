//==============================================================
//
// Interfaces of BPSK-AWGN channel.
//  
// Copyright 2022 and onwards Guodong Li
// 
// Licensed under the Apache License, Version 2.0.
// 
//==============================================================


#ifndef CHANNEL_INCLUDE
#define CHANNEL_INCLUDE

#include "common.h"
#include "Random.h"
#include <math.h>

// Standard deviation of Gauss channel.
const double stddev = sqrt(0.5);

double snr_sqrt_linear(double snr, double rate);

// transmit a codeword through a BPSK-AWGN channel
void transmit(const int* codeword, double* llr, double snr_sqrt_linear, int n);

// return transition matrix after discretization of Gaussian channel,
// where u is the upper bound of output alphabet size.
double** discretization(double snr_sqrt_linear, uint u);


//==================================channel.c===========================================

double snr_sqrt_linear(double snr, double rate){
    return pow(10.0, snr/20.0)*sqrt(rate);
}

void transmit(const int* codeword, double* llr, double snr_sqrt_linear, int n){
    std::random_device rd;
    std::mt19937_64       rng{ rd() };
    std::normal_distribution<double>       std_Gauss_d(0.0f, 1.0f);
    double coefficient = 4.0 * snr_sqrt_linear;
    for (int i = 0; i < n; i++){
        llr[i] = coefficient * Gauss_Distribution((codeword[i]?(-1.0):(1.0))*snr_sqrt_linear,  stddev);
    }
}



//==================================================================
//  For discretization of continuous channels, refer to 
//  I. Tal and A. Vardy, "How to construct polar codes,"" 
//  IEEE Transactions on Information Theory, vol. 59, no. 10,
//  pp. 6562–6582, 2013. 
//      https://ieeexplore.ieee.org/document/6557004
//  Section VI. Channels with Continuous Output.
//==================================================================



// standart normal distribution CDF
probability PHI(double x);
// CDF of N(u, 0.5)
probability   F(double x, double u);
double     capacity_lambda(double y, double snr_sqrt_linear);
double arc_capacity_lambda(double C, double snr_aqrt_linear);

double** discretization(double snr_sqrt_linear, uint u){
    uint v = u/2;
    double *y = MALLOC(v,double);
    for (uint i = 0; i < v; i++){
        y[i]=arc_capacity_lambda(((double)i)/((double)v),snr_sqrt_linear);
    }
    double *F0 = MALLOC(v+1,probability);
    double *F1 = MALLOC(v+1,probability);
    for (uint i=0;i<v;i++){
        F0[i] = F(y[i], snr_sqrt_linear);
        F1[i] = F(y[i],-snr_sqrt_linear);
    }
    F0[v] = 1.0;
    F1[v] = 1.0;

    probability** w = MALLOC(u,probability*);
    for(uint i=0;i<v;i++){
        w[i] = MALLOC(2,probability);
        w[i][0] = F0[i+1]-F0[i];
        w[i][1] = F1[i+1]-F1[i];
        w[u-1-i] = MALLOC(2,probability);
        w[u-1-i][0] = w[i][1];
        w[u-1-i][1] = w[i][0];
    }
    FREE(y);
    FREE(F0);
    FREE(F1);
    return w;
}

// standart normal distribution CDF
probability PHI(double x){
    return 0.5 + 0.5 * erf(x/sqrt(2.0));
}

probability F(double x, double u){
    return PHI((x-u)/stddev);
}

double capacity_lambda(double y, double snr_sqrt_linear){
    double c = 2.0*snr_sqrt_linear/pow(stddev,2);
    double b = exp(c*y);
    return 1.0-log2(b+1.0) + (c/log(2))*y/(1.0+1.0/b);
}

double arc_capacity_lambda(double C, double snr_aqrt_linear){
    if (EQUAL(C,0.0)) return 0.0;
    double a = 0.0;
    double b = 1.0;
    while(LT(capacity_lambda(b,snr_aqrt_linear),C))b*=2;
    if (!EQUAL(b,1.0)) a = 0.5 * b;
    double mid;
    while(!EQUAL(a,b)){
        mid = 0.5*(a+b);
        if(LT(capacity_lambda(mid,snr_aqrt_linear),C)) a = mid;
        else b = mid;
    }
    return b;
}

#endif // #ifndef CHANNEL_INCLUDE