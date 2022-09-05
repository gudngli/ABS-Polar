//==============================================================
//
// Interfaces of transforms used in decoder.
//  
// Copyright 2022 and onwards Guodong Li
// 
// Licensed under the Apache License, Version 2.0.
// 
//==============================================================

#ifndef NEW_TRANSFORM_INCLUDE
#define NEW_TRANSFORM_INCLUDE

#include "common.h"

const double pro_lower_threshold = 1e-15;
const double pro_upper_threshold = 1 - pro_lower_threshold;

const double llr_upper_threshold =  1e8;
const double llr_lower_threshold = -llr_upper_threshold;

#define pro_clip(p)({\
    (p) = GT(p, pro_upper_threshold) ? (pro_upper_threshold)\
                                     : (LT(p,pro_lower_threshold)?(pro_lower_threshold)\
                                                                 :(p));\
})

#define LLR_P0(llr) ({\
    GT(llr, llr_upper_threshold) ? ( pro_upper_threshold )\
                                 : ( LT(llr, llr_lower_threshold) ? (pro_lower_threshold)        \
                                                                  : (1.0 / (1.0 + exp(-llr))) ); \
})

#define p0_minus(d)     ((d[0])+(d[1]))
#define logproba(idx, p0) (   (idx)? (log(1-p0)) : (log(p0) ) )

// p = Pr[x0 = 0 | y0], q = Pr[x1 = 0 | y1]
void dis_init(double* d, double p, double q);
void oria(double* d, double* d0, double* d1);
void orib(double* d, double* d0, double* d1, int u0);
void oric(double* d, double* d0, double* d1, int u0, int u1);
void swpa(double* d, double* d0, double* d1);
void swpb(double* d, double* d0, double* d1, int u0);
void swpc(double* d, double* d0, double* d1, int u0, int u1);
void adda(double* d, double* d0, double* d1);
void addb(double* d, double* d0, double* d1, int u0);
void addc(double* d, double* d0, double* d1, int u0, int u1);
//================================transform.c=============================

void dis_init(double* d, double p, double q){
    d[0] =    p   *   q  ;
    d[1] =  (1-p) * (1-q);
    d[2] =  (1-p) *   q  ;
    d[3] =    p   * (1-q);
    pro_clip(d[0]);pro_clip(d[1]);pro_clip(d[2]);pro_clip(d[3]);
}

void oria(double* d, double* d0, double* d1){
    d[0] = d0[0]*d1[0] + d0[1]*d1[1] + d0[1]*d1[0] + d0[0]*d1[1];
    d[1] = d0[2]*d1[2] + d0[3]*d1[3] + d0[3]*d1[2] + d0[2]*d1[3];
    d[2] = d0[2]*d1[0] + d0[3]*d1[1] + d0[3]*d1[0] + d0[2]*d1[1];
    d[3] = d0[0]*d1[2] + d0[1]*d1[3] + d0[1]*d1[2] + d0[0]*d1[3];
    pro_clip(d[0]);pro_clip(d[1]);pro_clip(d[2]);pro_clip(d[3]);
}

void orib(double* d, double* d0, double* d1, int u0){
    if (u0){// u0 = 1
        d[0] = d0[2]*d1[0] + d0[3]*d1[1];
        d[1] = d0[3]*d1[0] + d0[2]*d1[1];
        d[2] = d0[0]*d1[2] + d0[1]*d1[3];
        d[3] = d0[1]*d1[2] + d0[0]*d1[3];
    }else{  // u0 = 0
        d[0] = d0[0]*d1[0] + d0[1]*d1[1];
        d[1] = d0[1]*d1[0] + d0[0]*d1[1];
        d[2] = d0[2]*d1[2] + d0[3]*d1[3];
        d[3] = d0[3]*d1[2] + d0[2]*d1[3];
    }
    pro_clip(d[0]);pro_clip(d[1]);pro_clip(d[2]);pro_clip(d[3]);
    double factor = d[0] + d[1] + d[2] + d[3];
    d[0]/=factor;
    d[1]/=factor;
    d[2]/=factor;
    d[3]/=factor;
}

void oric(double* d, double* d0, double* d1, int u0, int u1){
    int u = (u0<<1) + u1;
    if (u == 0){
        d[0] = d0[0]*d1[0];
        d[1] = d0[1]*d1[1];
        d[2] = d0[1]*d1[0];
        d[3] = d0[0]*d1[1];
    }else if (u==1){
        d[0] = d0[2]*d1[2];
        d[1] = d0[3]*d1[3];
        d[2] = d0[3]*d1[2];
        d[3] = d0[2]*d1[3];
    }else if (u==2){
        d[0] = d0[2]*d1[0];
        d[1] = d0[3]*d1[1];
        d[2] = d0[3]*d1[0];
        d[3] = d0[2]*d1[1];
    }else{ // u==3
        d[0] = d0[0]*d1[2];
        d[1] = d0[1]*d1[3];
        d[2] = d0[1]*d1[2];
        d[3] = d0[0]*d1[3];
    }
    pro_clip(d[0]);pro_clip(d[1]);pro_clip(d[2]);pro_clip(d[3]);
    double factor = d[0] + d[1] + d[2] + d[3];
    d[0] /= factor;
    d[1] /= factor;
    d[2] /= factor;
    d[3] /= factor;
}

void swpa(double* d, double* d0, double* d1){
    d[0] = d0[0]*d1[0] + d0[1]*d1[1] + d0[2]*d1[2] + d0[3]*d1[3];
    d[1] = d0[1]*d1[0] + d0[0]*d1[1] + d0[3]*d1[2] + d0[2]*d1[3];
    d[2] = d0[2]*d1[0] + d0[3]*d1[1] + d0[0]*d1[2] + d0[1]*d1[3];
    d[3] = d0[3]*d1[0] + d0[2]*d1[1] + d0[1]*d1[2] + d0[0]*d1[3];
    pro_clip(d[0]);pro_clip(d[1]);pro_clip(d[2]);pro_clip(d[3]);
}

void swpb(double* d, double* d0, double* d1, int u0){
    if (u0){// u0 = 1
        d[0] = d0[2]*d1[0] + d0[3]*d1[1];
        d[1] = d0[0]*d1[2] + d0[1]*d1[3];
        d[2] = d0[3]*d1[0] + d0[2]*d1[1];
        d[3] = d0[1]*d1[2] + d0[0]*d1[3];
    }else{  // u0 = 0
        d[0] = d0[0]*d1[0] + d0[1]*d1[1];
        d[1] = d0[2]*d1[2] + d0[3]*d1[3];
        d[2] = d0[1]*d1[0] + d0[0]*d1[1];
        d[3] = d0[3]*d1[2] + d0[2]*d1[3];
    }
    pro_clip(d[0]);pro_clip(d[1]);pro_clip(d[2]);pro_clip(d[3]);
    double factor = d[0] + d[1] + d[2] + d[3];
    d[0] /= factor;
    d[1] /= factor;
    d[2] /= factor;
    d[3] /= factor;
}


void swpc(double* d, double* d0, double* d1, int u0, int u1){
    int u = (u0<<1) + u1;
    if (u == 0){
        d[0] = d0[0]*d1[0];
        d[1] = d0[1]*d1[1];
        d[2] = d0[2]*d1[2];
        d[3] = d0[3]*d1[3];
    }else if (u==1){
        d[0] = d0[1]*d1[0];
        d[1] = d0[0]*d1[1];
        d[2] = d0[3]*d1[2];
        d[3] = d0[2]*d1[3];
    }else if (u==2){
        d[0] = d0[2]*d1[0];
        d[1] = d0[3]*d1[1];
        d[2] = d0[0]*d1[2];
        d[3] = d0[1]*d1[3];
    }else{ // u==3
        d[0] = d0[3]*d1[0];
        d[1] = d0[2]*d1[1];
        d[2] = d0[1]*d1[2];
        d[3] = d0[0]*d1[3];
    }
    pro_clip(d[0]);pro_clip(d[1]);pro_clip(d[2]);pro_clip(d[3]);
    double factor = d[0] + d[1] + d[2] + d[3];
    d[0] /= factor;
    d[1] /= factor;
    d[2] /= factor;
    d[3] /= factor;
}


void adda(double* d, double* d0, double* d1){
    d[0] = d0[0]*d1[0] + d0[1]*d1[1] + d0[3]*d1[2] + d0[2]*d1[3];
    d[1] = d0[2]*d1[2] + d0[3]*d1[3] + d0[1]*d1[0] + d0[0]*d1[1];
    d[2] = d0[2]*d1[0] + d0[3]*d1[1] + d0[1]*d1[2] + d0[0]*d1[3];
    d[3] = d0[0]*d1[2] + d0[1]*d1[3] + d0[3]*d1[0] + d0[2]*d1[1];
    pro_clip(d[0]);pro_clip(d[1]);pro_clip(d[2]);pro_clip(d[3]);
}

void addb(double* d, double* d0, double* d1, int u0){
    if (u0){// u0 = 1
        d[0] = d0[2]*d1[0] + d0[3]*d1[1];
        d[1] = d0[1]*d1[2] + d0[0]*d1[3];
        d[2] = d0[0]*d1[2] + d0[1]*d1[3];
        d[3] = d0[3]*d1[0] + d0[2]*d1[1];
    }else{  // u0 = 0
        d[0] = d0[0]*d1[0] + d0[1]*d1[1];
        d[1] = d0[3]*d1[2] + d0[2]*d1[3];
        d[2] = d0[2]*d1[2] + d0[3]*d1[3];
        d[3] = d0[1]*d1[0] + d0[0]*d1[1];
    }
    pro_clip(d[0]);pro_clip(d[1]);pro_clip(d[2]);pro_clip(d[3]);
    double factor = d[0] + d[1] + d[2] + d[3];
    d[0] /= factor;
    d[1] /= factor;
    d[2] /= factor;
    d[3] /= factor;
}

void addc(double* d, double* d0, double* d1, int u0, int u1){
    int u = (u0<<1) + u1;
    if (u == 0){
        d[0] = d0[0]*d1[0];
        d[1] = d0[1]*d1[1];
        d[2] = d0[3]*d1[2];
        d[3] = d0[2]*d1[3];
    }else if (u==1){
        d[0] = d0[2]*d1[2];
        d[1] = d0[3]*d1[3];
        d[2] = d0[1]*d1[0];
        d[3] = d0[0]*d1[1];
    }else if (u==2){
        d[0] = d0[2]*d1[0];
        d[1] = d0[3]*d1[1];
        d[2] = d0[1]*d1[2];
        d[3] = d0[0]*d1[3];
    }else{ // u==3
        d[0] = d0[0]*d1[2];
        d[1] = d0[1]*d1[3];
        d[2] = d0[3]*d1[0];
        d[3] = d0[2]*d1[1];
    }
    pro_clip(d[0]);pro_clip(d[1]);pro_clip(d[2]);pro_clip(d[3]);
    double factor = d[0] + d[1] + d[2] + d[3];
    d[0] /= factor;
    d[1] /= factor;
    d[2] /= factor;
    d[3] /= factor;
}

#endif // #ifndef NEW_TRANSFORM_INCLUDE