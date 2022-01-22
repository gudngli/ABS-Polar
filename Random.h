//==============================================================
//  
// Copyright 2022 and onwards Guodong Li
// 
// Licensed under the Apache License, Version 2.0.
// 
//==============================================================

#ifndef DISTRIBUTION_INCLUDE
#define DISTRIBUTION_INCLUDE

#include <random>

std::random_device rd;
std::mt19937       rng{ rd() };
std::normal_distribution<double>       std_Gauss_d(0.0f, 1.0f); 
std::uniform_real_distribution<double> urd;
#define Standard_Gauss_Distribution() (std_Gauss_d(rng))

//p in [0.0, 1.0] 
#define Binary_Distribution(p) ((urd(rng)<p)?(1):(0))
#define Gauss_Distribution(mean, stddev) (  (mean) + (stddev) * Standard_Gauss_Distribution()  )

void get_message(int* message, int k);    

//==================================Random.h========================================

void get_message(int* message, int k){
    for (int i = 0; i < k; i++){
        message[i] = Binary_Distribution(0.5);
    }
}

#endif