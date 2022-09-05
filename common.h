//==============================================================
//
// Copyright 2022 and onwards Guodong Li
// 
// Licensed under the Apache License, Version 2.0.
// 
//==============================================================

#ifndef COMMON_INCLUDE
#define COMMON_INCLUDE

typedef double probability;
typedef unsigned uint;


//memory allocation
#include <stdlib.h>
    #define MALLOC(len, type) ((type *)(malloc((len) * sizeof(type))))
    #define FREE(pointer) ( free((void*)pointer) )


// some float operations
#include <math.h>
    #define LOG2(n) ((int)log2((double)n))

    #define MAX(a, b) (((a) > (b)) ? (a) : (b))
    #define MIN(a, b) (((a) < (b)) ? (a) : (b))
    #define  LT(a, b)  ((a) < (b))   // a<b
    #define  GT(a, b)  ((a) > (b))   // a>b

    #define NEAR_ZERO  1e-15
    #define ABS(x) (((x) < 0.0) ? (-(x)) : (x))
    #define EQUAL(x, y) (ABS((x) - (y)) < NEAR_ZERO * ABS(x))

// debug
#define __ASSERT__
#ifdef __ASSERT__
    #include <assert.h>
    #define ASSERT(value) assert(value)
#else
    #define ASSERT(value)
#endif



#endif // #ifndef COMMON_INCLUDE

