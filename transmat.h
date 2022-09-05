//==============================================================
// 
// Interfaces of Double-Bit transition matrix used
// in code construction
//  
// Copyright 2022 and onwards Guodong Li
// 
// Licensed under the Apache License, Version 2.0.
// 
//==============================================================


#ifndef TRANSMAT_INCLUDE
#define TRANSMAT_INCLUDE

#include "interval_merge.h"

// Double bit transition matrix.
typedef struct transmat{
    int   size;     // size of output alphabat
    double** w;     // size x 4 transition matrix
                    // (00,01,10,11) --> (0,1,2,3)
}transmat;

// w is a size x 2 transition matrix of real channel.
transmat* trmat_init(double** w, int size, merge_map* map);
void      trmat_dele(transmat* trmat);

transmat*   ori_trans0(transmat* trmat, merge_map* map);
transmat*   ori_trans1(transmat* trmat, merge_map* map);
transmat*   ori_trans2(transmat* trmat, merge_map* map);
transmat*   swp_trans0(transmat* trmat, merge_map* map);
transmat*   swp_trans1(transmat* trmat, merge_map* map);
transmat*   swp_trans2(transmat* trmat, merge_map* map);
transmat*   add_trans0(transmat* trmat, merge_map* map);
transmat*   add_trans1(transmat* trmat, merge_map* map);
transmat*   add_trans2(transmat* trmat, merge_map* map);
double        profit(transmat* help, transmat* p_help);
double       gamma(transmat* trmat);
double           conditional_entropy(transmat* trmat);
double     minus_conditional_entropy(transmat* trmat);
double      plus_conditional_entropy(transmat* trmat);

//===========================transmat.c====================================

transmat* trmat_init(double** w, int size, merge_map* map){
    map_clean(map);
    double* y = MALLOC(4, double);

    for(int y0 = 0; y0 < size; y0++){
        for(int y1 = 0; y1 < size; y1++){
            y[0] = w[y0][0]*w[y1][0];
            y[1] = w[y0][1]*w[y1][1];
            y[2] = w[y0][1]*w[y1][0];
            y[3] = w[y0][0]*w[y1][1];
            map_push(y, map);
        }
    }

    transmat* trmat = MALLOC(1, transmat);
    trmat->size = map->cur_size;
    trmat->w = merge(map);
    

#ifdef  INCREASE_RATE
    increase_rate = (increase_rate)/(conditional_entropy(trmat) - increase_rate);
#endif

    return trmat;
}

void trmat_dele(transmat* trmat){
    if (!trmat) return;
    int      size = trmat->size;
    double**    w = trmat->w;
    for (int i = 0 ; i < size; i++) FREE(w[i]);
    FREE(w);
    FREE(trmat);
}


transmat*   ori_trans0(transmat* trmat, merge_map* map){
    int      size = trmat->size;
    double**    w = trmat->w;

    map_clean(map);
    double* y = MALLOC(4, double);

    for(int y0 = 0; y0 < size; y0++){
        for(int y1 = 0; y1 < size; y1++){
            y[0] = 0.25*(w[y0][0]*w[y1][0] + w[y0][1]*w[y1][1] + w[y0][1]*w[y1][0] + w[y0][0]*w[y1][1]);
            y[1] = 0.25*(w[y0][2]*w[y1][2] + w[y0][3]*w[y1][3] + w[y0][3]*w[y1][2] + w[y0][2]*w[y1][3]);
            y[2] = 0.25*(w[y0][2]*w[y1][0] + w[y0][3]*w[y1][1] + w[y0][3]*w[y1][0] + w[y0][2]*w[y1][1]);
            y[3] = 0.25*(w[y0][0]*w[y1][2] + w[y0][1]*w[y1][3] + w[y0][1]*w[y1][2] + w[y0][0]*w[y1][3]);
            map_push(y, map);
        }
    }

    transmat* ori0 = MALLOC(1, transmat);
    ori0->size = map->cur_size;
    ori0->w = merge(map);
    

#ifdef  INCREASE_RATE
    increase_rate = (increase_rate)/(conditional_entropy(trmat) - increase_rate);
#endif

    return ori0;
}

transmat*   ori_trans1(transmat* trmat, merge_map* map){
    int      size = trmat->size;
    double**    w = trmat->w;

    map_clean(map);
    double* y = MALLOC(4, double);

    for(int y0 = 0; y0 < size; y0++){
        for(int y1 = 0; y1 < size; y1++){
            //(u1, u2) --> (0, y1, y2)
            y[0] = 0.25*(w[y0][0]*w[y1][0] + w[y0][1]*w[y1][1]);
            y[1] = 0.25*(w[y0][1]*w[y1][0] + w[y0][0]*w[y1][1]);
            y[2] = 0.25*(w[y0][2]*w[y1][2] + w[y0][3]*w[y1][3]);
            y[3] = 0.25*(w[y0][3]*w[y1][2] + w[y0][2]*w[y1][3]);
            map_push(y, map);

            //(u1, u2) --> (1, y1, y2)
            y[0] = 0.25*(w[y0][2]*w[y1][0] + w[y0][3]*w[y1][1]);
            y[1] = 0.25*(w[y0][3]*w[y1][0] + w[y0][2]*w[y1][1]);
            y[2] = 0.25*(w[y0][0]*w[y1][2] + w[y0][1]*w[y1][3]);
            y[3] = 0.25*(w[y0][1]*w[y1][2] + w[y0][0]*w[y1][3]);
            map_push(y, map);
        }
    }

    transmat* ori1 = MALLOC(1, transmat);
    ori1->size = map->cur_size;
    ori1->w = merge(map);

#ifdef  INCREASE_RATE
    increase_rate = (increase_rate)/(conditional_entropy(trmat) - increase_rate);
#endif

    return ori1;
}

transmat*   ori_trans2(transmat* trmat, merge_map* map){
    int      size = trmat->size;
    double**    w = trmat->w;

    map_clean(map);
    double* y = MALLOC(4, double);

    for(int y0 = 0; y0 < size; y0++){
        for(int y1 = 0; y1 < size; y1++){
            //(u2, u3) --> (0, 0, y1, y2)
            y[0] = 0.25*w[y0][0]*w[y1][0];
            y[1] = 0.25*w[y0][1]*w[y1][1];
            y[2] = 0.25*w[y0][1]*w[y1][0];
            y[3] = 0.25*w[y0][0]*w[y1][1];
            map_push(y, map);

            //(u2, u3) --> (0, 1, y1, y2)
            y[0] = 0.25*w[y0][2]*w[y1][2];
            y[1] = 0.25*w[y0][3]*w[y1][3];
            y[2] = 0.25*w[y0][3]*w[y1][2];
            y[3] = 0.25*w[y0][2]*w[y1][3];
            map_push(y, map);

            //(u2, u3) --> (1, 0, y1, y2)
            y[0] = 0.25*w[y0][2]*w[y1][0];
            y[1] = 0.25*w[y0][3]*w[y1][1];
            y[2] = 0.25*w[y0][3]*w[y1][0];
            y[3] = 0.25*w[y0][2]*w[y1][1];
            map_push(y, map);

            //(u2, u3) --> (1, 1, y1, y2)
            y[0] = 0.25*w[y0][0]*w[y1][2];
            y[1] = 0.25*w[y0][1]*w[y1][3];
            y[2] = 0.25*w[y0][1]*w[y1][2];
            y[3] = 0.25*w[y0][0]*w[y1][3];
            map_push(y, map);
        }
    }

    transmat* ori2 = MALLOC(1, transmat);
    ori2->size = map->cur_size;
    ori2->w = merge(map);

#ifdef  INCREASE_RATE
    increase_rate = (increase_rate)/(conditional_entropy(trmat) - increase_rate);
#endif

    return ori2;
}


transmat*   swp_trans0(transmat* trmat, merge_map* map){
    int      size = trmat->size;
    double**    w = trmat->w;

    map_clean(map);
    double* y = MALLOC(4, double);

    for(int y0 = 0; y0 < size; y0++){
        for(int y1 = 0; y1 < size; y1++){
            y[0] = 0.25*(w[y0][0]*w[y1][0] + w[y0][1]*w[y1][1] + w[y0][2]*w[y1][2] + w[y0][3]*w[y1][3]);
            y[1] = 0.25*(w[y0][1]*w[y1][0] + w[y0][0]*w[y1][1] + w[y0][3]*w[y1][2] + w[y0][2]*w[y1][3]);
            y[2] = 0.25*(w[y0][2]*w[y1][0] + w[y0][3]*w[y1][1] + w[y0][0]*w[y1][2] + w[y0][1]*w[y1][3]);
            y[3] = 0.25*(w[y0][3]*w[y1][0] + w[y0][2]*w[y1][1] + w[y0][1]*w[y1][2] + w[y0][0]*w[y1][3]);
            map_push(y, map);
        }
    }

    transmat* swp0 = MALLOC(1, transmat);
    swp0->size = map->cur_size;
    swp0->w = merge(map);

#ifdef  INCREASE_RATE
    increase_rate = (increase_rate)/(conditional_entropy(trmat) - increase_rate);
#endif

    return swp0;
}


transmat*   swp_trans1(transmat* trmat, merge_map* map){
    int      size = trmat->size;
    double**    w = trmat->w;

    map_clean(map);
    double* y = MALLOC(4, double);

    for(int y0 = 0; y0 < size; y0++){
        for(int y1 = 0; y1 < size; y1++){
            //(u1, u2) --> (0, y1, y2)
            y[0] = 0.25*(w[y0][0]*w[y1][0] + w[y0][1]*w[y1][1]);
            y[1] = 0.25*(w[y0][2]*w[y1][2] + w[y0][3]*w[y1][3]);
            y[2] = 0.25*(w[y0][1]*w[y1][0] + w[y0][0]*w[y1][1]);
            y[3] = 0.25*(w[y0][3]*w[y1][2] + w[y0][2]*w[y1][3]);
            map_push(y, map);

            //(u1, u2) --> (1, y1, y2)
            y[0] = 0.25*(w[y0][2]*w[y1][0] + w[y0][3]*w[y1][1]);
            y[1] = 0.25*(w[y0][0]*w[y1][2] + w[y0][1]*w[y1][3]);
            y[2] = 0.25*(w[y0][3]*w[y1][0] + w[y0][2]*w[y1][1]);
            y[3] = 0.25*(w[y0][1]*w[y1][2] + w[y0][0]*w[y1][3]);
            map_push(y, map);
        }
    }

    transmat* swp1 = MALLOC(1, transmat);
    swp1->size = map->cur_size;
    swp1->w = merge(map);

#ifdef  INCREASE_RATE
    increase_rate = (increase_rate)/(conditional_entropy(trmat) - increase_rate);
#endif

    return swp1;
}

;
transmat*   swp_trans2(transmat* trmat, merge_map* map){
    int      size = trmat->size;
    double**    w = trmat->w;

    map_clean(map);
    double* y = MALLOC(4, double);

    for(int y0 = 0; y0 < size; y0++){
        for(int y1 = 0; y1 < size; y1++){
            //(u2, u3) --> (0, 0, y1, y2)
            y[0] = 0.25*w[y0][0]*w[y1][0];
            y[1] = 0.25*w[y0][1]*w[y1][1];
            y[2] = 0.25*w[y0][2]*w[y1][2];
            y[3] = 0.25*w[y0][3]*w[y1][3];
            map_push(y, map);

            //(u2, u3) --> (0, 1, y1, y2)
            y[0] = 0.25*w[y0][1]*w[y1][0];
            y[1] = 0.25*w[y0][0]*w[y1][1];
            y[2] = 0.25*w[y0][3]*w[y1][2];
            y[3] = 0.25*w[y0][2]*w[y1][3];
            map_push(y, map);

            //(u2, u3) --> (1, 0, y1, y2)
            y[0] = 0.25*w[y0][2]*w[y1][0];
            y[1] = 0.25*w[y0][3]*w[y1][1];
            y[2] = 0.25*w[y0][0]*w[y1][2];
            y[3] = 0.25*w[y0][1]*w[y1][3];
            map_push(y, map);

            //(u2, u3) --> (1, 1, y1, y2)
            y[0] = 0.25*w[y0][3]*w[y1][0];
            y[1] = 0.25*w[y0][2]*w[y1][1];
            y[2] = 0.25*w[y0][1]*w[y1][2];
            y[3] = 0.25*w[y0][0]*w[y1][3];
            map_push(y, map);
        }
    }

    transmat* swp2 = MALLOC(1, transmat);
    swp2->size = map->cur_size;
    swp2->w = merge(map);

#ifdef  INCREASE_RATE
    increase_rate = (increase_rate)/(conditional_entropy(trmat) - increase_rate);
#endif

    return swp2;
}

;
transmat*   add_trans0(transmat* trmat, merge_map* map){
    int      size = trmat->size;
    double**    w = trmat->w;

    map_clean(map);
    double* y = MALLOC(4, double);

    for(int y0 = 0; y0 < size; y0++){
        for(int y1 = 0; y1 < size; y1++){
            y[0] = 0.25*(w[y0][0]*w[y1][0] + w[y0][1]*w[y1][1] + w[y0][3]*w[y1][2] + w[y0][2]*w[y1][3]);
            y[1] = 0.25*(w[y0][2]*w[y1][2] + w[y0][3]*w[y1][3] + w[y0][1]*w[y1][0] + w[y0][0]*w[y1][1]);
            y[2] = 0.25*(w[y0][2]*w[y1][0] + w[y0][3]*w[y1][1] + w[y0][1]*w[y1][2] + w[y0][0]*w[y1][3]);
            y[3] = 0.25*(w[y0][0]*w[y1][2] + w[y0][1]*w[y1][3] + w[y0][3]*w[y1][0] + w[y0][2]*w[y1][1]);
            map_push(y, map);
        }
    }

    transmat* add0 = MALLOC(1, transmat);
    add0->size = map->cur_size;
    add0->w = merge(map);

#ifdef  INCREASE_RATE
    increase_rate = (increase_rate)/(conditional_entropy(trmat) - increase_rate);
#endif

    return add0;
}

transmat*   add_trans1(transmat* trmat, merge_map* map){
    int      size = trmat->size;
    double**    w = trmat->w;

    map_clean(map);
    double* y = MALLOC(4, double);

    for(int y0 = 0; y0 < size; y0++){
        for(int y1 = 0; y1 < size; y1++){
            //(u1, u2) --> (0, y1, y2)
            y[0] = 0.25*(w[y0][0]*w[y1][0] + w[y0][1]*w[y1][1]);
            y[1] = 0.25*(w[y0][3]*w[y1][2] + w[y0][2]*w[y1][3]);
            y[2] = 0.25*(w[y0][2]*w[y1][2] + w[y0][3]*w[y1][3]);
            y[3] = 0.25*(w[y0][1]*w[y1][0] + w[y0][0]*w[y1][1]);
            map_push(y, map);

            //(u1, u2) --> (1, y1, y2)
            y[0] = 0.25*(w[y0][2]*w[y1][0] + w[y0][3]*w[y1][1]);
            y[1] = 0.25*(w[y0][1]*w[y1][2] + w[y0][0]*w[y1][3]);
            y[2] = 0.25*(w[y0][0]*w[y1][2] + w[y0][1]*w[y1][3]);
            y[3] = 0.25*(w[y0][3]*w[y1][0] + w[y0][2]*w[y1][1]);
            map_push(y, map);
        }
    }

    transmat* add1 = MALLOC(1, transmat);
    add1->size = map->cur_size;
    add1->w = merge(map);

#ifdef  INCREASE_RATE
    increase_rate = (increase_rate)/(conditional_entropy(trmat) - increase_rate);
#endif

    return add1;
}

transmat*   add_trans2(transmat* trmat, merge_map* map){
    int      size = trmat->size;
    double**    w = trmat->w;

    map_clean(map);
    double* y = MALLOC(4, double);

    for(int y0 = 0; y0 < size; y0++){
        for(int y1 = 0; y1 < size; y1++){
            //(u2, u3) --> (0, 0, y1, y2)
            y[0] = 0.25*w[y0][0]*w[y1][0];
            y[1] = 0.25*w[y0][1]*w[y1][1];
            y[2] = 0.25*w[y0][3]*w[y1][2];
            y[3] = 0.25*w[y0][2]*w[y1][3];
            map_push(y, map);

            //(u2, u3) --> (0, 1, y1, y2)
            y[0] = 0.25*w[y0][2]*w[y1][2];
            y[1] = 0.25*w[y0][3]*w[y1][3];
            y[2] = 0.25*w[y0][1]*w[y1][0];
            y[3] = 0.25*w[y0][0]*w[y1][1];
            map_push(y, map);

            //(u2, u3) --> (1, 0, y1, y2)
            y[0] = 0.25*w[y0][2]*w[y1][0];
            y[1] = 0.25*w[y0][3]*w[y1][1];
            y[2] = 0.25*w[y0][1]*w[y1][2];
            y[3] = 0.25*w[y0][0]*w[y1][3];
            map_push(y, map);

            //(u2, u3) --> (1, 1, y1, y2)
            y[0] = 0.25*w[y0][0]*w[y1][2];
            y[1] = 0.25*w[y0][1]*w[y1][3];
            y[2] = 0.25*w[y0][3]*w[y1][0];
            y[3] = 0.25*w[y0][2]*w[y1][1];
            map_push(y, map);
        }
    }

    transmat* add2 = MALLOC(1, transmat);
    add2->size = map->cur_size;
    add2->w = merge(map);

#ifdef  INCREASE_RATE
    increase_rate = (increase_rate)/(conditional_entropy(trmat) - increase_rate);
#endif

    return add2;
}

double profit(transmat* help, transmat* p_help){
    double m0 = minus_conditional_entropy(help);
    double p0 =  plus_conditional_entropy(help);
    if (!LT(m0, p0)){
        return -10.0;
    }
    
    double m1 =minus_conditional_entropy(p_help);
    double p1 = plus_conditional_entropy(p_help);
    
    return (m0)*(1-m0)+(p0)*(1-p0) - (m1)*(1-m1) - (p1)*(1-p1);
}

double conditional_entropy(transmat* trmat){
    double** w = trmat->w;
    int size = trmat->size;
    double entropy = 0.0;

    for (int y = 0; y < size; y++){
        entropy+=singlecal4(w[y][0], w[y][1], w[y][2], w[y][3]);
    }
    return -0.25 * entropy;
}

double minus_conditional_entropy(transmat* trmat){
    double result = 0.0;
    double** w = trmat->w;
    
    for (int y = 0; y < trmat->size; y++){
        result+=singlecal2(w[y][0] + w[y][1], w[y][2] + w[y][3]);
    }
    return -0.25 * result;
}

double plus_conditional_entropy(transmat* trmat){
    double result = 0.0;
    double** w = trmat->w;
    
    for (int y = 0; y < trmat->size; y++){
        result+=singlecal2(w[y][0], w[y][1]);
        result+=singlecal2(w[y][2], w[y][3]);
    }
    return -0.25 * result;
}

double gamma(transmat* trmat){
    double minus_ce = minus_conditional_entropy(trmat);
    double  plus_ce =  plus_conditional_entropy(trmat);
    return minus_ce * (1-minus_ce) + plus_ce * (1-plus_ce);
}

#endif // #ifndef TRANSMAT_INCLUDE