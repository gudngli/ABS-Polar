//==============================================================
//
// Quantization of the output alphabet.
//  
// Copyright 2022 and onwards Guodong Li
// 
// Licensed under the Apache License, Version 2.0.
// 
//==============================================================

#ifndef INTERVAL_MERGE_INCLUDE
#define INTERVAL_MERGE_INCLUDE

#include "common.h"

double singlecal2(double s0, double s1);
double singlecal4(double s0, double s1, double s2, double s3);

//#define INCREASE_RATE
#ifdef  INCREASE_RATE
    double increase;
    double increase_rate;
    #define xlogx(x) ( (LT(x, threshold))? 0.0 : (x)*log2(x) )
    double calDeltaH(const double* y0, const double* y1);
    #define INFO ({printf("increase rate = %.8f%%\n", increase_rate*100.0);})    
#else
    #define INFO 
#endif


typedef struct stack_node{ int index; struct stack_node* next; }stack_node;

stack_node* snode_init(int index);
void        snode_dele(stack_node* snode);
void        snode_push(stack_node* root, int index);
int         snode_pop( stack_node* root);

typedef struct merge_map{
    int d; // number of invervals of every dimension;
    int u; // number of grids
    double **buf;
    stack_node* root;
    int cur_size;
}merge_map;

merge_map* map_init(int d);
void map_dele( merge_map* map);
void map_clean(merge_map* map);
void map_push(const double* y, merge_map* map);
double** merge( merge_map* map);

//================================interval_merge.c=======================================

const double threshold = 1e-40;

#define INDEX(i, j, k, d) (  (d)*((d)*(i) + (j)) + (k)  )

double singlecal2(double s0, double s1){
    double s  =  s0 + s1;
    double result = 0.0;
        result += (LT(s0, threshold)?(0.0):(s0*log2(s0/s)));
        result += (LT(s1, threshold)?(0.0):(s1*log2(s1/s)));
    return  result; 
}

double singlecal4(double s0, double s1, double s2, double s3){
    double s  =  s0 + s1 + s2 + s3;
    double result = 0.0;
        result += (LT(s0, threshold)?(0.0):(s0*log2(s0/s)));
        result += (LT(s1, threshold)?(0.0):(s1*log2(s1/s)));
        result += (LT(s2, threshold)?(0.0):(s2*log2(s2/s)));
        result += (LT(s3, threshold)?(0.0):(s3*log2(s3/s)));
    return  result; 
}

merge_map* map_init(int d){
    merge_map* map = MALLOC(1, merge_map);
    map->d = d;
    map->u = (d+1)*(d+1)*(d+1);
    
    int u = map->u;
    map->buf = MALLOC(u, double*);
    double** buf = map->buf;
    // clean
    for(int i = 0; i < u; buf[i++] = NULL);
    map->root = snode_init(-1);
    map->cur_size = 0;
#ifdef  INCREASE_RATE
    increase = 0.0;
#endif
    return map;
}

void map_dele(merge_map* map){
    map_clean(map);
    FREE(map->root);
    FREE(map->buf);
    FREE(map);
}

void map_clean(merge_map* map){
    if(map->cur_size){
        double    **buf  = map->buf;
        stack_node *root = map->root;
        for(int size = map->cur_size;size;size--){
            buf[snode_pop(root)] = NULL;
        }
        map->cur_size = 0;
#ifdef  INCREASE_RATE
        increase = 0.0;
#endif
    }
}

void map_push(const double* y, merge_map* map){
    int d = map->d;
    //coefficient
    double c = y[0]+y[1]+y[2]+y[3];
    if (LT(c, 4*threshold)) return;

    double temp = (double)d;
    int index = INDEX((int)(y[0]*temp/c), (int)(y[1]*temp/c), (int)(y[2]*temp/c), d);
    
    double* p = map->buf[index];
    if(p){
#ifdef  INCREASE_RATE
        increase += calDeltaH(p, y);
#endif
        p[0]+=y[0];
        p[1]+=y[1];
        p[2]+=y[2];
        p[3]+=y[3];
    }else{
        snode_push(map->root, index);
        p = map->buf[index] =  MALLOC(4, double);
        ASSERT(p);
        p[0]=y[0];
        p[1]=y[1];
        p[2]=y[2];
        p[3]=y[3];
        map->cur_size++;
    }
}

double** merge(merge_map* map){
    double    **buf  = map->buf;
    stack_node *root = map->root;
    int         size = map->cur_size;
    int index;

    double** w = MALLOC(map->cur_size, double*);
    for(int j = 0; j < size; j++){
        index = snode_pop(root);
        w[j] = buf[index];
        buf[index] = NULL;
    }
    map->cur_size = 0;
#ifdef  INCREASE_RATE
        increase_rate = increase;
        increase = 0.0;
#endif
    return w;
}

stack_node* snode_init(int index){
    stack_node* snode = MALLOC(1, stack_node);
    snode->index = index;
    snode->next  = NULL;
    return snode;
}

void snode_dele(stack_node* snode){
    FREE(snode);
}

void snode_push(stack_node* root, int index){
    stack_node* snode = snode_init(index);
    snode->next = root->next;
    root->next = snode;
}

int snode_pop(stack_node* root){
    if (!root || !(root->next))return -1;
    stack_node* snode = root->next;
    int result = snode->index;
    root->next = snode->next;
    FREE(snode);
    return result;
}


#ifdef INCREASE_RATE
double calDeltaH(const double* y0, const double* y1){
    double* y = MALLOC(4, double);
    y[0] = y0[0] + y1[0];
    y[1] = y0[1] + y1[1];
    y[2] = y0[2] + y1[2];
    y[3] = y0[3] + y1[3];

    double s0 = y0[0]+y0[1]+y0[2]+y0[3];
    double s1 = y1[0]+y1[1]+y1[2]+y1[3];
    double s  = y[0]+y[1]+y[2]+y[3];

    double result = 0.0;

    result+=(xlogx(y[0])+xlogx(y[1])+xlogx(y[2])+xlogx(y[3])-xlogx(s));
    result-=(xlogx(y0[0])+xlogx(y0[1])+xlogx(y0[2])+xlogx(y0[3])-xlogx(s0));
    result-=(xlogx(y1[0])+xlogx(y1[1])+xlogx(y1[2])+xlogx(y1[3])-xlogx(s1));

    FREE(y);
    return -0.25 * result;
}
#endif

#endif // #ifndef INTERVAL_MERGE_INCLUDE