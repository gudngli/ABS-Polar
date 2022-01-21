#ifndef PRIORITY_QUEUE_INCLUDE
#define PRIORITY_QUEUE_INCLUDE

#include "common.h"

typedef struct node{ double val; int   iter; int x0; int x1; }heap_node;

heap_node* node_creat(double val, int iter, int idx, int sign);
int node_cmp(heap_node* a, heap_node* b);

typedef struct priority_quene{
    int capacity;
    int cur_size;
    heap_node** arr;
}PriorityQueue;

PriorityQueue* pq_init(int capacity);
int   pq_is_empty(PriorityQueue* pq);
int   pq_is_full( PriorityQueue* pq);
heap_node* pq_top(PriorityQueue* pq);
void  pq_push(heap_node* node, PriorityQueue*pq);
void  pq_pop(   PriorityQueue* pq);
void  pq_dele(  PriorityQueue* pq);
void  pq_update(PriorityQueue* pq, int index);

//====================================priority_queue.c=========================

heap_node* node_creat(double val, int iter, int x0, int x1){
    heap_node* nd = MALLOC(1,heap_node);
    nd->val  = val;
    nd->iter = iter;
    nd->x0  = x0;
    nd->x1 = x1;
    return nd;
}

PriorityQueue* pq_init(int capacity){
    PriorityQueue* pq = MALLOC(1,PriorityQueue);
    pq->capacity = capacity;
    pq->cur_size = 0;
    pq->arr = MALLOC(capacity+1,heap_node*);
    pq->arr[0]=NULL;
    return pq;
}

int pq_is_empty(PriorityQueue* pq){
    return pq->cur_size==0;
}

int pq_is_full(PriorityQueue* pq){
    return pq->cur_size==pq->capacity;
}

heap_node* pq_top(PriorityQueue* pq){
    if(pq_is_empty(pq)) return NULL;
    return pq->arr[1];
}

int  node_cmp(heap_node* a, heap_node* b){
    if (LT(a->val,b->val)) return 1;
    else return 0;
}

void pq_update(PriorityQueue* pq, int index){
    int i = index;
    heap_node* node = pq->arr[i];
    heap_node** arr = pq->arr;
    int size = pq->cur_size;

    while (i!=1 && node_cmp(node, arr[i>>1])){
        arr[i] = arr[i>>1];
        i>>=1;
    }

    int smaller_son_index;
    while ( ((i<<1)<=size &&  node_cmp(arr[i<<1], node)) ||  (((i<<1)+1)<=size && node_cmp(arr[(i<<1)+1],node))){
        if(((i<<1)+1)<=size && node_cmp(arr[(i<<1)+1],arr[i<<1])) smaller_son_index = (i<<1)+1;
        else smaller_son_index = i<<1;
        arr[i] = arr[smaller_son_index];
        i = smaller_son_index;
    }

    arr[i] = node;
}


void  pq_push(heap_node* node,PriorityQueue*pq){
    pq->cur_size+=1;
    pq->arr[pq->cur_size]=node;
    pq_update(pq,pq->cur_size);
}

void pq_pop(PriorityQueue* pq){
    FREE(pq->arr[1]);

    pq->arr[1] = pq->arr[pq->cur_size];
    pq->arr[pq->cur_size]=NULL;
    pq->cur_size--;
    pq_update(pq,1);
}

void pq_dele(PriorityQueue* pq){
    for(int i=0;i<pq->cur_size;i++)FREE(pq->arr[i]);
    FREE(pq->arr);
    FREE(pq);
}

#endif // #ifndef PRIORITY_QUEUE_INCLUDE