#ifndef DECODING_INCLUDE
#define DECODING_INCLUDE

#include "crc.h"
#include "encoding.h"
#include "priority_queue.h"

typedef struct decoder{
    int   m;
    int   n;
    int   L;
    int  *I;
    int **swap;

    int  cur_size;
    double* slist;

    double* mem_buf;

    // distribution
    double**   ditems;
    int*      dfrindp;
    double**    dlist;
    double** oridlist;
    double** newdlist;
    // result
    int**      ritems;
    int*      rfrindp;
    int**       rlist;
    int**    orirlist;
    int**    newrlist;
    // helper bits
    int**      hitems;
    int*      hfrindp;
    int**       hlist;
    int**    orihlist;
    int**    newhlist;

    
    PriorityQueue*     pq;
    int*         codeword;
}decoder;

decoder*  dec_init(int m, int L, int* I, int** swap);
void      dec_dele(decoder* dec);
void        decode(const double* llr, int* deResult, CRC* crc, encoder* enc, decoder* dec);
void         error(const int* codeword, const int* deResult, const double* llr, int n, int* err, int* ML_err);

//=============================decoding.c==============================

#include "transform.h"
#include <string.h>

#define D_LIST_POP(        cur_m, dec)  ( (dec)->ditems[cur_m] + (1 << (cur_m + 2)) * ((dec)->dfrindp[cur_m]++))
#define ORI_D_LISTP( iter, cur_m, dec)  ( (dec)->oridlist[ (iter)*((dec)->m) + (cur_m) ] )
#define NEW_D_LISTP( iter, cur_m, dec)  ( (dec)->newdlist[ (iter)*((dec)->m) + (cur_m) ] )
#define ORI_D_LISTPP(iter, cur_m, dec)  ( (dec)->oridlist+ (iter)*((dec)->m) + (cur_m) )
#define NEW_D_LISTPP(iter, cur_m, dec)  ( (dec)->newdlist+ (iter)*((dec)->m) + (cur_m) )

#define R_LIST_POP(        cur_m, dec)  ( (dec)->ritems[cur_m] + (1 << (cur_m + 1)) * ((dec)->rfrindp[cur_m]++))
#define ORI_R_LISTP( iter, cur_m, dec)  ( (dec)->orirlist[ (iter)*((dec)->m) + (cur_m) ] )
#define NEW_R_LISTP( iter, cur_m, dec)  ( (dec)->newrlist[ (iter)*((dec)->m) + (cur_m) ] )
#define ORI_R_LISTPP(iter, cur_m, dec)  ( (dec)->orirlist+ (iter)*((dec)->m) + (cur_m) )
#define NEW_R_LISTPP(iter, cur_m, dec)  ( (dec)->newrlist+ (iter)*((dec)->m) + (cur_m) )

#define H_LIST_POP(        cur_m, dec)  ( (dec)->hitems[cur_m] + (1 << (cur_m    )) * ((dec)->hfrindp[cur_m]++))
#define ORI_H_LISTP( iter, cur_m, dec)  ( (dec)->orihlist[ (iter)*((dec)->m) + (cur_m) ] )
#define NEW_H_LISTP( iter, cur_m, dec)  ( (dec)->newhlist[ (iter)*((dec)->m) + (cur_m) ] )
#define ORI_H_LISTPP(iter, cur_m, dec)  ( (dec)->orihlist+ (iter)*((dec)->m) + (cur_m) )
#define NEW_H_LISTPP(iter, cur_m, dec)  ( (dec)->newhlist+ (iter)*((dec)->m) + (cur_m) )

#define PQ_NODE_PUSH(temp, iter, x0, x1, pq, dec)({\
        if((pq)->cur_size<(dec)->L){\
            pq_push(node_creat(temp,iter,x0,x1), pq);\
        }else if (GT(temp, pq_top(pq)->val)){\
            pq_pop(pq);\
            pq_push(node_creat(temp,iter,x0,x1), pq);\
        }\
    })

void recursive_list_decode(int cur_m, int branch, decoder* dec);
int  parity_check(const int* res, CRC* crc, encoder* enc, decoder* dec);
void res_to_code( const int* res, int* codeword, int n);
void swap_pointers(decoder* dec);

decoder* dec_init(int m, int L, int* I, int** swap){
    decoder* dec = MALLOC(1, decoder);
    dec->m = m;
    dec->n = 1 << m;
    dec->L = L;
    dec->I    = I;
    dec->swap = swap;

    dec->slist = MALLOC(L, double);

    int ds = 4, rs = 8, hs = 2;

    size_t mem = (dec->n * L)*(ds * 4*sizeof(double) + rs * 2*sizeof(int) + hs * sizeof(int));
    dec->mem_buf = (double*)malloc(mem);
    
    dec->ditems = MALLOC(m, double*);
    dec->ritems = MALLOC(m, int*);
    dec->hitems = MALLOC(m, int*);
    dec->ditems[0] = dec->mem_buf;
    dec->ritems[0] = (int*)(dec->mem_buf + dec->n * L * ds * 4);
    dec->hitems[0] =      dec->ritems[0] + dec->n * L * rs * 2;
    for(int i = 1; i< m; i++){
        dec->ditems[i] = dec->ditems[i-1] + (1<<(i-1)) * L * ds * 4;
        dec->ritems[i] = dec->ritems[i-1] + (1<<(i-1)) * L * rs * 2;
        dec->hitems[i] = dec->hitems[i-1] + (1<<(i-1)) * L * hs    ;        
    }
    dec->dfrindp = MALLOC(3 * m, int);
    dec->rfrindp = dec->dfrindp + m;
    dec->hfrindp = dec->rfrindp + m;
    
    dec->dlist = MALLOC(2 * L * m, double*);
    dec->rlist = MALLOC(2 * L * m, int*   );
    dec->hlist = MALLOC(2 * L * m, int*   );
    dec->oridlist = dec->dlist;
    dec->orirlist = dec->rlist;
    dec->orihlist = dec->hlist;
    dec->newdlist = dec->dlist + dec->m * dec->L;
    dec->newrlist = dec->rlist + dec->m * dec->L;
    dec->newhlist = dec->hlist + dec->m * dec->L;

    dec->pq = pq_init(L);
    dec->codeword = MALLOC(dec->n, int);

    return dec;
}

void dec_dele(decoder* dec){
    FREE(dec->codeword);
    pq_dele(dec->pq);

    FREE(dec->hlist);
    FREE(dec->rlist);
    FREE(dec->dlist);
    
    FREE(dec->dfrindp);

    FREE(dec->hitems);
    FREE(dec->ritems);
    FREE(dec->ditems);

    FREE(dec->mem_buf);    
    FREE(dec->slist);
    FREE(dec);
}

void decode(const double* llr, int* deResult, CRC* crc, encoder* enc, decoder* dec){
    //clear
    int len = 3 * dec->m;
    for(int i = 0; i < len; i++) dec->dfrindp[i] = 0;

    // init
    dec->slist[0] = 1.0;
    dec->cur_size = 1;
    double* dis = ORI_D_LISTP(0, dec->m-1, dec) = D_LIST_POP(dec->m-1, dec);
    int length = 1 << (dec->m-1);
    for(int i = 0; i < length; i++){
        dis_init(dis+(i<<2), LLR_P0(llr[i]), LLR_P0(llr[i+length]));
    }
    
    recursive_list_decode(dec->m-1, 0, dec);
    
    double curMax;
    int    maxidx = -1;

    for(int iter = 0, *res; iter < dec->cur_size; iter++){
        res = ORI_R_LISTP(iter, dec->m-1, dec);
        if(maxidx==-1){ if(!parity_check(res, crc, enc, dec))continue;  }
        else          { if(!GT(dec->slist[iter], curMax) || !parity_check(res, crc, enc, dec)) continue; }
        curMax = dec->slist[iter];
        maxidx = iter;
    }
    
    if(maxidx == -1){
        for(int i = 0; i < dec->n; i++) deResult[i] = 0;
    }else{
        int* res = ORI_R_LISTP(maxidx, dec->m-1, dec);
        res_to_code(res, deResult, dec->n);
    }
}

void  error(const int* codeword, const int* deResult, const double* llr, int n, int* err, int* ML_err){
    *err = 0;
    for(int i = 0; i < n; i++){
        if (deResult[i]!=codeword[i]){ *err = 1; break; }
    }

    if(*err){
        double true_llr_sum = 0.0;
        double  dec_llr_sum = 0.0;
        for(int i = 0; i < n; i++){
            true_llr_sum+=(llr[i])*(1.0-2.0*codeword[i]);
             dec_llr_sum+=(llr[i])*(1.0-2.0*deResult[i]);
        }
        *ML_err = LT(true_llr_sum, dec_llr_sum)?(1):(0);
    }else{
        *ML_err = 0;
    }
}

void recursive_list_decode(int cur_m, int branch, decoder* dec){
    if (cur_m==0){
        PriorityQueue* pq = dec->pq;
        double temp;
        int flag;
        if(branch!=dec->n-2){    // only decode the first bit.
            if(!dec->I[branch]){         // frozen bit
                flag = 0;
                double* dis;
                int*    res;
                for(int iter = 0; iter<dec->cur_size; iter++){
                    dis = ORI_D_LISTP(iter, 0, dec);
                    dec->slist[iter]+=logproba(0, p0_minus(dis));
                    res = ORI_R_LISTP(iter, 0, dec) = R_LIST_POP(0, dec);
                    res[0] = 0;
                }
            }else{                       // information bit
                flag = 1;
                double* dis;
                for(int iter = 0; iter<dec->cur_size; iter++){
                    dis = ORI_D_LISTP(iter, 0, dec);
                    temp = dec->slist[iter] + logproba(0, p0_minus(dis));
                    PQ_NODE_PUSH(temp, iter, 0, -1, pq, dec);
                    temp = dec->slist[iter] + logproba(1, p0_minus(dis));
                    PQ_NODE_PUSH(temp, iter, 1, -1, pq, dec);
                }
            }
        }else{                   // need to decode two bits
            flag = (dec->I[branch]<<1) + dec->I[branch+1];
            if(flag==0){                 // both are frozen bits
                double* dis;
                int*    res;
                for(int iter = 0; iter<dec->cur_size; iter++){
                    dis = ORI_D_LISTP(iter, 0, dec);
                    dec->slist[iter]+=log(dis[0]);
                    res = ORI_R_LISTP(iter, 0, dec) = R_LIST_POP(0, dec);
                    res[0] = res[1] = 0;
                }
            }else if (flag == 1){        //the branch-th bit is frozen bit, and the (branch+1)-th bit is information bit
                double* dis;
                for(int iter = 0; iter<dec->cur_size; iter++){
                    dis = ORI_D_LISTP(iter, 0, dec);
                    temp = dec->slist[iter] + log(dis[0]);
                    PQ_NODE_PUSH(temp, iter, 0, 0, pq, dec);
                    temp = dec->slist[iter] + log(dis[1]);
                    PQ_NODE_PUSH(temp, iter, 0, 1, pq, dec);
                }
            }else if (flag == 2){        //the branch-th bit is information bit, and the (branch+1)-th bit is frozen bit
                double* dis;
                for(int iter = 0; iter<dec->cur_size; iter++){
                    dis = ORI_D_LISTP(iter, 0, dec);
                    temp = dec->slist[iter] + log(dis[0]);
                    PQ_NODE_PUSH(temp, iter, 0, 0, pq, dec);
                    temp = dec->slist[iter] + log(dis[2]);
                    PQ_NODE_PUSH(temp, iter, 1, 0, pq, dec);
                }
            }else{                       // both are information bits
                double* dis;
                for(int iter = 0; iter<dec->cur_size; iter++){
                    dis = ORI_D_LISTP(iter, 0, dec);
                    temp = dec->slist[iter] + log(dis[0]);
                    PQ_NODE_PUSH(temp, iter, 0, 0, pq, dec);
                    temp = dec->slist[iter] + log(dis[1]);
                    PQ_NODE_PUSH(temp, iter, 0, 1, pq, dec);
                    temp = dec->slist[iter] + log(dis[2]);
                    PQ_NODE_PUSH(temp, iter, 1, 0, pq, dec);
                    temp = dec->slist[iter] + log(dis[3]);
                    PQ_NODE_PUSH(temp, iter, 1, 1, pq, dec);
                }
            }
        }
        if(flag!=0){
            dec->cur_size = pq->cur_size;
            
            for(int count = dec->cur_size-1; !pq_is_empty(pq); count--){
                heap_node* node = pq_top(pq);
                dec->slist[count] = node->val;
                int iter = node->iter;
                memcpy(NEW_D_LISTPP(count, 0, dec), ORI_D_LISTPP(iter, 0, dec), dec->m*sizeof(double*));
                memcpy(NEW_R_LISTPP(count, 0, dec), ORI_R_LISTPP(iter, 0, dec), dec->m*sizeof(int*));
                memcpy(NEW_H_LISTPP(count, 0, dec), ORI_H_LISTPP(iter, 0, dec), dec->m*sizeof(int*));
                
                int* res = NEW_R_LISTP(count, 0, dec) = R_LIST_POP(0, dec);
                res[0] = node->x0;
                res[1] = node->x1;
                
                pq_pop(pq);
            }
            swap_pointers(dec);
        }
        if((branch&1)  && dec->swap[1][branch>>1]){
            for(int iter=0; iter<dec->cur_size; iter++){
                int* helper = ORI_H_LISTP(iter, 0, dec) = H_LIST_POP(0, dec);
                int* res    = ORI_R_LISTP(iter, 0, dec);
                helper[0] = res[0];
            }
        }
        dec->dfrindp[0] = 0;
        return;
    }
    
    int number = (1<<(dec->m-cur_m)) - 1;
    int length = (1<<cur_m);
    int halflength = (length>>1);

    int exchange = dec->swap[cur_m][branch];
    int store = (branch&1 && dec->swap[cur_m+1][branch>>1])?(1):(0);

    // mode of every subbranch
    int s0 = (branch==0||!(dec->swap[cur_m][branch-1]))?(1):(-1);
    int s2 = (branch==number-1)?(2):((exchange)?(1):(0));

    double *d,  *dm, *dp;
    int    *r0, *r1, *r2, *rm, *rp;

    //w0
    if(s0==1){
        for(int iter = 0; iter < dec->cur_size; iter++){
            dm = ORI_D_LISTP(iter, cur_m  , dec);
            dp = dm + (halflength<<2);
            d  = ORI_D_LISTP(iter, cur_m-1, dec) = D_LIST_POP(cur_m-1, dec);
            if(exchange){ for(int i = 0; i < halflength; i++){ p_dis0(d+(i<<2), dm+(i<<2), dp+(i<<2)); }}
            else        { for(int i = 0; i < halflength; i++){   dis0(d+(i<<2), dm+(i<<2), dp+(i<<2)); }}
        }
        recursive_list_decode(cur_m-1, branch<<1, dec);
    }else{
        for(int iter = 0, *helper; iter < dec->cur_size; iter++){
            r0     = ORI_R_LISTP(iter, cur_m-1, dec) = R_LIST_POP(cur_m-1, dec);
            helper = ORI_H_LISTP(iter, cur_m-1, dec);
            for(int i = 0; i < halflength; i++){ r0[i<<1] = helper[i]; }
        }
        dec->hfrindp[cur_m-1] = 0;
    }
    
    for(int iter = 0; iter < dec->cur_size; iter++){
        r0 = ORI_R_LISTP(iter, cur_m  , dec) = ORI_R_LISTP(iter, cur_m-1, dec);
        d  = ORI_D_LISTP(iter, cur_m-1, dec) = D_LIST_POP(cur_m-1, dec);
        dm = ORI_D_LISTP(iter, cur_m  , dec);
        dp = dm + (halflength<<2);

        if(exchange){ for(int i = 0; i < halflength; i++){ p_dis1(d+(i<<2), dm+(i<<2), dp+(i<<2), r0[(i<<1)]); }}
        else        { for(int i = 0; i < halflength; i++){   dis1(d+(i<<2), dm+(i<<2), dp+(i<<2), r0[(i<<1)]); }}
    }

    // w1
    recursive_list_decode(cur_m-1, (branch<<1)+1, dec);
    
    if(s2==0){
        for(int iter = 0; iter < dec->cur_size; iter++){
            r0 = ORI_R_LISTP(iter, cur_m  , dec);
            r1 = ORI_R_LISTP(iter, cur_m-1, dec);
            rm = ORI_R_LISTP(iter, cur_m  , dec) = R_LIST_POP(cur_m, dec);
            rp = rm + (halflength<<1);
            for(int i = 0; i < halflength; i++){
                rm[i<<1] = r0[i<<1] ^ r1[i<<1];
                rp[i<<1] = r1[i<<1];
            }
        }
    }else{
        for(int iter = 0; iter < dec->cur_size; iter++){
            r0 = ORI_R_LISTP(iter, cur_m  , dec);
            r1 = ORI_R_LISTP(iter, cur_m-1, dec);
            rm = ORI_R_LISTP(iter, cur_m  , dec) = R_LIST_POP(cur_m, dec);
            rp = rm + (halflength<<1);
            for(int i = 0; i < halflength; i++){
                rm[i<<1] = r0[i<<1];
                rp[i<<1] = r1[i<<1];
            }
            d  = ORI_D_LISTP(iter, cur_m-1, dec) = D_LIST_POP(cur_m-1, dec);
            dm = ORI_D_LISTP(iter, cur_m  , dec);
            dp = dm + (halflength<<2);
            if(exchange){ for(int i = 0; i < halflength; i++){ p_dis2(d+(i<<2), dm+(i<<2), dp+(i<<2), r0[(i<<1)], r1[(i<<1)]); }}
            else        { for(int i = 0; i < halflength; i++){   dis2(d+(i<<2), dm+(i<<2), dp+(i<<2), r0[(i<<1)], r1[(i<<1)]); }}
        }
        // w2
        recursive_list_decode(cur_m-1, (branch<<1)+2, dec);
        if(s2==1){
            for(int iter = 0; iter < dec->cur_size; iter++){
                r0 = ORI_R_LISTP(iter, cur_m  , dec);
              //r1 = r0 + (halflength<<1);
              //r1[i<<1] has been writed in ORI_H_LISTP(iter, cur_m-1, dec)  
                r2 = ORI_R_LISTP(iter, cur_m-1, dec); // r2
                rm = ORI_R_LISTP(iter, cur_m  , dec) = R_LIST_POP(cur_m, dec);
                rp = rm + (halflength<<1);
                for(int i = 0; i < halflength; i++){
                    rm[i<<1] = r0[i<<1] ^ r2[i<<1];
                    rp[i<<1] = r2[i<<1];
                }
            }
        }else{
            for(int iter = 0; iter < dec->cur_size; iter++){
                r0 = ORI_R_LISTP(iter, cur_m  , dec);
                r1 = r0 + (halflength<<1);
                r2 = ORI_R_LISTP(iter, cur_m-1, dec);
                rm = ORI_R_LISTP(iter, cur_m  , dec) = R_LIST_POP(cur_m, dec);
                rp = rm + (halflength<<1);
                if(exchange){
                    for(int i = 0; i < halflength; i++){
                        rm[ i<<1   ] = r0[ i<<1   ] ^ r2[ i<<1   ];
                        rm[(i<<1)+1] = r1[ i<<1   ] ^ r2[(i<<1)+1];
                        rp[ i<<1   ] = r2[ i<<1   ];
                        rp[(i<<1)+1] = r2[(i<<1)+1];
                    }
                }else{
                    for(int i = 0; i < halflength; i++){
                        rm[ i<<1   ] = r0[ i<<1   ] ^ r1[ i<<1   ];
                        rm[(i<<1)+1] = r2[ i<<1   ] ^ r2[(i<<1)+1];
                        rp[ i<<1   ] = r1[ i<<1   ];
                        rp[(i<<1)+1] = r2[(i<<1)+1];
                    }
                }
            }
        }
    }

   
    if(store){
        for(int iter = 0, *helper; iter < dec->cur_size; iter++){
            helper = ORI_H_LISTP(iter, cur_m, dec) = H_LIST_POP(cur_m, dec);
                rm = ORI_R_LISTP(iter, cur_m, dec);
            for(int i = 0; i < length; i++){ helper[i] = rm[i<<1]; }
        }
    }

    dec->dfrindp[cur_m] = dec->rfrindp[cur_m-1] = 0;
}

int parity_check(const int* res, CRC* crc, encoder* enc, decoder* dec){
    if (crc==NULL) return 1;
    res_to_code(res, dec->codeword, dec->n);

    int* message_crc = dec->codeword;
    inverse_encode(message_crc, dec->codeword, enc);
    
    return crc_check(message_crc, crc);
}

void res_to_code(const int* res, int* codeword, int n){
    int length = (n >> 1);
    for(int i = 0; i < length; i++){ 
        codeword[i       ] = res[(i<<1)  ];
        codeword[i+length] = res[(i<<1)+1];
        codeword[i]^=codeword[i+length];
    }
}

void swap_pointers(decoder* dec){
    double** dlp = dec->oridlist;
             dec->oridlist = dec->newdlist;
             dec->newdlist = dlp;
    int   ** rlp = dec->orirlist;
             dec->orirlist = dec->newrlist;
             dec->newrlist = rlp;
    int   ** hlp = dec->orihlist;
             dec->orihlist = dec->newhlist;
             dec->newhlist = hlp;
}

#undef PQ_NODE_PUSH

#undef D_LIST_POP 
#undef ORI_D_LISTP
#undef NEW_D_LISTP
#undef ORI_D_LISTPP
#undef NEW_D_LISTPP

#undef R_LIST_POP
#undef ORI_R_LISTP
#undef NEW_R_LISTP
#undef ORI_R_LISTPP
#undef NEW_R_LISTPP

#undef H_LIST_POP
#undef ORI_H_LISTP
#undef NEW_H_LISTP
#undef ORI_H_LISTPP
#undef NEW_H_LISTPP

#endif // #ifndef DECODING_INCLUDE