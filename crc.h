#ifndef CRC_INCLUDE
#define CRC_INCLUDE

#include "common.h"
#include "Random.h"

typedef struct CRC{
    int k;
    int c;
    int** H; // k x c parity check matrix
}CRC;

CRC* crc_init(int k, int c);
void crc_dele(CRC* crc);
void   set_crc(      int* message_crc, CRC* crc);
int  crc_check(const int* message_crc, CRC* crc);

//=================================crc.c=======================

CRC* crc_init(int k, int c){
    if (k<=0||c<=0) return NULL;
    CRC* crc = MALLOC(1, CRC);
    crc->k = k;
    crc->c = c;
    crc->H = MALLOC(k, int*);
    for(int i = 0; i < k; i++){
        crc->H[i] = MALLOC(c, int);
        for (int j = 0; j < c; j++){
            crc->H[i][j] = Binary_Distribution(0.5);
        }
    }
    return crc;
}

void crc_dele(CRC* crc){
    if (crc==NULL) return;
    for(int i = 0; i < crc->k; i++){
        FREE(crc->H[i]);
    }
    FREE(crc->H);
    FREE(crc);
}

void set_crc(int* message_crc, CRC* crc){
    if (crc==NULL) return;
    int* check = message_crc + crc->k;
    for(int i = 0; i < crc->c; i++){
        check[i] = 0;
        for(int j = 0; j < crc->k; j++){
            check[i] ^= (message_crc[j] * crc->H[j][i]);
        }
    }
}

int crc_check(const int* message_crc, CRC* crc){
    if (crc==NULL) return 1;
    int  check_bit;
    for(int i = 0, check = crc->k; i < crc->c; i++, check++){
        check_bit = 0;
        for (int j = 0; j < crc->k; j++){
            check_bit ^= (message_crc[j] * crc->H[j][i]);
        }
        if (check_bit != message_crc[check]){
            return 0;
        }
    }
    return 1;
}

#endif // #ifndef CRC_INCLUDE