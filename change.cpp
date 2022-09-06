#include <stdio.h>
#include "construction.h"

int main(){
    double snr = 2.00;
    int u = 250000;
    
    
    int n;
    int k;
    double rate;

    char* family = (char*)"ABS_Plus";
    scanf("%d %d %lf", &n, &k, &rate);
    printf("n = %d, k = %d, R = %f\n", n, k, rate);
    
    char* oldfile = MALLOC(100, char);
    sprintf(oldfile, "oldfile/%s_BPSK-AWGN_%.3fdB_n%d_k%d_u%d.txt", family, snr, n, k, u);
    char* newfile = MALLOC(100, char);
    sprintf(newfile, "consfile/%s_BPSK-AWGN_%.3fdB_n%d_R%.1f_u%d.txt",family, snr, n, rate, u);
    
    FILE* fp = fopen(oldfile, "r");
    
    if(fp){
        
        int m = LOG2(n);
        int** state = MALLOC(m, int*);
        state[0] = NULL;
        for(int cur_m = 1; cur_m < m  ; cur_m++){
            int length = ((1<<(m-cur_m))-1);
            state[cur_m] = MALLOC(length, int);
            for(int branch = 0; branch < length; branch++){
                fscanf(fp, "%d", &state[cur_m][branch]);
            }
        }
        
        double* conditional_entropy = MALLOC(n, double);
        for(int i = 0; i < n; i++){
            fscanf(fp, "%lf", &conditional_entropy[i]);
        }
        
        index_value* iv = MALLOC(n, index_value);
        for (int i = 0; i < n; i++){
            iv[i].index = i;
            iv[i].value = conditional_entropy[i];
        }
        qsort(iv, n, sizeof(index_value), index_cmp);

        FILE* nfp = fopen(newfile, "a");
            for(int cur_m = 1; cur_m < m  ; cur_m++){
                int length = ((1<<(m-cur_m))-1);
                for(int branch = 0; branch < length; branch++){
                    fprintf(nfp, "%d ", state[cur_m][branch]);
                }
                fprintf(nfp, "\n");
            }
            
            for(int i = 0; i < n; i++){
                fprintf(nfp, "%4d  %.16e\n", iv[i].index, iv[i].value);
            }
        fclose(nfp);
        FREE(iv);

        for(int i = 1; i < m; i++)free(state[i]);
        free(state);
        free(conditional_entropy);
    }else{
        printf("There is no file named %s!\n", oldfile);
        exit(0);
    }

    fclose(fp);
    FREE(oldfile);
    FREE(newfile);

    return 0;
}