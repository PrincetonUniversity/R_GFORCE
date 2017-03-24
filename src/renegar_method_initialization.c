#include "math.h"
// #include "R.h"
#include "convex_kmeans.h"

// R Entrypoint
void full_rank_feasible_R(int* d, int* K, double* E) {
    full_rank_feasible(*d,*K,E);
}

/* Return Value stored as E -- full rank feasible matrix for Peng-Wei SDP */
void full_rank_feasible(int d, int K, double* E) {
    int itmp1, itmp2, itmp3;
    int c_low = d/(K-1);
    int c_high = ceil(((double) d)/(K - 1));
    double big_group_val = 1.0/((double) (d-K+1));

    // Zero E Matrix
    for(int i=0; i < d*d; i++){
        E[i] = 0.0;
    }

    // Average B(G) for various G
    for(int i=0; i < c_low; i++){
        // entries for size 1 groups
        int start_1_groups_idx = i*(K-1);
        int end_1_groups_idx = (i+1)*(K-1) - 1;
        for(int j=0; j < K-1; j++){
            itmp1 = start_1_groups_idx + j; //which diagonal element to choose
            itmp2 = itmp1*d + itmp1; //index of E_(itmp1,itmp1)
            E[itmp2]++;
        }
        // iterate over all pairs of indices within the big group
        // do this by skipping the K-1 indices that are their own group 
        for(int a=0; a < d; a++) {
            //skip singleton groups
            if(a == start_1_groups_idx){
                a = end_1_groups_idx;
                continue;
            }
            for(int b=0; b < d; b++) {
                //skip singleton groups
                if(b == start_1_groups_idx){
                    b = end_1_groups_idx;
                    continue;
                }
                itmp1 = b*d + a;
                E[itmp1] = E[itmp1] + big_group_val;
            }
        }
    }

    if(c_low - c_high < 0){
       // entries for size 1 groups
       int start_1_groups_idx = d-K+1;
       int end_1_groups_idx = d-1;
       for(int j=0; j < K-1; j++){
           itmp1 = start_1_groups_idx + j; //which diagonal element to choose
           itmp2 = itmp1*d + itmp1; //index of E_(itmp1,itmp1)
           E[itmp2]++;
       }
       // iterate over all pairs of indices within the big group
       // do this by skipping the K-1 indices that are their own group 
       for(int a=0; a < d; a++) {
           //skip singleton groups
           if(a == start_1_groups_idx){
               a = end_1_groups_idx;
               continue;
           }
           for(int b=0; b < d; b++) {
               //skip singleton groups
               if(b == start_1_groups_idx){
                   b = end_1_groups_idx;
                   continue;
               }
               itmp1 = b*d + a;
               E[itmp1] = E[itmp1] + big_group_val;
           }
       }
    }

    // RESCALE ALL ENTRIES
    for(int a=0; a < d; a++){
        for(int b=0; b < d; b++){
            itmp1 = b*d + a;
            E[itmp1] = E[itmp1] / c_high;
        }
    }
}