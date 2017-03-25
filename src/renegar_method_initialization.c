#include "math.h"
#include "convex_kmeans.h"
#include "convex_kmeans_util.h"
#include "R.h"
#include "R_ext/BLAS.h"

// CONSTANTS
static const int INC1 = 1;


// R Entrypoint
void FORCE_initialization_R(double* D, double* s, int* d, int* K, double* opt_estimate,
                                        double* E, double* X0, double* E_obj, double* X0_obj) {
    int d0 = *d;
    int* iwork = (int *) R_alloc(d0,sizeof(int));
    double* dwork = (double *) R_alloc(d0*d0,sizeof(double));
    FORCE_initialization(D,*s,d0,*K,opt_estimate,E,X0,E_obj,X0_obj,dwork,iwork);
}


// average random permutations of fr_base
// assumes initialization of E to zero
void add_random_shuffle(int d, int num_shuffles, double* E, double* fr_base, int* shuffled){
    double dtmp1;
    volatile int itmp1;
    volatile int itmp2;

    for(int i=0; i < num_shuffles; i++){
        random_shuffle(d,shuffled);
        for(int b=0; b < d; b++){
            for(int a=0; a < d; a++) {
                itmp1 = shuffled[a];
                itmp2 = shuffled[b];
                itmp1 = itmp2*d + itmp1;
                itmp2 = b*d + a;
                dtmp1 = fr_base[itmp1] + E[itmp2];
                E[itmp2] = dtmp1;
            }
        }
    }

    //rescale all entries by 1/d
    for(int i=0; i < d*d; i++){
        E[i] = E[i] / d;
    }

}


/* Return Values are Matrices E, X0 and their objective values */
// iwork should have length d
// fr_base must be at least length d^2
void FORCE_initialization(double* D, double s, int d, int K, double* opt_estimate,
                            double* E, double* X0, double* E_obj, double* X0_obj, double* fr_base, int* iwork) {
    double X0_obj_l, E_obj_l, dtmp1;
    int d2 = d*d;

    full_rank_feasible(d,K,fr_base);

    /* zero out E and X0 arrays */
    for(int i=0; i < d2; i++){
        E[i] = 0;
        X0[i] = 0;
    }

    /* average random shuffles of indices */
    add_random_shuffle(d,d,E,fr_base,iwork);
    add_random_shuffle(d,d,X0,fr_base,iwork);

    /* get corresponding objective values */
    X0_obj_l = F77_NAME(ddot)(&d2,D,&INC1,X0,&INC1);
    E_obj_l = F77_NAME(ddot)(&d2,D,&INC1,E,&INC1);

    // X0 should have smaller objective value
    if(X0_obj_l > E_obj_l) {
        memcpy(fr_base,E,d2*sizeof(double));
        memcpy(E,X0,d2*sizeof(double));
        memcpy(X0,fr_base,d2*sizeof(double));
        E_obj_l = X0_obj_l;
    }

    // mix with opt_estimate
    daxpby(s,opt_estimate,1-s,X0,d2);
    X0_obj_l = F77_NAME(ddot)(&d2,D,&INC1,X0,&INC1);

    // Copy return objective values
    *X0_obj = X0_obj_l;
    *E_obj = E_obj_l;
}