#include "math.h"
#include "convex_kmeans.h"
#include "convex_kmeans_util.h"
#include "omp.h"
#include "R.h"
#include "R_ext/BLAS.h"

// CONSTANTS
static const int INC1 = 1;


// R Entrypoints
void FORCE_initialization_R(double* D, double* s, int* d, int* K, double* opt_estimate,int* clusters, int* cluster_representation,
                                        double* E, double* X0, double* E_obj, double* X0_obj) {
    int d0 = *d;
    int* iwork = (int *) R_alloc(d0,sizeof(int));
    double* dwork = (double *) R_alloc(d0*d0,sizeof(double));
    FORCE_initialization(D,*s,d0,*K,opt_estimate,clusters,*cluster_representation,E,X0,E_obj,X0_obj,dwork,iwork);
}

void FORCE_initialization_par_R(double* D, double* s, int* d, int* K, double* opt_estimate,int* clusters, int* cluster_representation,
                                        double* E, double* X0, double* E_obj, double* X0_obj) {
    int d0 = *d;
    int* iwork = (int *) R_alloc(d0,sizeof(int));
    double* dwork = (double *) R_alloc(d0*d0,sizeof(double));
    FORCE_initialization_par(D,*s,d0,*K,opt_estimate,clusters,*cluster_representation,E,X0,E_obj,X0_obj,dwork,iwork);
}


// average random permutations of fr_base
// assumes initialization of E to zero
void add_random_shuffle(int d, int num_shuffles, double* E, double* fr_base, int* shuffled){
    double dtmp1;
    int itmp1, itmp2;

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

// average random permutations of fr_base
// assumes initialization of E to zero
void add_random_shuffle_par(int d, int num_shuffles, double* E, double* fr_base, int* shuffled){
    double dtmp1;
    int itmp1, itmp2;
    int d2 = d*d;

    for(int i=0; i < num_shuffles; i++){
        random_shuffle(d,shuffled);

        #pragma omp parallel for
        for(int j=0; j < d2; j++) {
            itmp1 = j % d;
            itmp2 = j / d;
            itmp1 = shuffled[itmp1];
            itmp2 = shuffled[itmp2];
            itmp1 = itmp2*d + itmp1;
            dtmp1 = fr_base[itmp1] + E[j];
            E[j] = dtmp1;
        }
    }

    //rescale all entries by 1/d
    #pragma omp parallel for
    for(int i=0; i < d2; i++){
        E[i] = E[i] / d;
    }
}


/* Return Values are Matrices E, X0 and their objective values */
// iwork should have length d
// fr_base must be at least length d^2
void FORCE_initialization(double* D, double s, int d, int K, double* opt_estimate, int* clusters, int cluster_representation,
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
    if(cluster_representation == 0){
        daxpby(s,opt_estimate,1-s,X0,d2);
    } else {
        dgxpby(s,clusters,K,1-s,X0,d,iwork);
    }
    
    X0_obj_l = F77_NAME(ddot)(&d2,D,&INC1,X0,&INC1);

    // Copy return objective values
    *X0_obj = X0_obj_l;
    *E_obj = E_obj_l;
}

/* Return Values are Matrices E, X0 and their objective values */
// iwork should have length d
// fr_base must be at least length d^2
void FORCE_initialization_par(double* D, double s, int d, int K, double* opt_estimate, int* clusters, int cluster_representation,
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
    add_random_shuffle_par(d,d,E,fr_base,iwork);
    add_random_shuffle_par(d,d,X0,fr_base,iwork);

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
    if(cluster_representation == 0){
        daxpby(s,opt_estimate,1-s,X0,d2);
    } else {
        dgxpby(s,clusters,K,1-s,X0,d,iwork);
    }
    
    X0_obj_l = F77_NAME(ddot)(&d2,D,&INC1,X0,&INC1);

    // Copy return objective values
    *X0_obj = X0_obj_l;
    *E_obj = E_obj_l;
}


// add matrices with B(G) given as G
// d signifies dimension of matrices X and Y
// X is given in cluster representation form
// Y = a*B(G) + b*Y
// iwork needs length at least K+1
// clusters are numbered either 1..K or 0..K-1
void dgxpby(double a, int* restrict G, int K, double b, double* restrict Y, int d, int* group_sizes) {
    // get group sizes
    // Local Vars
    int itmp1,itmp2,itmp3;
    double dtmp1;

    // Zero out
    for(int i=0; i < K+1; i++){
        group_sizes[i] = 0;
    }

    // Group Sizes
    for(int i=0; i < d; i++){
        itmp1 = G[i];
        group_sizes[itmp1] = group_sizes[itmp1] + 1;
    }

    //update Y
    for(int j=0; j < d; j++){
        for(int i=0; i < d; i++){
            itmp1 = G[i];
            itmp2 = G[j];
            if(itmp1 == itmp2){
                itmp3 = group_sizes[itmp1];
                dtmp1 = a/itmp3 + b * (*Y);
            } else {
                dtmp1 = b * (*Y);
            }
            *Y = dtmp1;

            Y++;
        }
    }
}
