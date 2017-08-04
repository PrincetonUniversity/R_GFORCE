#include "float.h"
#include "math.h"
#include "time.h"
#include "R.h"
#include "R_ext/BLAS.h"
#include "R_ext/Lapack.h"
#include "convex_kmeans_util.h"


// CONSTANTS
const double DUAL_EPS1_DEFAULT = 0.01;
const double DUAL_EPS2_DEFAULT = 0.000001;
const double DUAL_Y_T_MIN_DEFAULT = 0.01;
static const char JOBZV = 'V';
static const char JOBZN = 'N';
static const char UPLO = 'U';
static const int INC1 = 1;


// GENERAL UTILITIES, MEMORY MANAGEMENT
double min_array(int n, double* V){
    double min_val = DBL_MAX;
    for(int i=0; i < n; i++) min_val = (V[i] < min_val) ? V[i] : min_val;
    return min_val;
}

void mem_pool_insert(mem_pool* pool, void* mem_ptr){
    int tmp_idx = (pool->end_idx);
    if(tmp_idx == (pool->length) -1) {
       tmp_idx = -1;
    }
    *((pool->base) + ++tmp_idx) = mem_ptr;
    pool->end_idx = tmp_idx;
}

void* mem_pool_remove(mem_pool* pool){
    int tmp_idx = pool->start_idx;
    void* mem_ptr = *((pool->base) + tmp_idx++);
    if(tmp_idx == pool->length){
        tmp_idx = 0;
    }
    pool->start_idx = tmp_idx;
    return mem_ptr;
}

void random_shuffle(int n,int* shuffled){
    int itmp1,itmp2;
    double dtmp1;

    for(int i=0; i<n; i++) {
        shuffled[i] = i;
    }
    GetRNGstate();
    for(int i=0; i < n; i++){
        dtmp1 = unif_rand();
        dtmp1 = dtmp1*(n-0.02)  - 0.49;
        itmp1 = round(dtmp1);
        itmp2 = shuffled[i];
        shuffled[i] = shuffled[itmp1];
        shuffled[itmp1] = itmp2;
    }

    PutRNGstate();
}

void random_shuffle_threadsafe(int n,int* shuffled,threadsafe_rng tsrng){
    int itmp1,itmp2;
    double dtmp1 = 0;

    for(int i=0; i<n; i++) {
        shuffled[i] = i;
    }
    for(int i=0; i < n; i++){
        dtmp1 = threadsafe_rng_next(tsrng);
        // Rprintf("Next Number -- %f\n",dtmp1);
        dtmp1 = dtmp1*(n-0.02)  - 0.49;
        itmp1 = round(dtmp1);
        itmp2 = shuffled[i];
        shuffled[i] = shuffled[itmp1];
        shuffled[itmp1] = itmp2;
    }
}

// Allocates properly sized workspace for primal_dual_adar
void allocate_workspace_pd(int d, int K, workspace* work){
    // local vars
    int liwork = -1;
    int liwork_N = -1;
    int ldwork = -1;
    int ldwork_N = -1;
    int tmp1;

    // USE CASE 0 -- Find Space Needed for both types of DSYEVD Calls
    double* X;
    double dwork0[2];
    int iwork0[2];
    double dvec[2];
    int lapack_info = 0;
    F77_CALL(dsyevd)(&JOBZN,&UPLO,&d,X,&d,dvec,dwork0,&ldwork_N,iwork0,&liwork_N,&lapack_info);
    liwork_N = *iwork0;
    ldwork_N = (int) *dwork0;
    work -> dsyevd_liwork_N = liwork_N;
    work -> dsyevd_ldwork_N = ldwork_N;
    F77_CALL(dsyevd)(&JOBZV,&UPLO,&d,X,&d,dvec,dwork0,&ldwork,iwork0,&liwork,&lapack_info);
    liwork = *iwork0;
    ldwork = (int) *dwork0;
    work -> dsyevd_liwork = liwork;
    work -> dsyevd_ldwork = ldwork;


    // USE CASE 1 -- SMOOTHED_GRADIENT COMPUTATION
    ldwork = ldwork + d; //add space for w to hold eigvals

    // Workspace size for smoothed gradient (ignoring) syevd is 2d^2+d
    tmp1 = 2*d*d + d;
    ldwork = ldwork > tmp1 ? ldwork : tmp1;


    // USE CASE 2 -- C_perp_update and clust_to_opt_val
    // Get workspace size needed for C_perp_update and clust_to_opt_val
    // These require ldwork >= d^2+5d+2 and liwork >= d+3K+3
    tmp1 = d*d + 5*d + 2;
    ldwork = ldwork > tmp1 ? ldwork : tmp1;
    tmp1 = d + 3*K + 3;
    liwork = liwork > tmp1 ? liwork : tmp1;
    

    // USE CASE 3 -- kmeans_pp
    // STRICTLY DOMINATED IN REQUIREMENTS BY CASE 2


    // USE CASE 4 -- smoothed_objective
    // STRICTLY DOMINATED BY REQUIREMENTS IN USE CASE 2 (should be)
    tmp1 = ldwork_N + 2*d*d + d;
    ldwork = ldwork > tmp1 ? ldwork : tmp1;
    liwork = liwork > liwork_N ? liwork : liwork_N;


    // Initialize workspace
    work -> dwork = (double *) R_alloc(ldwork,sizeof(double));
    work -> ldwork = ldwork;
    work -> iwork = (int *) R_alloc(liwork,sizeof(int));
    work -> liwork = liwork;
}

// Takes input and precomputes values for problem instance
// computes the fields D_rsums, TrD, DTD
// Assumes that these are un-initialized
void initialize_problem_instance(double* D, double* E, double* ESI, double mu, int d,
                                    int K, problem_instance* prob){
    int d2 = d*d;
    double* D_rsums = (double *) R_alloc(d,sizeof(double));
    dcsum(D,d,D_rsums);

    prob -> d = d;
    prob -> D = D;
    prob -> E = E;
    prob -> ESI = ESI;
    prob -> D_rsums = D_rsums;
    prob -> mu_n = -1*mu;
    prob -> K = K;
    prob -> TrD = dtrace(D,d);
    prob -> DTD = F77_NAME(ddot)(&d2,D,&INC1,D,&INC1);
}

void initialize_identity_matrix(double* restrict I,int d){
    int d2 = d*d;
    int dp1 = d+1;
    double* ptmp1 = I;
    for (int i=0; i < d2; i++){
        *ptmp1 = 0;
        ptmp1++;
    }
    ptmp1 = I;
    for(int i=0; i < d; i++){
        *ptmp1 = 1;
        ptmp1 = ptmp1 + dp1;
    }
}

// VECTOR, MATRIX OPS
void daps(double* restrict A, int inc_A, double c, int d){
    double dtmp;
    if(inc_A == 1){
        for(int i=0; i < d; i++){
            A[i] = A[i] + c;
        }
    } else{
        for(int i=0; i < d; i++){
            dtmp = *A + c;
            *A = dtmp;
            A = A + inc_A;
        }
    }
}

// performs y = a*x + b*y
void daxpby(double a, double* restrict X, double b, double* restrict Y, int d){
    double dtmp1;
    for(int i=0; i < d; i++){
        dtmp1 = a * (*X) + b * (*Y);
        *Y = dtmp1;
        X++;
        Y++;
    }
}

void dsmtd(double* restrict A, double* restrict B, int d, const char side){
    int stride = 1;
    if(side == 'R'){
        for(int i=0; i <d; i++){
            F77_NAME(dscal)(&d,B+i,A + d*i,&stride);
        }
    }else {
        stride = d;
        for(int i=0; i <d; i++){
            F77_NAME(dscal)(&d,B+i,A + i,&stride);
        }
    }
}

void dvexp(double* restrict A, int d) {
    double tmp;
    for(int i = 0; i < d; i++){
        tmp = *A;
        *A = exp(tmp);
        A++;
    }
}

double dsumv(double* restrict A, int d){
    double acc = 0;
    for(int i=0; i < d; i++){
        acc = acc + *A;
        A++;
    }
    return acc;
}

double dtrace(double* restrict A, int d){
    double trace = 0;
    for(int i=0; i < d; i++){
        trace = trace + *A;
        A = A + d + 1;
    }
    return trace;
}


void dcsum(double* restrict A, int d, double* restrict A_csums){
    // zero A_csums
    double* tmp_ptr = A_csums;
    double dtmp;
    for(int i=0; i <d; i++){
        *tmp_ptr = 0;
        tmp_ptr++;
    }
    
    for(int i=0; i < d; i++){
        for(int j=0; j < d; j++){
            dtmp = *A_csums + *A;
            *A_csums = dtmp;
            A++;
        }
        A_csums++;
    }
}

void dxpyez(int d, double* restrict X, double* restrict Y, double* restrict Z){
    double dtmp1;
    for(int i=0; i < d; i++){
        dtmp1 = *X + *Y;
        *Z = dtmp1;
        X++;
        Y++;
        Z++;
    }
}

// Compute x = tr(A * B(G))
// Assumes clusters labeled 1..K or 0..K-1
// REQUIRES work -> dwork be of length K
// REQUIRES work -> iwork be of length d+3K+3
double dabgtp(double* restrict A, int* restrict G, int d, int K,int* iwork,double* dwork) {
    // Local Vars
    double opt_val = 0.0;
    int tmp1,tmp2,tmp3,tmp4;
    double dtmp1;
    int Kp1 = K+1;

    // Get Memory From Workspace
    double* group_sums = dwork;
    int* group_sizes = iwork;
    int* group_tailp1_idx = group_sizes + Kp1;
    int* group_start_idx = group_tailp1_idx + Kp1; 
    int* group_idxs = group_start_idx + Kp1; // stores by group all idxs in that group

    // Zero out
    for(int i=0; i < Kp1; i++){
        group_sizes[i] = 0;
        group_sums[i] = 0;
    }

    // Group Sizes, zero out Y_a_base
    for(int i=0; i < d; i++){
        tmp1 = G[i];
        group_sizes[tmp1] = group_sizes[tmp1] + 1;
    }

    // Group Bounds Initialization -- points to current last group member, first group member
    tmp1 = 0;
    for(int i=0; i < Kp1; i++){
        group_tailp1_idx[i] = tmp1;
        group_start_idx[i] = tmp1;
        tmp1 += group_sizes[i];
    }

    // Group Membership
    for(int i=0; i < d; i++){
        tmp1 = G[i]; //store group membership of i
        tmp2 = group_tailp1_idx[tmp1];
        group_tailp1_idx[tmp1] = tmp2 + 1; //store location to add this index, then increment
        group_idxs[tmp2] = i;
    }

    // Group Sums
    for(int i=0; i < d; i++){
        tmp1 = G[i]; // group membership of i
        tmp2 = group_start_idx[tmp1]; //start idx of group
        tmp3 = group_tailp1_idx[tmp1]; // end idx

        for(int j=tmp2; j < tmp3; j++){
            tmp4 = (group_idxs[j])*d + i; //Access is in column major form
            dtmp1 = group_sums[tmp1] + A[tmp4];
            group_sums[tmp1] = dtmp1;
        }
    }

    // Objective Value
    opt_val = 0;
    for(int i=0; i < Kp1; i++){
        tmp1 = group_sizes[i];
        if(tmp1 > 0){
            opt_val = opt_val + (group_sums[i] / ((double) tmp1));
        }
    }

    return opt_val;

}

void relabel_clusters(int* restrict clusters,int K) {
    ;
}


double time_difference_ms(clock_t start, clock_t end){
    double diff;
    int msec_diff = end-start;
    diff = (((double) msec_diff))/CLOCKS_PER_SEC;
    return diff; 
}

double euclidean_distance(double const* restrict p1, double const* restrict p2, int m, double* restrict euclidean_distance_tmp) {
    double acc = 0.0;
    double dtmp1;

    #pragma omp simd
    for(int i=0; i < m; i++){
        //it seems to vectorize better this way...
        euclidean_distance_tmp[i] = p1[i] - p2[i];
        acc = acc + euclidean_distance_tmp[i] * euclidean_distance_tmp[i] ;
    } 

    return acc;
}