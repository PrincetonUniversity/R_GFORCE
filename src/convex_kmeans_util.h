#ifndef CONVEX_KMEANS_UTIL_H
#define CONVEX_KMEANS_UTIL_H

#include "parallel_rng.h"
#include "time.h"

// CONSTANTS
extern const double DUAL_EPS1_DEFAULT, DUAL_EPS2_DEFAULT, DUAL_Y_T_MIN_DEFAULT;


// NEW TYPES
typedef struct mem_pool {
    void** base;
    int length;
    int start_idx;
    int end_idx;
} mem_pool;

typedef struct workspace {
    double* dwork;
    int* iwork;
    int ldwork;
    int liwork;
    int dsyevd_liwork;
    int dsyevd_liwork_N;
    int dsyevd_ldwork;
    int dsyevd_ldwork_N;

} workspace;

typedef struct problem_instance {
    double* D;
    double* E;
    double* ESI;
    double* D_rsums;
    double TrD;
    double DTD;
    double mu_n;
    int d;
    int K;
} problem_instance;

// External Functions


// GENERAL UTILITIES, MEMORY MANAGEMENT
double min_array(int d,double* V);
void* mem_pool_remove(mem_pool* pool);
void mem_pool_insert(mem_pool* pool, void* mem_ptr);
void allocate_workspace_pd(int d, int K, workspace* work);
void initialize_problem_instance(double* D, double* E, double* ESI, double mu,
                                int d,int K, problem_instance* prob);
double time_difference_ms(clock_t start, clock_t end);
void random_shuffle(int n,int* shuffled);
void random_shuffle_threadsafe(int n,int* shuffled,threadsafe_rng rng_func);
void initialize_identity_matrix(double* restrict I, int d);

// VECTOR, MATRIX OPS
// Computes A = A + c, where A is a vector of length d, c is scalar
void daps(double* restrict A,int inc_A, double c, int d);
// Computes A = AB, where A is a square matrix, and B is diagonal if side=='R'
// Computes A = BA, where A is a square matrix, and B is diagonal if side=='L'
// B is the diagonal of matrix B
// column major layout means 'R' should be FASTER
void dsmtd(double* restrict A,double* restrict B, int d, const char side);
// Computes entrywise exponential function of A
void dvexp(double* restrict A, int d);
// Computes sum of all entries in vector A, of length d
double dsumv(double* restrict A, int d);
// Computes trace of square matrix A, which is dxd
double dtrace(double* restrict A, int d);
// Compute row or column sums of dXd matrix A, store result in A_rsums
void dcsum(double* restrict A, int d, double* restrict A_rsums);
// Compute X + Y and store in Z
void dxpyez(int d, double* restrict X, double* restrict Y, double* restrict Z);
// Compute y = a*x + b*y
void daxpby(double a, double* restrict X, double b, double* restrict Y, int d);
// Compute x = tr(A * B(G))
double dabgtp(double* restrict A, int* restrict G, int d, int K,int* iwork,double* dwork);


// PGD HELPERS
double clust_to_opt_val(problem_instance* prob, int* ga_hat, workspace* work);
void smoothed_gradient(problem_instance* prob, double* X, double* GX_t, double* GS_t,
                        workspace* work);
void smoothed_gradient_nok(problem_instance* prob, double* X, double* GX_t, double* GS_t,
                        workspace* work);
void smoothed_objective(problem_instance* prob, double* X, double* lambda_min,
                        double* obj_val, workspace* work);
void smoothed_objective_nok(problem_instance* prob, double* X, double* lambda_min,
                        double* obj_val, workspace* work);
void C_perp_update(problem_instance* prob, double alpha, double* X_t, double* GX_t,
                    double* GS_t, workspace* work);
void C_perp_update_nok(problem_instance* prob, double alpha, double* X_t, double* GX_t,
                    double* GS_t, workspace*work);
void project_E(problem_instance* prob, double* Z, double lmin, double* Z_proj);


// CLUSTERING
void kmeans_dual_solution_impl(int* restrict ga_hat, problem_instance* restrict prob, double eps1, double eps2,
                                double Y_T_min, double* restrict Y_a_r, double* restrict Y_T_r, int* restrict feasible_r,
                                workspace* restrict work);
void kmeans_dual_solution_nok_impl(int* restrict ga_hat, problem_instance* restrict prob, double eps1,
                                   double* restrict Y_a_r, int* restrict feasible_r, workspace* restrict work);
void kmeans_pp_impl(double const* D, int K, int n, int m, int* cluster_assignment_r,
                        double* centers_r, int* num_iters_R, double* time_R, workspace* work);
void relabel_clusters(int* restrict clusters,int K);


#endif
