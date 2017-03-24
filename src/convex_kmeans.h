#ifndef CONVEX_KMEANS_H
#define CONVEX_KMEANS_H

// Type for PGD Options and return values
typedef struct {
    int verbosity;
    int kmeans_iter;
    int dual_frequency;
    int max_iter;
    int finish_pgd;
    int number_restarts;
    int* restarts;
    double alpha;
    double eps_obj;
} pgd_opts;

typedef struct {
    double* Z_T;
    double* B_Z_T;
    double Z_T_lmin;
    double B_Z_T_opt_val;
    double* Z_best;
    double* B_Z_best;
    double Z_best_lmin;
    double B_Z_best_opt_val;
    int* kmeans_best;
    double kmeans_opt_val;
    double kmeans_best_time;
    int kmeans_iter_best;
    int kmeans_iter_total;
    int dc;
    double dc_time;
    int dc_grad_iter;
    int grad_iter_best;
    double grad_iter_best_time;
    double total_time;
} pgd_results;


// C Access points
void kmeans(double* D, int K, int n, int m, int* centers_init, int* cluster_assignment_r, double* centers_r);
void kmeans_pp(double* D, int K, int n, int m, int* cluster_assignment_r, double* centers_r);
void kmeans_dual_solution_primal_min(int* ga_hat, double* D, int K, int d, double eps1, 
                                        double eps2, double Y_T_min, double* Y_a_r,
                                        double* Y_T_r, int* feasible_r);
void primal_dual_adar(double* D, double* D_kmeans, double* E, double* ESI, double* X0, 
                        int d, int K, pgd_opts* opts, pgd_results* results);
void gamma_alternative_estimator(double* restrict IPS,double* restrict ips_diag_sqrt, int d, double scaling,
                                int* restrict nes, double* restrict gamma_hat, double* restrict ne_meas);
void gamma_alternative_estimator_par(double* restrict IPS,double* restrict ips_diag_sqrt, int d, double scaling,
                                        int* restrict nes, double* restrict gamma_hat, double* restrict ne_meas,int* restrict par_idxs);


// R Access Points
void kmeans_pp_R(double* D, int* K0, int* n0, int* m0, int* cluster_assignment_r, double* centers_r);
void primal_dual_adar_R(double* D, double* D_kmeans, double* E, double* ESI, double* X0, int* d, int* K, 
    int* in_verbosity, int* in_kmeans_iter, int* in_dual_frequency, int* in_max_iter,
    int* in_finish_pgd, int* in_number_restarts, int* in_restarts, double* in_alpha, double* in_eps_obj,
    double* out_Z_T, double* out_B_Z_T, double* out_Z_T_lmin, double* out_Z_best, double* out_B_Z_best, double* out_Z_best_lmin,
    double* out_B_Z_T_opt_val, double* out_B_Z_best_opt_val, double* out_kmeans_opt_val,  int* out_kmeans_best, double* out_kmeans_best_time, 
    int* out_kmeans_iter_best, int* out_kmeans_iter_total, int* out_dc, double* out_dc_time, 
    int* out_dc_grad_iter, int* out_grad_iter_best, double* out_grad_iter_best_time, double* out_total_time);
void kmeans_dual_solution_primal_min_R(int* ga_hat, double* D, int* K_0, int *dimension, 
                                        double* eps1_0, double* eps2_0, double* Y_T_min_0, 
                                        double* Y_a_r, double* Y_T_r, int* feasible_r);
void v_measure(double* restrict IPS,double* restrict n_xc_xd, int* dimension, double* restrict vm);
void v_measure_par(double* restrict IPS,double* restrict n_xc_xd, int* dimension, double* restrict vm);
void gamma_alternative_estimator_R(double* restrict IPS,double* restrict ips_diag_sqrt, int* dimension,
                                    double* scaling, int* restrict nes, double* restrict gamma_hat);
void gamma_alternative_estimator_par_R(double* restrict IPS,double* restrict ips_diag_sqrt, int* dimension, 
                                        double* scaling, int* restrict nes, double* restrict gamma_hat);


#endif