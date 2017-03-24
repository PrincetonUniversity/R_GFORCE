#include "math.h"
#include "convex_kmeans.h"

// R Entrypoint
void FORCE_initialization_R(double* D, double* s, int* d, int* K, int* km_estimate,
                                        double* E, double* X0, double* E_obj, double* X0_obj) {
    FORCE_initialization(D,*s,*d,*K,km_estimate,E,X0,E_obj,X0_obj);
}

/* Return Values are Matrices E, X0 and their objective values */
void FORCE_initialization(double* D, double s, int d, int K, int* km_estimate,
                            double* E, double* X0, double* E_obj, double* X0_obj) {
    ;
}