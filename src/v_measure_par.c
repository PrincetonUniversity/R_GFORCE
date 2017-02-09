#include "R.h"
#include "math.h"
#include "omp.h"

/* return UPPER diagonal matrix of V measures b/c in COLMAJOR ORDER */
/* a is row and b is column */
void pecok_gamma_par(double *IPS,double *n_xc_xd, int *dimension, double *vm)
{
    int p = *dimension;
    int num_entries = p*(p-1) / 2;

    #pragma omp parallel for
    for(int entry_number=0; entry_number < num_entries; entry_number++)
    {
        double vmax = 0;
        int a = entry_number % p;
        int b = entry_number / p;
        for(int d=0; d < p; d++)
        {
            for(int c = 0; c < d; c++)
            {
                if(c!=a && c != b && d!=a && d!=b && n_xc_xd[d*p+ c] != 0)
                {
                    double vcd = (IPS[a*p+ c] + IPS[b*p+ d] - IPS[b*p+ c] - IPS[a*p+ d]) / n_xc_xd[d*p+ c];
                    double vcd_abs = (vcd) > 0 ? vcd : -1*vcd;
                    vmax = (vcd_abs) > (vmax) ? vcd_abs : vmax;
                }
            }
        }
        vm[b*p+ a] = vmax;
    }
}
