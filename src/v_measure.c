#include "R.h"
#include "math.h"

void pecok_gamma(double *IPS,double *n_xc_xd, int *dimension, double *vm)
{
    int p = *dimension;
    /* return UPPER diagonal matrix of V measures b/c in COLMAJOR ORDER */
    for(int b=0; b < p; b++)
    {
        for(int a=0; a < b; a++)
        {
            double vmax = 0;
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
}
