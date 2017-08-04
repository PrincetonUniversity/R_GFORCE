#include "R.h"
#include "math.h"
#include "string.h"
#include "float.h"
#include "R_ext/Lapack.h"
#include "R_ext/BLAS.h"
#include "FORCE.h"
// #include "convex_kmeans_util.h"


// void hclust(double* dist,int n,int m,int* agglomerate_idx_1, int* agglomerate_idx_2, double* agglomerate_dmin,double* dwork,int* iwork);

void hclust_agglomerate_R(double* data, int* n0, int* m0, int* agglomerate_idx_1, int* agglomerate_idx_2, double* agglomerate_dmin){
    hclust_agg_t hclust_sol;
    int n = *n0;
    int m = *m0;
    int* iwork;
    double* dwork;
    hclust_sol.agg_idx_min = agglomerate_idx_1;
    hclust_sol.agg_idx_max = agglomerate_idx_2;
    hclust_sol.agg_dist = agglomerate_dmin;
    dwork = (double *) R_alloc(n,sizeof(double));
    iwork = (int *) R_alloc(2*n,sizeof(int));

    hclust_agglomerate(data,n,m,&hclust_sol,dwork,iwork);
}


void hclust_agg2clust();


// returns hclust_t
// 
void hclust();




//requires dwork of length at least n
//requires iwork of length at least 2n
void hclust_agglomerate(double* dist,int n,int m,hclust_agg_t* hclust_sol,double* dwork,int* iwork) {
    // Step 0 - Declarations
    int* nn_idx; //nearest neighbor idx
    int* active;
    int* agglomerate_idx_1;
    int* agglomerate_idx_2;
    double* agglomerate_dmin;
    double* nn_dist; //nearest neighbor distance
    int n2,itmp1,itmp2,idx_min,idx_max,dmin_idx1,dmin_idx2,num_clust;
    double dtmp1,dtmp2,dmin;

    // Step 1 - Init
    nn_dist = dwork;
    nn_idx = iwork;
    active = iwork + n;
    n2 = n*n;
    agglomerate_idx_1 = hclust_sol -> agg_idx_min;
    agglomerate_idx_2 = hclust_sol -> agg_idx_max;
    agglomerate_dmin = hclust_sol -> agg_dist;

    // Step 2 - Active Flags

    for(int i=0; i < n; i++){
        active[i] = 1;
    }


    // Step 3 - Initial nearest neighbors
    // This finds nearest neighbor with larger index
    for(int i=0; i < n-1;i++) {
        dtmp1 = DBL_MAX; //min distance
        itmp1 = -1; //min dist idx
        for(int j=i+1; j < n; j++) {
            itmp2 = i*n + j; //offset
            if(dtmp1 > dist[itmp2]) {
                dtmp1 = dist[itmp2];
                itmp1 = j;
            }

        }
        nn_idx[i] = itmp1;
        nn_dist[i] = dtmp1;
    }

    // Step 4 - 
    num_clust = n;

    while(num_clust > 1) {

        // find closest distance
        dmin = DBL_MAX;
        dmin_idx1 = -1;
        dmin_idx2 = -1;   
        for(int i=0; i < n-1; i++){
            if(active[i] == 1 && nn_dist[i] < dmin) {
                dmin = nn_dist[i];
                dmin_idx1 = i;
                dmin_idx2 = nn_idx[i];
            }
        }
        

        // add to list of agglomerations
        idx_min = dmin_idx1 > dmin_idx2 ? dmin_idx2 : dmin_idx1; // smaller idx
        idx_max = dmin_idx1 > dmin_idx2 ? dmin_idx1 : dmin_idx2; // larger idx
        
        agglomerate_idx_1[n-num_clust] = idx_min;
        agglomerate_idx_2[n-num_clust] = idx_max;
        agglomerate_dmin[n-num_clust] = dmin;

        // deactivate larger idx (bc we found nearest neighbors to the right)
        active[idx_max] = 0;

        // update distances from the new merged cluster
        dmin = DBL_MAX;
        for(int i=0; i < n; i++) {
            if(active[i] == 1 && i != idx_min) {
                //update dist matrix
                dtmp1 = dist[idx_min*n + i];
                dtmp2 = dist[idx_max*n + i];
                dtmp1 = dtmp1 > dtmp2 ? dtmp1 : dtmp2;
                dist[idx_min*n+ i] = dtmp1;
                dist[i*n+ idx_min] = dtmp1;

                //update dmin for idx_min to the right
                //check if need to update nn_idx,nn_dist for i
                dtmp1 = dist[idx_min*n + i];
                dtmp2 = nn_dist[i];
                if(idx_min < i && dtmp1 < dmin) {
                    dmin = dtmp1;
                    dmin_idx2 = i;
                } else if(idx_min > i && dtmp1 < dtmp2) {
                    nn_dist[i] = dtmp1;
                    nn_idx[i] = idx_min;
                }
            }
        }

        // update nearest neighbor for idx_min
        nn_idx[idx_min] = dmin_idx2;
        nn_dist[idx_min] = dmin;

        //update nearest neighbor to the right as required
        for(int i=0; i < n-1; i++){
            if(active[i] == 1 && (nn_idx[i] == idx_min || nn_idx[i] == idx_max)) {
                dmin = DBL_MAX;
                dmin_idx1 = -1;
                for(int j = i+1; j < n; j++) {
                    dtmp1 = dist[i*n + j];
                    if(active[j] == 1 && dtmp1 < dmin) {
                        dmin = dtmp1;
                        dmin_idx1 = j;
                    }
                }
                nn_idx[i] = dmin_idx1;
                nn_dist[i] = dmin;
            }
        }


        // decrement cluster number for next agglomeration
        num_clust = num_clust - 1;
    }

}