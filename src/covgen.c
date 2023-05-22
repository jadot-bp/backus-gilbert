#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

/*
 * covgen
 *
 * Calculates the mean and covariance of the entire set of Nc configurations
 * of the Euclidean correlator.
 *
 * ============================================================================
 */

int main(){
    
    int Nt;     //Number of time points
    int Nc;     //Number of configurations
    
    scanf("%d",&Nt);
    scanf("%d",&Nc);

    double G[Nc][Nt];           //Matrix of all correlators
    double GMean[Nt];           //Mean correlator value at each time point
    double GCov[Nt][Nt];        //Covariance of matrix of correlators

    /* Initialise empty arrays */

    for (int i=0; i<Nt; i++){GMean[i] = 0.0;}

    for (int i=0; i<Nt; i++){
        for (int j=0; j<Nt; j++){
            GCov[i][j] = 0.0;
        }
    }

    /* Load in G across all configurations*/

    for (int i=0; i<Nc; i++){
        for (int j=0; j<Nt; j++){
            double tmp;
            scanf("%lf",&tmp);
            G[i][j] = tmp;
        }
    }
    
    /* Calculate GMean */

    for (int i=0; i<Nt; i++){
        for (int j=0; j<Nc; j++){
            GMean[i] += G[j][i]/Nc;
        }
    }

    /* Calculate Covariance in G */

    for (int i=0; i<Nt; i++){
        double gi_mean = 0.0;
        for(int gi=0; gi<Nc; gi++){
            gi_mean += G[gi][i];
        }
        gi_mean /= Nc;

        for (int j=0; j<Nt; j++){
            double gj_mean = 0.0;
            for(int gj=0; gj<Nc; gj++){
                gj_mean += G[gj][j];
            }
            gj_mean /= Nc;

            GCov[i][j] = 0.0;
            for (int k=0; k<Nc; k++){
                GCov[i][j] += (G[k][i]-gi_mean)*(G[k][j]-gj_mean);
            }
            GCov[i][j] /= (Nc-1);
        }
    }

    /* Output Mean */
    
    for (int i=0; i<Nt; i++){
        printf("%g\n",GMean[i]);
    }
    
    /* Output variance terms (diagonal elements)*/

    for (int i=0; i<Nt; i++){
        printf("%g\n",GCov[i][i]);
    }
    
    /* Output covariance terms (off-diagonals) for full error */

    for (int i=0; i<Nt; i++){
        for (int j=i+1; j<Nt; j++){
            printf("%.17g\n",GCov[i][j]);
        }
    }

    return 0;
}
