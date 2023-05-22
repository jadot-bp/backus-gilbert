#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <time.h>
#include <math.h>

/*
 * covgend
 *
 * Calculates the mean and the covariance of a bootstrap sample of the 
 * configuration set for the Euclidean correlator.
 *
 * ============================================================================
 */

int main(){
   
    srand(time(NULL));  //Initialise pseudo-RNG
 
    int Nt;     //Number of time points
    int Nc;     //Number of configurations
    int Nb;     //Number of (bootstrap) samples

    char line[255];    
    char tmp[255]; 
    char *strptr;
    FILE *fptr;

    fptr = fopen("seed","a");
    fprintf(fptr,"seed:\n");

    scanf("%d",&Nt);
    scanf("%d",&Nc);
    scanf("%d",&Nb);

    double G[Nc][Nt];           //Matrix of all correlators
    double Gboot[Nb][Nt];       //Matrix of bootstrap samples of correlators
    double GMean[Nt];           //Mean correlator value for bootstrap samples
    double GCov[Nt][Nt];        //Covariance matrix of boostrap samples

    /* Initialise arrays */

    for (int i=0; i<Nt; i++){GMean[i] = 0.0;}
    for (int i=0; i<Nt; i++){
        for (int j=i+1; j<Nt; j++){
            GCov[i][j] = GCov[j][i] = 0.0;
        }
    }

    /* Load in G */

    for (int i=0; i<Nc; i++){
        for (int j=0; j<Nt; j++){
            double tmp;
            scanf("%lf",&tmp);
            G[i][j] = tmp;
        }
    }

    /* Bootstrap sample G */

    int count = 0;

    while(count<Nb){
        int index = rand()%Nc;
        fprintf(fptr,"%d;",index);      //Record choices
        for (int i=0; i<Nt; i++){
            Gboot[count][i] = G[index][i];
        }
        count++;
    }

    fclose(fptr);

    /* Calculate GMean */

    for (int i=0; i<Nt; i++){
        for (int j=0; j<Nb; j++){
            GMean[i] += Gboot[j][i]/Nb;
        }
    }

    /* Calculate Covariance in G */

    for (int i=0; i<Nt; i++){
        double gi_mean = 0.0;
        for(int gi=0; gi<Nb; gi++){
            gi_mean += Gboot[gi][i];
        }
        gi_mean /= Nb;

        for (int j=0; j<Nt; j++){
            double gj_mean = 0.0;
            for(int gj=0; gj<Nb; gj++){
                gj_mean += Gboot[gj][j];
            }
            gj_mean /= Nb;

            GCov[i][j] = 0.0;
            for (int k=0; k<Nc; k++){
                GCov[i][j] += (Gboot[k][i]-gi_mean)*(Gboot[k][j]-gj_mean);
            }
            GCov[i][j] /= (Nb-1);
        }
    }

    /* Output Mean */
    
    for (int i=0; i<Nt; i++){
        printf("%g\n",GMean[i]);
    }
   
    int isdynamic = 1;

    if (isdynamic == 1){

        /* Output Variances (diagonal elements)*/

        for (int i=0; i<Nt; i++){
            printf("%g\n",GCov[i][i]);
        }
        
        /* Output off-diagonals for full error */

        for (int i=0; i<Nt; i++){
            for (int j=i+1; j<Nt; j++){
                printf("%.17g\n",GCov[i][j]);
            }
        }
    }

    return 0;
}
