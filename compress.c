//Script to process bootstrapped output of the Backus-Gilbert code

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <time.h>
#include <math.h>

void main(){
   
    srand(time(NULL));
 
    int Ns;     //Frequency sample
    int Nb;     //Bootstrap sample

    scanf("%d",&Ns);
    scanf("%d",&Nb);

    double Rho[Nb][Ns];
    double RhoMean[Ns];
    double RhoCov[Ns];

    double RhoErr[Nb][Ns];
    double RhoErrMean[Ns];
    double RhoErrCov[Ns];

    double FreqWidth[Nb][Ns];
    double FreqWidthMean[Ns];
    double FreqWidthCov[Ns];

    double Score[Nb][Ns];
    double ScoreMean[Ns];
    double ScoreCov[Ns];

    // Initialise output arrays

    for (int i=0; i<Ns; i++){RhoMean[i] = 0.0;}
    for (int i=0; i<Ns; i++){RhoCov[i] = 0.0;}

    for (int i=0; i<Ns; i++){RhoErrMean[i] = 0.0;}
    for (int i=0; i<Ns; i++){RhoErrCov[i] = 0.0;}

    for (int i=0; i<Ns; i++){FreqWidthMean[i] = 0.0;}
    for (int i=0; i<Ns; i++){FreqWidthCov[i] = 0.0;}

    for (int i=0; i<Ns; i++){ScoreMean[i] = 0.0;}
    for (int i=0; i<Ns; i++){ScoreCov[i] = 0.0;}

    // Read in bootstrapped results from tempfile

    for (int i=0; i<Nb; i++){
            scanf("%*s");
        for (int j=0; j<Ns; j++){
            double tmp;
            scanf("%lf",&tmp);
            Rho[i][j] = tmp;
        }
        for (int j=0; j<Ns; j++){
            double tmp;
            scanf("%lf",&tmp);
            RhoErr[i][j] = tmp;
        }
        for (int j=0; j<Ns; j++){
            double tmp;
            scanf("%lf",&tmp);
            FreqWidth[i][j] = tmp;
        }
        for (int j=0; j<Ns; j++){
            double tmp;
            scanf("%lf",&tmp);
            Score[i][j] = tmp;
        }
 
    }

    // Compute mean and covariance for bootstrap samples       

    for (int i=0; i<Ns; i++){
        for (int j=0; j<Nb; j++){
            RhoMean[i] += Rho[j][i]/Nb;
        }
    }
    for (int i=0; i<Ns; i++){
        for (int j=0; j<Nb; j++){
            RhoCov[i] += pow(Rho[j][i]-RhoMean[i],2.0)/(Nb-1);
        }
    }
        
    for (int i=0; i<Ns; i++){
        for (int j=0; j<Nb; j++){
            RhoErrMean[i] += RhoErr[j][i]/Nb;
        }
    }
    for (int i=0; i<Ns; i++){
        for (int j=0; j<Nb; j++){
            RhoErrCov[i] += pow(RhoErr[j][i]-RhoErrMean[i],2.0)/(Nb-1);
        }
    }
        
    for (int i=0; i<Ns; i++){
        for (int j=0; j<Nb; j++){
            FreqWidthMean[i] += FreqWidth[j][i]/Nb;
        }
    }
    for (int i=0; i<Ns; i++){
        for (int j=0; j<Nb; j++){
            FreqWidthCov[i] += pow(FreqWidth[j][i]-FreqWidthMean[i],2.0)/(Nb-1);
        }
    }
        
    for (int i=0; i<Ns; i++){
        for (int j=0; j<Nb; j++){
            ScoreMean[i] += Score[j][i]/Nb;
        }
    }
    for (int i=0; i<Ns; i++){
        for (int j=0; j<Nb; j++){
            ScoreCov[i] += pow(Score[j][i]-ScoreMean[i],2.0)/(Nb-1);
        }
    }

    //Print to output   

    for (int i=0; i<Ns; i++){
        printf("%.6f\n",RhoMean[i]);
    }
    for (int i=0; i<Ns; i++){
        printf("%.6f\n",RhoCov[i]);
    }

     for (int i=0; i<Ns; i++){
        printf("%.6f\n",RhoErrMean[i]);
    }
    for (int i=0; i<Ns; i++){
        printf("%.6f\n",RhoErrCov[i]);
    }
    
    for (int i=0; i<Ns; i++){
        printf("%.6f\n",FreqWidthMean[i]);
    }
    for (int i=0; i<Ns; i++){
        printf("%.6f\n",FreqWidthCov[i]);
    }
 
    for (int i=0; i<Ns; i++){
        printf("%.6f\n",ScoreMean[i]);
    }
    for (int i=0; i<Ns; i++){
        printf("%.6f\n",ScoreCov[i]);
    }
 
}
