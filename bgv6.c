#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <float.h>
#include <omp.h>
#include <gmp.h>
#include <mpfr.h>

#include "interface.h"

/*
    BACKUS-GILBERT ANALYSIS SCRIPT
    v6.0
    =================================================================
*/
double wmin;
double wmax;
float scaling;
int prec = 128;          //MPFR precision

#define NCORES 4         //Number of cores for MPFR - currently unused

void KFunc(double w, double tau, mpfr_t kfunc, mpfr_t factor);
void WFunc(double w, double tau1, double tau2, mpfr_t wfunc, mpfr_t work);

int main(int argc, char *argv[]){

    clock_t main_start = clock();

    int Nt;                  //Number of time points in lattice (the temperature)
    int Ns;                  //Number of sample slices to calculate omega
    int t1;                  //Initial Euclidean time
    int t2;                  //Final Euclidean time

    FILE *fptr;
    FILE *avgf;
    fptr = fopen("out.report","w");
    avgf = fopen("out.avgf","w");
    fprintf(fptr,"%s","OUTPUT REPORT\n====================\n\n");

    if (argc == 1){
        printf("Argument Error: Must supply at least 1 command line argument! (debug mode)\n");
        return(22);
    }
    if (argc > 7){
        printf("Argument Error: Too many arguments!\n");
        return(7); 
    }

    int g = atoi(argv[1]);       //Debug mode for pre-calculating averaging functions
    int width = 0;               //Apply unimodularity constriant (1) or use pure least-squares (0) -- Least squares has unregulated averaging functions
    int emode = 1;               //Full error mode (1) or partial error mode (0) -- Full errors uses covariance matrix

    scanf("%d",&Nt);
    scanf("%d",&Ns);
    
    /*
    t1 = 1;                      //Default start time (inclusive)
    t2 = 64;                     //Default end time (t+1, exclusive)
    */

    double alpha;                //Whitening parameter

    if (argc == 2){              //Set default parameters
        wmin = -0.1;              //Scanning range
        wmax = 1.0;
        alpha = 0.8;   
        t1 = 0;
        t2 = Nt;
    }else{                       //Read from cmd args
        wmin = atof(argv[2]);
        wmax = atof(argv[3]);
        alpha = atof(argv[4]);
        t1 = atoi(argv[5]);
        t2 = atoi(argv[6]);
        if (t2 == -1){
            t2 = Nt;
        }
    }

    scaling = 0.5;               //Kernel rescaling - 0.5 for S-wave and 1.5 for P-wave

    double lmin = wmin;
    double lmax = 10.0;          //Integration range
    double n = 25*Ns;            //Integral precision for composite-trapezium rule

    double w0s[Ns+1];            //Frequency probe (proxy for position in energy space)
    double tau[t2-t1];           //Temporal position 

    //double KWeight[Nt][Nt];       //Kernel weighting matrix
    //double KConst[Nt];           //Constraint vector
    //double AvgCoeff[Nt];    //Averaging coefficient vector
    double AvgCoeffs[Ns+1][t2-t1];

    double G[t2-t1];               //Correlator data
    double Cov[t2-t1][t2-t1];             //Correlator covariance
    double rho[Ns+1];           //Spectral density estimate
    double errs[Ns+1];          //Spectral density error
    double widths[Ns+1];        //Spectral density frequency width/resolution
    double score;        //Inversion deviation
    double condition;     //Kernel weight matrix condition number
    
    /*
     * =======================================================================
                            BACKUS-GILBERT CODE
                                    v6
     * =======================================================================
    */

    /* Initialise sampling array w0*/
    
    double delta = (double)(wmax-wmin)/Ns;

    for (int i=0; i<=Ns; i++){w0s[i] = wmin+i*delta;}
    
    /* Parse G from input*/ 

    for (int i=0; i<Nt; i++){
        double tmp;
        scanf("%lf",&tmp);
        if (i >= t1 && i < t2){
            G[i-t1] = tmp;
        }
    }

    /* Parse Cov from input */
    
    for (int i=0; i<Nt; i++){
        double tmp;
        scanf("%lf",&tmp);
        if (i >= t1 && i < t2){
            Cov[i-t1][i-t1] = tmp;
        }
    }

    for (int i=0; i<Nt; i++){
        for (int j=i+1; j<Nt; j++){
            double tmp;
            scanf("%lf",&tmp);
            if (i >= t1 && i < t2){
                if(emode==1){
                    Cov[i-t1][j-t1] = Cov[j-t1][i-t1] = tmp;
                }else{
                    Cov[i-t1][j-t1] = Cov[j-t1][i-t1] = 0.0;
                }
            }
        }   
    }

    /* Construct Euclidean time vector */

    for (int i=t1; i<t2; i++){tau[i-t1] = (double)i;}   //Initialise lattice time vector

    /* Constructing the kernel weighting matrix*/

    mpfr_t KWeight[t2-t1][t2-t1];

    mpfr_t wfunc;              //Variable to capture WFunc output
    mpfr_t trapz;              //Running trapezium sum
    mpfr_t work;               //Work variable
    mpfr_t temp;               //Temporary (work) variable
 
    mpfr_init2(trapz,prec);
    mpfr_init2(work,prec);
    mpfr_init2(temp,prec);
    mpfr_init2(wfunc,prec);

    clock_t kweight_start = clock(); 

    for (int i=0; i<t2-t1; i++){       //Diagonal, equal-time elements
        
        mpfr_init2(KWeight[i][i],prec);
        //mpfr_set_d(trapz,0.0,MPFR_RNDN);
        WFunc(lmin,tau[i],tau[i],wfunc,work);
        mpfr_set(trapz,wfunc,MPFR_RNDF);
        mpfr_div_d(trapz,trapz,2.0,MPFR_RNDF);
        //result += WFunc(lmin,w0,tau[i],tau[i])/2.0;

        for (int k=1; k<n; k++){
            //result += WFunc(lmin+(lmax-lmin)*(double)k/n,w0,tau[i],tau[i]);
            WFunc(lmin+(lmax-lmin)*k/n,tau[i],tau[i],wfunc,work),
            mpfr_add(trapz,trapz,wfunc,MPFR_RNDF);
        }
        WFunc(lmax,tau[i],tau[i],wfunc,work);
        mpfr_set(work,wfunc,MPFR_RNDF);
        mpfr_div_d(work,work,2.0,MPFR_RNDF);
        mpfr_add(trapz,trapz,work,MPFR_RNDF);

        //result += WFunc(lmax,w0,tau[i],tau[i])/2.0;
        //result *= (lmax-lmin)/n;
        mpfr_mul_d(trapz,trapz,(lmax-lmin)/n,MPFR_RNDF);

        if (width == 1){
            mpfr_mul_d(trapz,trapz,2.0,MPFR_RNDF);
        }

        mpfr_set(KWeight[i][i],trapz,MPFR_RNDF);
    }
    for (int i=0; i<t2-t1; i++){       //Off-diagonal, time-symmetric elements
        for (int j=i+1; j<t2-t1; j++){
                            
            mpfr_init2(KWeight[i][j],prec);
            mpfr_init2(KWeight[j][i],prec);
            WFunc(lmin,tau[i],tau[j],wfunc,work);
            mpfr_set(trapz,wfunc,MPFR_RNDF);
            mpfr_div_d(trapz,trapz,2.0,MPFR_RNDF);
            for (int k=1; k<n; k++){
                WFunc(lmin+(lmax-lmin)*k/n,tau[i],tau[j],wfunc,work);
                mpfr_add(trapz,trapz,wfunc,MPFR_RNDF);
            }
            WFunc(lmax,tau[i],tau[j],wfunc,work);
            mpfr_set(work,wfunc,MPFR_RNDF);
            mpfr_div_d(work,work,2.0,MPFR_RNDF);
            mpfr_add(trapz,trapz,work,MPFR_RNDF);
            mpfr_mul_d(trapz,trapz,(lmax-lmin)/n,MPFR_RNDF);
           
            if (width == 1){
                mpfr_mul_d(trapz,trapz,2.0,MPFR_RNDF);
            }

            mpfr_set(KWeight[i][j],trapz,MPFR_RNDF);
            mpfr_set(KWeight[j][i],trapz,MPFR_RNDF);
        }
    }
    
    /* Whitening kernel weighting matrix*/
    
    for (int i=0; i<t2-t1; i++){
        for (int j=0; j<t2-t1; j++){
            mpfr_mul_d(KWeight[i][j],KWeight[i][j],alpha,MPFR_RNDF);
            mpfr_add_d(KWeight[i][j],KWeight[i][j],(1-alpha)*Cov[i][j],MPFR_RNDF);
        }
    }

    /* Invert weighting matrix */
            
    mpfr_t KCopy[t2-t1][t2-t1];
    mpfr_t KInverse[t2-t1][t2-t1];

    for(int i=0; i<t2-t1; i++){
        for(int j=0; j<t2-t1; j++){
            mpfr_init2(KCopy[i][j],prec);
            mpfr_set(KCopy[i][j],KWeight[i][j],MPFR_RNDF);   //Initialise MPFR copy of KWeight for SVD
            mpfr_init2(KInverse[i][j],prec);
        }
    }

    mpfr_t S[t2-t1];           //Singular elements

    for(int i=0; i<t2-t1; i++){
        mpfr_init2(S[i],prec);
    }
    c_zcall(prec,t2-t1,*KCopy,*KInverse,S); //Pipe matrices to ZKCM interface

    condition = fabsl(mpfr_get_ld(S[0],MPFR_RNDN)/mpfr_get_ld(S[t2-t1-1],MPFR_RNDN));         //Save condition number

    /* Check quality of inversion */
    /*
    {   
        mpfr_t Identity[t2-t1][t2-t1];    
        double score = 0.0;                 
        double dev;
        //printf("KxK^-1:\n");
        for (int i=0; i<t2-t1; i++){
            for (int j=0; j<t2-t1; j++){
                mpfr_init2(Identity[i][j],prec);
                mpfr_set_d(Identity[i][j],0.0,MPFR_RNDN);
                for (int k=0; k<t2-t1; k++){
                    mpfr_mul(temp,KInverse[i][k],KWeight[k][j],MPFR_RNDF);
                    mpfr_add(Identity[i][j],Identity[i][j],temp,MPFR_RNDF);
                }
                dev = mpfr_get_d(Identity[i][j],MPFR_RNDF); 
                if(i==j){
                    score += pow(1.0-dev,2);
                }else{
                    
                    score += dev*dev;
                }
            }
        }
        score = score/((t2-t1)*(t2-t1));
        mpfr_clear(temp);
    
    }//Inversion container 
    */
    mpfr_clear(temp);
    mpfr_clear(trapz);
    mpfr_clear(work);
    mpfr_clear(wfunc);

    mpfr_t SConst[t2-t1];        //Backus-Gilbert spread constraint

    if (width == 1){
        /* Construct Backus-Gilbert spread constraint vector*/
     
        mpfr_t kfunc;              //Variable to capture KFunc output
        mpfr_t trapz;              //Running trapezium sum
        mpfr_t work;               //Work variable
 
        mpfr_init2(trapz,prec);
        mpfr_init2(work,prec);
        mpfr_init2(kfunc,prec);

        for (int i=0; i<t2-t1; i++){
            mpfr_init2(SConst[i],prec);
            KFunc(lmin,tau[i],kfunc,work);
            mpfr_set(trapz,kfunc,MPFR_RNDF);
            mpfr_div_d(trapz,trapz,2.0,MPFR_RNDF);

            for (int k=1; k<n; k++){
                KFunc(lmin+(lmax-lmin)*k/n,tau[i],kfunc,work);
                mpfr_add(trapz,trapz,kfunc,MPFR_RNDF);
            }//Trapezium sum loop
            KFunc(lmax,tau[i],kfunc,work);
            mpfr_set(work,kfunc,MPFR_RNDF);
            mpfr_div_d(work,work,2.0,MPFR_RNDF);
            mpfr_add(trapz,trapz,work,MPFR_RNDF);
            mpfr_mul_d(trapz,trapz,(lmax-lmin)/n,MPFR_RNDF);

            mpfr_set(SConst[i],trapz,MPFR_RNDF);
        }//SConst loop

        mpfr_clear(kfunc);
        mpfr_clear(trapz);
        mpfr_clear(work);

    }//Spread constraint

    /* Construct Dirichlet constraint vectors */
    {
    omp_set_num_threads(NCORES);
    //printf("NUM THREADS: %d\tNUM CORES: %d\n", omp_get_num_threads(),NCORES);

    #pragma omp parallel shared(G,Cov,rho,errs,widths,AvgCoeffs,KWeight)
    #pragma omp for
    for (int w=0; w<=Ns; w++){    

        //int tid = omp_get_thread_num();
        //printf("ID: %d\n",tid);
        //fflush(stdout);
        /* Constructing constraint vector */

        mpfr_t KConst[t2-t1];

        double w0 = w0s[w];

        mpfr_t work;               //Work variable
        mpfr_t temp;               //Temporary (work) variable
        mpfr_t kfunc;              //Variable to capture KFunc output
     
        mpfr_init2(work,prec);
        mpfr_init2(temp,prec);
        mpfr_init2(kfunc,prec);

        for (int i=0; i<t2-t1; i++){
            mpfr_init2(KConst[i],prec);
            KFunc(w0,tau[i],KConst[i],work);
            mpfr_mul_d(KConst[i],KConst[i],2.0,MPFR_RNDF);
        }

        mpfr_t AvgCoeff[t2-t1];
        
        if (width == 0){
            /* Calculating AvgCoeff = K^-1 x C */
            
            mpfr_t temp2;
            mpfr_init2(temp2,prec);

            for (int i=0; i<t2-t1; i++){
                mpfr_init2(AvgCoeff[i],prec);
                mpfr_set_d(AvgCoeff[i],0.0,MPFR_RNDN);
                for (int j=0; j<t2-t1; j++){
                    mpfr_mul(temp2,KInverse[i][j],KConst[j],MPFR_RNDF);
                    mpfr_add(AvgCoeff[i],AvgCoeff[i],temp2,MPFR_RNDF);
                }
            }
            
            if(g==1){
                for (int i=0; i<t2-t1; i++){
                    AvgCoeffs[w][i] = mpfr_get_ld(AvgCoeff[i],MPFR_RNDN);
                }
            }
        }else{
            /* Calculate AvgCoeff = (K^-1 x D)/(D^T K^-1 D) */
            
            mpfr_t DConst[t2-t1];
            mpfr_t KInvD[t2-t1];
            mpfr_t norm;
            mpfr_t temp2;

            mpfr_init2(temp2,prec);
            mpfr_init2(norm,prec);
            
            mpfr_set_d(norm,0.0,MPFR_RNDN);
            for (int i=0; i<t2-t1; i++){
                mpfr_init2(DConst[i],prec);
                mpfr_init2(KInvD[i],prec);

                mpfr_set(DConst[i],SConst[i],MPFR_RNDF);
            }
            for (int i=0; i<t2-t1; i++){
                mpfr_set_d(KInvD[i],0.0,MPFR_RNDN);
                for (int j=0; j<t2-t1; j++){
                    mpfr_mul(temp2,KInverse[i][j],DConst[j],MPFR_RNDF);
                    mpfr_add(KInvD[i],KInvD[i],temp2,MPFR_RNDF);
                }
                mpfr_mul(temp2,DConst[i],KInvD[i],MPFR_RNDF);
                mpfr_add(norm,norm,temp2,MPFR_RNDF);
            }//KInvC and norm loop
            for (int i=0; i<t2-t1; i++){
                mpfr_init2(AvgCoeff[i],prec);
                mpfr_set_d(AvgCoeff[i],0.0,MPFR_RNDN);
                mpfr_div(AvgCoeff[i],KInvD[i],norm,MPFR_RNDF);
            }//AvgCoeff loop

            if(g==1){
                for (int i=0; i<t2-t1; i++){
                    AvgCoeffs[w][i] = mpfr_get_ld(AvgCoeff[i],MPFR_RNDN);
                }
            }
        }

        /* Adding spectral estimate */

        double rho_est = 0.0;
        
        for (int i=0; i<t2-t1; i++){
            rho_est += mpfr_get_d(AvgCoeff[i],MPFR_RNDN)*G[i];
        }
        rho[w] = rho_est;

        widths[w] = -1.0;
        
        /* Calculate error in spectral estimate */

        double err = 0.0;
        for (int i=0; i<t2-t1; i++){
            double tmp = 0.0;
            for (int j=0; j<t2-t1; j++){
                tmp += Cov[i][j]*mpfr_get_d(AvgCoeff[j],MPFR_RNDN);
            }
            err += tmp*mpfr_get_d(AvgCoeff[i],MPFR_RNDN);
        }
        errs[w] = err;

        for (int i=0; i<t2-t1; i++){
            mpfr_clear(KConst[i]);
        }
    }//End of w loop 
    }//End of pragma 

    /* Output metadata */
    
    printf("%d;%d;%f;%f;%g;%d;%g;%g\n",Nt,Ns,wmin,wmax,alpha,prec,score,condition);

    /* Output spectral density estimate */
    
    for (int i=0; i<=Ns; i++){
        printf("%g",rho[i]);
        if (i!=Ns){
            printf(",");
        }    
    }
    printf("\n");
    
    /* Output spectral density error */
    
    for (int i=0; i<=Ns; i++){
        printf("%g",errs[i]);
        if (i!=Ns){
            printf(",");
        }    
    }
    printf("\n");
    
    /* Output resolution estimate */
    
    for (int i=0; i<=Ns; i++){
        printf("%g",widths[i]);
        if (i!=Ns){
            printf(",");
        }
    }
    printf("\n");
    
    if (g==1){

        mpfr_t kfunc2;
        mpfr_t temp3;
        mpfr_init2(kfunc2,prec);   
        mpfr_init2(temp3,prec);

        FILE *avgc;

        avgc = fopen("out.avgc","w");        

        /* Compute averaging function */
        for (int n=0; n<=Ns; n++){
            for (int i=0; i<=Ns; i++){
                long double AvgFunc = 0.0;
                for (int j=0; j<t2-t1; j++){
                    //KFunc(w0s[i],tau[j],kfunc2,temp3);
                    //AvgFunc += AvgCoeffs[n][j]*mpfr_get_ld(kfunc2,MPFR_RNDF);
                    AvgFunc += AvgCoeffs[n][j]*expl(-w0s[i]*tau[j]);
                }
                fprintf(avgf,"%g",AvgFunc);
                if(i!=Ns){
                    fprintf(avgf,"%s",",");
                }
            }
            fprintf(avgf,"%s","\n");
        }
        /* Print Averaging Coefficients */
        for (int n=0; n<=Ns; n++){
            for(int i=0; i<t2-t1; i++){
                fprintf(avgc,"%g",AvgCoeffs[n][i]);
                if(i!=t2-t1-1){
                    fprintf(avgc,"%s",",");
                }
            }
            fprintf(avgc,"%s","\n");
        }      
    }

    fprintf(fptr,"Error Mode:\t%d ([1] - Full, [0] - Partial)\n\n",emode);
    fprintf(fptr,"Scaling:\tw^%f\n\n",scaling);
    fprintf(fptr,"Start (t1):\t%d\n",t1);
    fprintf(fptr,"End (t2):\t%d\n\n",t2);

    fclose(fptr);
    fclose(avgf);

    mpfr_free_cache();

    return 0;
}

void KFunc(double w, double tau, mpfr_t kfunc, mpfr_t factor){
    /**
     * Calculates the value of the kernel function
     *
     * (double) w       - energy/mass (of which the spectrum is a function)
     * (double) t       - Euclidean time 
     * (mpfr_t) kfunc   - empty variable which will contain the value of K(w,t)
     *                    calculated by the function
     * (mpfr_t) factor  - work variable used by the function during calculation
     *
     * returns: none (result contained in kfunc)
     */

    mpfr_set_d(factor,-w*tau,MPFR_RNDF);
    mpfr_exp(factor,factor,MPFR_RNDF);      
    if (scaling == 0.5){        //Multiply by sqrt(w), zeroed at wmin (s01 channel scaling)
        mpfr_mul_d(factor,factor,sqrt(w-wmin),MPFR_RNDF);
    }
    if (scaling == 1.5){        //Multiply by sqrt(w)**3, zeroed at wmin (p10 channel scaling)
        mpfr_mul_d(factor,factor,sqrt(pow(w-wmin,3)),MPFR_RNDF);
    }
    mpfr_set(kfunc,factor,MPFR_RNDF); 
}

void WFunc(double w, double tau1, double tau2, mpfr_t wfunc, mpfr_t work){
    /**
     * Calculates the value of the width connection (the elements of the kernel
     * width matrix).
     *
     * (double) w       - energy/mass (of which the spectrum is a function)
     * (double) tau1    - Euclidean time for the first kernel function
     * (double) tau2    - Euclidean time fot the second kernel function
     * (mpfr_t) wfunc   - empty variable which will contain the value of the
     *                    width connection calculated by the function
     * (mpfr_t) work    - work variable used by the function during calculation
     *
     * returns: none (result contained in wfunc)
     */

    KFunc(w,tau1,wfunc,work);
    KFunc(w,tau2,work,work);
    mpfr_mul(wfunc,wfunc,work,MPFR_RNDF);
    //mpfr_mul_d(wfunc,wfunc,24.0,MPFR_RNDF);
}
