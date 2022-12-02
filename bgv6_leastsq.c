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

    int Nt;                  //Number of time points in lattice (the temperature)
    int Ns;                  //Number of sample slices to calculate omega
    int t1;                  //Initial Euclidean time (inclusive)
    int t2;                  //Final Euclidean time (exclusive)

    FILE *fptr;
    fptr = fopen("out.report","w");
    fprintf(fptr,"%s","OUTPUT REPORT\n====================\n\n");

    // Check validity of cmd line args
    if (argc == 1){
        printf("Argument Error: Must supply at least 1 command line argument! (debug mode)\n");
        return(22);
    }
    if (argc > 7){
        printf("Argument Error: Too many arguments!\n");
        return(7); 
    }

    int g = atoi(argv[1]);       //Debug mode for outputting averaging coefficients
    int emode = 1;               //Full error mode (1) or partial error mode (0) -- Full errors uses covariance matrix

    scanf("%d",&Nt);
    scanf("%d",&Ns);
    
    /*
    t1 = 1;                      //Default start time (inclusive)
    t2 = 64;                     //Default end time (t+1, exclusive)
    */

    double alpha;                //Whitening parameter
    double condition;            //Condition of the weight matrix

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

        //Check if t2 set to default to max
        if (t2 == -1){
            t2 = Nt;
        }
    }

    scaling = 0.5;               //Kernel rescaling - 0.5 for S-wave and 1.5 for P-wave

    double lmin = wmin;          //Start of integration range
    double lmax = 10.0;          //End of integration range
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

    //Read diagonal elements    
    for (int i=0; i<Nt; i++){
        double tmp;
        scanf("%lf",&tmp);
        if (i >= t1 && i < t2){
            Cov[i-t1][i-t1] = tmp;
        }
    }

    //Construct off-diagonal components
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

    for (int i=0; i<t2-t1; i++){       //Diagonal, equal-time elements
        
        mpfr_init2(KWeight[i][i],prec);
        WFunc(lmin,tau[i],tau[i],wfunc,work);
        mpfr_set(trapz,wfunc,MPFR_RNDF);
        mpfr_div_d(trapz,trapz,2.0,MPFR_RNDF);

        for (int k=1; k<n; k++){
            WFunc(lmin+(lmax-lmin)*k/n,tau[i],tau[i],wfunc,work),
            mpfr_add(trapz,trapz,wfunc,MPFR_RNDF);
        }
        WFunc(lmax,tau[i],tau[i],wfunc,work);
        mpfr_set(work,wfunc,MPFR_RNDF);
        mpfr_div_d(work,work,2.0,MPFR_RNDF);
        mpfr_add(trapz,trapz,work,MPFR_RNDF);

        mpfr_mul_d(trapz,trapz,(lmax-lmin)/n,MPFR_RNDF);
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

    mpfr_t S[t2-t1];           //Vector of singular elements

    for(int i=0; i<t2-t1; i++){
        mpfr_init2(S[i],prec);
    }

    /* Start of C++ block */
    c_zcall(prec,t2-t1,*KCopy,*KInverse,S); //Pipe matrices to ZKCM interface
    /* End of C++ block */

    condition = fabsl(mpfr_get_ld(S[0],MPFR_RNDN)/mpfr_get_ld(S[t2-t1-1],MPFR_RNDN));         //Save condition number

    mpfr_clear(temp);
    mpfr_clear(trapz);
    mpfr_clear(work);
    mpfr_clear(wfunc);

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

        /* Calculate spectral estimate */

        double rho_est = 0.0;
        
        for (int i=0; i<t2-t1; i++){
            rho_est += mpfr_get_d(AvgCoeff[i],MPFR_RNDN)*G[i];
        }

        rho[w] = rho_est;

        widths[w] = 1/w0;   //Oldenburg estimate of least-squares width
        
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

        //Explicitly flush KConst for next loop
        for (int i=0; i<t2-t1; i++){
            mpfr_clear(KConst[i]);
        }

    }//End of w loop 
    }//End of pragma wrapper

    /*
     * =======================================================================
                            End of Inversion Container
                                    
     * =======================================================================
    */

    /* Output metadata */
    
    printf("%d;%d;%f;%f;%g;%d;%g\n",Nt,Ns,wmin,wmax,alpha,prec,condition);

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
   
    /* Note that this resolution output uses the reciprocal estimate [Oldenburg]
       and as such is not as accurate as the resolution estimate for the Backus-
       Gilbert spread function.
    */ 
    for (int i=0; i<=Ns; i++){
        printf("%g",widths[i]);
        if (i!=Ns){
            printf(",");
        }
    }
    printf("\n");

    /* Output averaging coefficients */

    if (g==1){

        FILE *avgc;

        avgc = fopen("out.avgc","w");        

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
    
    /* Print Inverse of Kernel Weight Matrix */
    if (g==2){

        FILE *kinv;

        kinv = fopen("out.kinv","w");        

        /* Print Averaging Coefficients */
        for (int i=0; i<t2-t1; i++){
            for(int j=0; j<t2-t1; j++){
                fprintf(kinv,"%0.16Lf",mpfr_get_ld(KInverse[i][j],MPFR_RNDN));
                if(j<t2-t1-1){
                    fprintf(kinv,"%s",",");
                }
            }
            fprintf(kinv,"%s","\n");
        }      
    }
    /* Output miscellaneous info to report file */

    fprintf(fptr,"Error Mode:\t%d ([1] - Full, [0] - Partial)\n\n",emode);
    fprintf(fptr,"Scaling:\tw^%f\n\n",scaling);
    fprintf(fptr,"Start (t1):\t%d\n",t1);
    fprintf(fptr,"End (t2):\t%d\n\n",t2);

    fclose(fptr);

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
