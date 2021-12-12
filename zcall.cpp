#include <iostream>
#include <string>
#include <math.h>
#include "zcall.h"
#include <gmp.h>
#include <mpfr.h>

ZKCM::ZKCM() {
}
int ZKCM::zcall(int prec,int width, mpfr_t *A, mpfr_t *AInv, mpfr_t *D)
{
    zkcm_set_default_prec(prec);        //Set zkcm precision to that specified by bgv*

    double tol = pow(2,-prec);  //Moore-Penrose tolerance

    zkcm_matrix zA(width);      //Initialise zkcm matrix to read in A
    zkcm_matrix zAInv(width);

    for (int i=0; i<width; i++){
        for (int j=0; j<width;j++){
            zA(i,j) = A[i+j*width];     //Load mpfr type into zkcm matrix
        }
    }

    zA = re(zA);

    zkcm_matrix zU(width);
    zkcm_matrix zV(width);
    zkcm_matrix zD(width);
    
    int result; 
    result = SVD(zD,zU,zV,zA);  //Perform SVD on A (zA)

    zkcm_matrix zUh(width);     //Conjugate transpose

    zUh = adjoint(zU);

    for (int i=0; i<width; i++){
        for (int j=0; j<width; j++){
            zAInv(i,j) = 0.0;
            for (int k=0; k<width; k++){
                if (getRinD(zD(k,k))>tol){
                    zAInv(i,j) += zV(i,k)*(1.0/zD(k,k))*zUh(k,j);
                }else{
                    zAInv(i,j) += zV(i,k)*zD(k,k)*zUh(k,j);
                }
            }
        }
    }

    mpf_class mpf_element;              //Create mpf_class to accept zkcm type conversion
    mpfr_t mpfr_element;                //Create mpfr_t placeholder for mpf to mpfr conversion

    /* Load SVD output into mpfr matrices */
    for (int i=0; i<width; i++){
        for(int j=0; j<width; j++){
            mpf_element = getRinmpf_class(zAInv(i,j));
            mpfr_set_f(AInv[i+j*width],mpf_element.get_mpf_t(),MPFR_RNDF);
        }
    }
    /* Save singular values */
    for (int i=0; i<width; i++){
        mpf_element = getRinmpf_class(zD(i,i));
        mpfr_set_f(D[i],mpf_element.get_mpf_t(),MPFR_RNDF);
    }
    //Free miscellaneous internal memory allocations.
    return 0;
}
