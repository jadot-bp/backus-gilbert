#include <cstdlib>
#include <gmp.h>
#include <mpfr.h>

#include "interface.h"
#include "zcall.h"

#ifdef __cplusplus
extern "C" {
#endif

static ZKCM *ZKCM_instance = NULL;

void initZKCM() {
    if (ZKCM_instance == NULL){
        ZKCM_instance = new ZKCM();
    }
}

void c_zcall(int prec,int width,mpfr_t *A, mpfr_t *AInv, mpfr_t *D){
    initZKCM();
    ZKCM_instance->zcall(prec,width,A,AInv,D);
}

#ifdef __cplusplus
}
#endif
