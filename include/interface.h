#ifndef INTERFACE_H
#define INTERFACE_H

#include <gmp.h>
#include <mpfr.h>

#ifdef __cplusplus
extern "C" {
#endif

void c_zcall(int prec,int width,mpfr_t *A, mpfr_t *AInv, mpfr_t *D);

#ifdef __cplusplus
}
#endif

#endif
