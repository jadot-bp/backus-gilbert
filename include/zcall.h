#ifndef ZCALL_H
#define ZCALL_H

#include <zkcm.hpp>

class ZKCM{
    public:
        ZKCM();
        int zcall(int prec,int width,mpfr_t *A, mpfr_t *AInv, mpfr_t *D);
};
#endif
