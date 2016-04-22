#ifndef DIRACDIRACDISTANCEC_H
#define DIRACDIRACDISTANCEC_H

#include "OpenBLAS-v0.2.12-Win64-int64/include/cblas.h"
#define _USE_MATH_DEFINES
#include <math.h>
#include <memory.h>
#include <mex.h>

    double DiracDiracDistanceC(double* density1, size_t numSamples1, double* mean1,double* density2,size_t numSamples2, double* mean2,size_t dim,double* bmax);

#endif