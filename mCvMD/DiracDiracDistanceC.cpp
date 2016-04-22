/*
	C++ Implementation of the modified Cramer--von Mises distance between two Dirac densities 
	as published in U.D. Hanebeck, "Optimal Reduction of Multivariate Dirac Mixture Densities", 
	at-Automatisierungstechnik, 2015.
	
	I compiled it against OpenBLAS (https://github.com/xianyi/OpenBLAS). However, every other BLAS
	library should work.
	
	Tested on i5-3320M, 8 GB RAM, Win10
	
	Limitations: 
		- only densities with equal numbers of components
		- no built in checks (for speed purposes)
	
	Only for academic use. Provided as is.
	
	April 21, 2016
	Maxim Dolgov, maxim.dolgov@kit.edu
*/


#include "DiracDiracDistanceC.h"

double Gamma = 0.5772156649015328606;

inline double xlog(double x) {
	if (x == 0) {
		return 0;
	}	
	else {
		return (x)*log(x);
	}
}

double DiracDiracDistanceC(double* density1, size_t numSamples1, double* mean1, double* density2, size_t numSamples2, double* mean2, size_t dim, double* bmax){
    double Df = 0;
    double Dfg = 0;
    double Dg = 0;

    double* tmp=(double*)mxMalloc(dim*sizeof(double));
    double* tmpPtr;
    
    double out;
    
    // compute Df, Dfg, Dg
    for (size_t i = 0; i < numSamples1; i++) {
        for (size_t j = 0; j < numSamples2; j++) {
            tmpPtr = density1 + j*dim;
            memcpy(tmp,density1 + i*dim,dim*sizeof(double));
            cblas_daxpy(dim, -1, tmpPtr, 1, tmp, 1);
            Df = Df + xlog(cblas_ddot(dim, tmp, 1, tmp, 1));

            tmpPtr = density2 + j*dim;
            memcpy(tmp, density1+i*dim, dim * sizeof(double));
            cblas_daxpy(dim, -1, tmpPtr, 1, tmp, 1);
            Dfg = Dfg + xlog(cblas_ddot(dim, tmp, 1, tmp, 1));

            //tmpPtr = density2[j*dim];
            memcpy(tmp, density2 + i*dim, dim * sizeof(double));
            cblas_daxpy(dim, -1, tmpPtr, 1, tmp, 1);            
            Dg = Dg + xlog(cblas_ddot(dim, tmp, 1, tmp, 1));
        }
    }
    
    // compute parts of De
    memcpy(tmp, mean1, dim * sizeof(double));
    cblas_daxpy(dim, -1, mean2, 1, tmp, 1);

    out =  pow(M_PI, (double)dim / 2)*((Df - 2 * Dfg + Dg) / (numSamples1*numSamples2) + 2 * (log(4*(*bmax)*(*bmax))-Gamma)*cblas_ddot(dim, tmp, 1, tmp, 1))/8;
    
    mxFree(tmp);
    return out;
}