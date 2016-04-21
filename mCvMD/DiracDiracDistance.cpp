/*
	Implementation of the modified Cramer--von Mises distance between two Dirac densities 
	as published in U.D. Hanebeck, "Optimal Reduction of Multivariate Dirac Mixture Densities", 
	at-Automatisierungstechnik, 2015. Much much faster than the implementation as a Matlab class.
	
	I compiled it against OpenBLAS (https://github.com/xianyi/OpenBLAS). However, every other BLAS
	library should work.
	
	Tested on i5-3320M, 8 GB RAM, Win10, Matlab 2015b x64
	
	Limitations: 
		- only densities with equal numbers of components
		- no built in checks (for speed purposes)
	
	Structure of the input parameters: see Matlab script with the example
	
	Only for academic use. Provided as is.
	
	April 21, 2016
	Maxim Dolgov, maxim.dolgov@kit.edu
*/

#include <mex.h>
#include "OpenBLAS-v0.2.12-Win64-int64/include/cblas.h"
#define _USE_MATH_DEFINES
#include <math.h>
#include <memory.h>

double Gamma = 0.5772156649015328606;

inline double xlog(double x) {
	if (x == 0) {
		return 0;
	}	
	else {
		return (x)*log(x);
	}
}

void mexFunction(int nlhs, mxArray* plhs[],			// output arguments
				int nrhs, const mxArray* prhs[]) {	// input arguments

	// get inputs
		// density1
		size_t dim1 = mxGetM(prhs[0]);
		size_t numSamples1 = mxGetN(prhs[0]);
		double* density1 = mxGetPr(prhs[0]);
		
		// mean1
		double* mean1 = mxGetPr(prhs[1]);

		// density2
		size_t dim2 = mxGetM(prhs[2]);
		size_t numSamples2 = mxGetN(prhs[2]);
		double* density2 = mxGetPr(prhs[2]);
		
		// mean2
		double* mean2 = mxGetPr(prhs[3]);

		double* bmax = mxGetPr(prhs[4]);

	// allocate output
		plhs[0] = mxCreateDoubleMatrix(1,1, mxREAL);
		double* out = mxGetPr(plhs[0]);
		
		double Df = 0;
		double Dfg = 0;
		double Dg = 0;

		double* tmp=(double*)malloc(dim1*sizeof(double));
		double* tmpPtr;
		
		// compute Df, Dfg, Dg
		for (size_t i = 0; i < numSamples1; i++) {
			for (size_t j = 0; j < numSamples2; j++) {
				tmpPtr = density1 + j*dim2;
				memcpy(tmp,density1+i*dim1,dim1*sizeof(double));
				cblas_daxpy(dim1, -1, tmpPtr, 1, tmp, 1);
				Df = Df + xlog(cblas_ddot(dim1, tmp, 1, tmp, 1));

				tmpPtr = density2 + j*dim2;
				memcpy(tmp, density1 + i*dim1, dim1 * sizeof(double));
				cblas_daxpy(dim1, -1, tmpPtr, 1, tmp, 1);
				Dfg = Dfg + xlog(cblas_ddot(dim1, tmp, 1, tmp, 1));

				tmpPtr = density2 + j*dim2;
				memcpy(tmp, density2 + i*dim1, dim1 * sizeof(double));
				cblas_daxpy(dim1, -1, tmpPtr, 1, tmp, 1);
				
				Dg = Dg + xlog(cblas_ddot(dim1, tmp, 1, tmp, 1));
			}
		}
		
		// compute parts of De
		memcpy(tmp, mean1, dim1 * sizeof(double));
		cblas_daxpy(dim1, -1, mean2, 1, tmp, 1);

		*out = pow(M_PI, (double)dim1 / 2)*((Df - 2 * Dfg + Dg) / (numSamples1*numSamples2) + 2 * (log(4*(*bmax)*(*bmax))-Gamma)*cblas_ddot(dim1, tmp, 1, tmp, 1))/8;
	return;
}