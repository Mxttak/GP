/*
	Computes the cross distance matrix between the first Dirac densitiy and the other densities. 
	The distances are as published in U.D. Hanebeck, "Optimal Reduction of Multivariate Dirac 
	Mixture Densities", at-Automatisierungstechnik, 2015.
	
	I compiled it against OpenBLAS (https://github.com/xianyi/OpenBLAS). However, every other BLAS
	library should work.
	
	Tested on i5-3320M, 8 GB RAM, Win10
	
	Limitations: 
		- only densities with equal numbers of components
		- no built in checks (for speed purposes); this might crash your Matlab if you pass something wrong
	
	Structure of the input parameters: see Matlab script with the example
	
	Only for academic use. Provided as is.
	
	April 21, 2016
	Maxim Dolgov, maxim.dolgov@kit.edu
*/

#include <mex.h>
#include "DiracDiracDistanceC.h"

void mexFunction(int nlhs, mxArray* plhs[],			// output arguments
				 int nrhs, const mxArray* prhs[]) {	// input arguments

	// get inputs
		// prhs[0]: density
		// prhs[1]: mean
		// prhs[2]: densities
		// prhs[3]: means
		// prhs[4]: bmax

	double* density = mxGetPr(prhs[0]);
	double* mean = mxGetPr(prhs[1]);

	mwSize numDims = mxGetNumberOfDimensions(prhs[2]);
	mwSize* dim = (mwSize*)mxGetDimensions(prhs[2]);
	double* densities = mxGetPr(prhs[2]);
	double* means = mxGetPr(prhs[3]);

	double* bmax = mxGetPr(prhs[4]);

	// allocate output
	plhs[0] = mxCreateDoubleMatrix(1, dim[2], mxREAL);
	double* out = mxGetPr(plhs[0]);
	
	// compute distances
	for (mwSize i = 0; i < dim[2]; i++) {
		out[i] = DiracDiracDistanceC(density,dim[1], mean,densities + dim[0]*dim[1]*i, dim[1], means + dim[0]*i, dim[0], bmax);
	}
	return;
}