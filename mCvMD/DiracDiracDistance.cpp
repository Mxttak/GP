/*
	Implementation of the modified Cramer--von Mises distance between two Dirac densities 
	as published in U.D. Hanebeck, "Optimal Reduction of Multivariate Dirac Mixture Densities", 
	at-Automatisierungstechnik, 2015. Much much faster than the implementation as a Matlab class.
	
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
		// prhs[0]: density1
		// prhs[1]: mean1
		// prhs[2]: density2
		// prhs[3]: mean2
		// prhs[4]: bmax
		
		mwSize numDims = mxGetNumberOfDimensions(prhs[2]);
		mwSize* dim = (mwSize*)mxGetDimensions(prhs[2]);
		double* density1 = mxGetPr(prhs[0]);

		double* mean1 = mxGetPr(prhs[1]);

		double* density2 = mxGetPr(prhs[2]);

		double* mean2 = mxGetPr(prhs[3]);

		double* bmax = mxGetPr(prhs[4]);

	// allocate output
		plhs[0] = mxCreateDoubleMatrix(1,1, mxREAL);
		double* out = mxGetPr(plhs[0]);
				
		*out = DiracDiracDistanceC(density1,dim[1],mean1,density2,dim[1],mean2,dim[0],bmax);
	return;
}