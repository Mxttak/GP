/*
	Computes the cross distance matrix between all supplied Dirac densities. The distances are 
	as published in U.D. Hanebeck, "Optimal Reduction of Multivariate Dirac Mixture Densities", 
	at-Automatisierungstechnik, 2015.
	
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
        // prhs[0]: densities
        // prhs[1]: means
        // prhs[2]: bmax    
	mwSize numDims = mxGetNumberOfDimensions(prhs[0]);
    mwSize* dim = (mwSize*) mxGetDimensions(prhs[0]);
    double* densities = mxGetPr(prhs[0]);
    double* means = mxGetPr(prhs[1]);
    
    double* bmax = mxGetPr(prhs[2]);
    
    // allocate output
    plhs[0] = mxCreateDoubleMatrix(dim[2],dim[2],mxREAL);
    double* out = mxGetPr(plhs[0]);
    
    for(mwSize i = 0; i < dim[2]; i++){
        for(mwSize j = i+1; j < dim[2]; j++){     
            out[j+dim[2]*i] = DiracDiracDistanceC(densities + dim[0]*dim[1]*i, dim[1], means + dim[0]*i, densities + dim[0]*dim[1]*j, dim[1], means + dim[0]*j, dim[0], bmax);
            out[i+dim[2]*j] = out[j+dim[2]*i];
        }
    }
    
    return;    
}