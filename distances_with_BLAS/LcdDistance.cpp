/*
    Implementation of the modified Cramer - von Mises distance for Dirac densities that
    was published in U.D. Hanebeck, "Optimal Reduction of Multivariate Dirac Mixture 
    Densities", at-Automatisierungstechnik, 2015.
    
    This implementation uses CBLAS from http://www.netlib.org/ 
    
    Build (on Windows): 
        mex LcdDistance.cpp libcblas.dll.a
    
        tested with Matlab R2015b on Win10 x64
        
    Requirements:
        You will need CBLAS that is available on http://www.netlib.org/ 
        You can get info on how to build CBLAS for Windows at http://icl.cs.utk.edu/lapack-for-windows/lapack/
        A prebuilt version of CBLAS for Win x64 is supplied with this implementation. However,
        you will need MinGWx64 in your path.
    
    Limitations:
        No builtin checks.
        Only densities with equal numbers of samples.
    
    Usage:
        [dist] = LcdDistance(bmax,dens1,mean1)
        [dist] = LcdDistance(bmax,dens1,mean1,dens2,mean2)
        
        - bmax:  distance parameter (see the paper)
        - dens1: set of Dirac densities arranged in a 3D array as follows: dens1(:,i,j): i-th sample of j-th density
        - mean1: means of the densities in dens1 computed, e.g., as follows mean(dens1,2)
        - dens2: second set of Dirac densities arranged in a 3D array as follows: dens2(:,i,j): i-th sample of j-th density
        - mean2: means of the densities in dens2 computed, e.g., as follows mean(dens2,2)
        - dist:  matrix with distances between the densities from dens1 or dens1 and dens2, respectively
        
    Maxim Dolgov
    May 12, 2016
    
    Provided as is. No warranty, no commercial usage.
*/

#define _USE_MATH_DEFINES
#include <math.h>
#include <memory.h>
#include "include/cblas.h"
#include <mex.h>

double Gamma = 0.5772156649015328606;
double one = 1;
double mone = -1;
ptrdiff_t stepone = 1;

/*-------------------------------------------------------------------------------------------*/

double DiracDiracDistanceC(double* density1, size_t numSamples1, double* mean1, double* density2, size_t numSamples2, double* mean2, size_t dim, double* bmax);

inline double xlog(double x) {
	if (x == 0) {
		return 0;
	}	
	else {
		return (x)*log(x);
	}
}

/*-------------------------------------------------------------------------------------------*/

/*            
inputs: prhs[0]: bmax
        prhs[1]: densities1 
        prhs[2]: means1
        prhs[3]: densities2
        prhs[4]: means2        
*/

void mexFunction(int nlhs, mxArray* plhs[],			// output arguments
				int nrhs, const mxArray* prhs[]) {	// input arguments
    
    // get inputs
    // bmax    
    double* bmax = mxGetPr(prhs[0]);
    
    // first density
    double* densities1 = mxGetPr(prhs[1]);
        mwSize numDimDensity1 = mxGetNumberOfDimensions(prhs[1]);
        mwSize* dimsDensity1 = (mwSize*)mxGetDimensions(prhs[1]);
    
    // first mean
    double* means1 = mxGetPr(prhs[2]);
    
    // assign dimensions
    mwSize m = dimsDensity1[0], nSamples1 = dimsDensity1[1];
    
    mwSize nDensities1;
    if(numDimDensity1 == 2) 
        nDensities1 = 1;
    else 
        nDensities1 = dimsDensity1[2];
        
    mwSize nDensities2 = nDensities1, nSamples2 = nSamples1;

    // second density and mean
    mwSize *dimsDensity2, *dimsMeans2, numDimDensity2;
    double *means2, *densities2;
    
    if(nrhs == 5){   
        densities2 = mxGetPr(prhs[3]);
            mwSize numDimDensity2 = mxGetNumberOfDimensions(prhs[1]);
            dimsDensity2 = (mwSize*) mxGetDimensions(prhs[3]);
        means2 = mxGetPr(prhs[4]);//               dimsMeans2 = (mwSize*) mxGetDimensions(prhs[4]);
        
        nSamples2 = dimsDensity2[1];// nDensities2 = dimsMeans2[1];
        
        if(numDimDensity2 == 2)
            nDensities2 = 1;
        else
            nDensities2 = dimsDensity2[2];
    }
    
               
    // create outputs 
    plhs[0] = mxCreateDoubleMatrix(nDensities1,nDensities2, mxREAL);
    double* out = mxGetPr(plhs[0]);
    
    // auxiliary variables
    double* tmp = (double*) mxMalloc(m*sizeof(double));
    
    // compute the distances
    if(nrhs == 5){
        // if two densities are provided
        for(mwSize i=0;i<nDensities1;i++){
            for(mwSize j=0;j<nDensities2;j++){
                out[i+j*nDensities1] = DiracDiracDistanceC(densities1+i*m*nSamples1, nSamples1, means1+i*m, densities2+j*m*nSamples2, nSamples2, means2+j*m, m,bmax);
            }
        }
    } else{
        // if one density is provided
        for(mwSize j=0;j<nDensities2;j++){
            for(mwSize i=j+1;i<nDensities1;i++){                                   
                out[j+i*nDensities1] = DiracDiracDistanceC(densities1+i*m*nSamples1, nSamples1, means1+i*m, densities1+j*m*nSamples1, nSamples1, means1+j*m, m,bmax);
                out[i+j*nDensities1] = out[j+i*nDensities1];
            }
        }
    }
    
    mxFree(tmp);
    return;
}

/*-------------------------------------------------------------------------------------------*/

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
            
            //tmpPtr = density2+j*dim;
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