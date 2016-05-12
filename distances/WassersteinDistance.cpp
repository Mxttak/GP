/*
    Implemenetation of the Wasserstein distance between Dirac densities. See J.R. Hoffman
    and R.P.S. Mahler, "Multitarget Miss Distance and its Applications", FUSION, 2002.
    
    This implementation uses CBLAS from http://www.netlib.org/ and LAPJV (Linear Assignment
    Problem) algorithm by R. Jonker and A. Volgenant (see R. Jonker and A. Volgenant "A 
    Shortest Augmenting Path Algorithm for Dense and Sparse Linear Assignment Problems," 
    Computing 38, 325-340, 1987)
        
    Build (on Windows): 
        mex WassersteinDistance.cpp lapjv/lap.cpp libcblas.dll.a
    
        tested with Matlab R2015b on Win10 x64
        
    Requirements:
        You will need CBLAS that is available on http://www.netlib.org/ 
        You can get info on how to build CBLAS for Windows at http://icl.cs.utk.edu/lapack-for-windows/lapack/
        A prebuilt version of CBLAS for Win x64 is supplied with this implementation. However,
        you will need MinGWx64 in your path. 
        Yout will need the implementation of the LAPJV algorithm. See, e.g., my other
        repository https://github.com/Mxttak/lapjv_matlab
        
    Usage:
        [dist] = WassersteinDistance(dens1)
        [dist] = WassersteinDistance(dens1,dens2)
        
        - dens1: set of Dirac densities arranged in a 3D array as follows: dens1(:,i,j): i-th sample of j-th density        
        - dens2: second set of Dirac densities arranged in a 3D array as follows: dens2(:,i,j): i-th sample of j-th density
        - dist:  matrix with distances between the densities from dens1 or dens1 and dens2, respectively
        
    Warning: no builtin checks 
        
    Maxim Dolgov
    May 12, 2016
    
    Provided as is. No warranty, no commercial usage.
*/

#include <mex.h>
#include <memory.h>
#include "lapjv/lap.h"
#include "include/cblas.h"
#include <math.h>

/*-------------------------------------------------------------------------------------------*/

double DistMat(unsigned int m, unsigned int n, double* dens1, double* dens2, double** dist,double *tmp,int* rowsol,int *colsol, double* u, double* v);
double DistMat(unsigned int m, unsigned int n, double* dens1, double** dist,double *tmp,int* rowsol,int *colsol, double* u, double* v);

/*-------------------------------------------------------------------------------------------*/


void mexFunction(int nlhs, mxArray* plhs[],			// output arguments
				int nrhs, const mxArray* prhs[]) {	// input arguments
    
    // get inputs
    // first density
    double* dens1 = mxGetPr(prhs[0]);
    mwSize* dims1 = (mwSize*) mxGetDimensions(prhs[0]);
    mwSize numDims1 = mxGetNumberOfDimensions(prhs[0]);

    mwSize m = *dims1, n = dims1[1], numDensities1,numDensities2;
    double* dens2;
    
    if(numDims1 == 2)
        numDensities1 = 1;
    else
        numDensities1 = dims1[2];
    
    // eventually the second density
    if(nrhs == 1){
        numDensities2 = numDensities1;
        dens2 = dens1;
    } else{
        dens2 = mxGetPr(prhs[1]);
        mwSize* dims2 = (mwSize*) mxGetDimensions(prhs[1]);
        mwSize numDims2 = mxGetNumberOfDimensions(prhs[1]);

        if(numDims2 == 2)
            numDensities2 = 1;
        else
            numDensities2 = dims2[2];
    }

    // create outputs 
    plhs[0] = mxCreateDoubleMatrix(numDensities1,numDensities2, mxREAL);
    double* out = mxGetPr(plhs[0]);
    
    double *tmp = (double*) mxMalloc(m*sizeof(double));
    double **dist = (double**) mxMalloc(n*sizeof(double**));
    for(unsigned int i=0; i<n; i++)
        dist[i] = (double*) mxMalloc(n*sizeof(double));
    int *rowsol = (int*) mxMalloc(n*sizeof(int));
    int *colsol = (int*) mxMalloc(n*sizeof(int));
    double* u = (double*) mxMalloc(n*sizeof(double));
    double* v = (double*) mxMalloc(n*sizeof(double));
    
    if(nrhs == 1){
        // compute distance for one density set
        for(unsigned int i=0; i<numDensities1; i++){
            for(unsigned int j=i+1; j<numDensities1; j++){
                out[i*numDensities1+j] = sqrt(DistMat(m,n,dens1+i*m*n,dens1+j*m*n,dist,tmp,rowsol,colsol,u,v)/n);
                out[j*numDensities1+i] = out[i*numDensities1+j];
            }
        }
    } else{
        // compute distance for two density sets
        for(unsigned int i=0; i<numDensities1; i++){
            for(unsigned int j=0; j<numDensities2; j++){
                out[i+j*numDensities1] = sqrt(DistMat(m,n,dens1+i*m*n,dens2+j*m*n,dist,tmp,rowsol,colsol,u,v)/n);
            }
        }
    }
     
    return;
}

/*-------------------------------------------------------------------------------------------*/

double DistMat(unsigned int m, unsigned int n, double* dens1, double* dens2, double** dist,double *tmp,int* rowsol,int *colsol, double* u, double* v){
    
    for(unsigned int i=0; i<n; i++){
        for(unsigned int j=0; j<n; j++){
            memcpy(tmp,dens1+i*m,m*sizeof(double));
            cblas_daxpy(m,-1,dens2+j*m,1,tmp,1);
            dist[i][j] = cblas_ddot(m,tmp,1,tmp,1);
        }
    }
    
    return lap(n,dist,rowsol,colsol,u,v);
}

double DistMat(unsigned int m, unsigned int n, double* dens1, double** dist,double *tmp,int* rowsol,int *colsol, double* u, double* v){
    for(unsigned int i=0; i<n; i++){
        for(unsigned int j=i+1; j<n; j++){
            memcpy(tmp,dens1+i*m,m*sizeof(double));
            cblas_daxpy(m,-1,dens1+j*m,1,tmp,1);
            dist[i][j] = cblas_ddot(m,tmp,1,tmp,1);
            dist[j][i] = dist[i][j]; 
        }
    }
    return lap(n,dist,rowsol,colsol,u,v);
}

/*-------------------------------------------------------------------------------------------*/