/*
    Implemenetation of the Wasserstein distance between Dirac densities. See J.R. Hoffman
    and R.P.S. Mahler, "Multitarget Miss Distance and its Applications", FUSION, 2002.
    
    This implementation uses Eigen available at https://eigen.tuxfamily.org and LAPJV (Linear Assignment
    Problem) algorithm by R. Jonker and A. Volgenant (see R. Jonker and A. Volgenant "A 
    Shortest Augmenting Path Algorithm for Dense and Sparse Linear Assignment Problems," 
    Computing 38, 325-340, 1987)
        
    Build (on Windows): 
        mex WassersteinDistance.cpp lapjv/lap.cpp -IEigen
    
        tested with Matlab R2015b on Win10 x64
        
    Requirements:
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
    May 13, 2016
    
    Provided as is. No warranty, no commercial usage.
*/

#include <mex.h>
#include <Eigen>
#include <cmath>
#include "lapjv/lap.h"

/*------------------------------------------------------------------------------------------------------*/

double DistMat(Eigen::Map<Eigen::MatrixXd> &dens1, Eigen::Map<Eigen::MatrixXd> &dens2,double **dist, int *rowsol, int *colsol, double *u, double *v);

/*------------------------------------------------------------------------------------------------------*/

void mexFunction(int nOnputs, mxArray* outputs[],			// output arguments
				int nInputs, const mxArray* inputs[]) {	// input arguments
    
    // get inputs
    // first density
    double* densities1 = mxGetPr(inputs[0]);
        mwSize numDimDensity1 = mxGetNumberOfDimensions(inputs[0]);
        mwSize* dimsDensity1 = (mwSize*)mxGetDimensions(inputs[0]);
    
    // assign dimensions
    mwSize m = dimsDensity1[0], n = dimsDensity1[1];
    
    mwSize nDensities1;
    if(numDimDensity1 == 2) 
        nDensities1 = 1;
    else 
        nDensities1 = dimsDensity1[2];   
        
    mwSize nDensities2 = nDensities1; 
    
    // second density and mean
    mwSize *dimsDensity2, *dimsMeans2, numDimDensity2;
    double *means2, *densities2;
    
    if(nInputs == 2){   
        densities2 = mxGetPr(inputs[1]);
            mwSize numDimDensity2 = mxGetNumberOfDimensions(inputs[1]);
            dimsDensity2 = (mwSize*) mxGetDimensions(inputs[1]);
        means2 = mxGetPr(inputs[4]);
                
        if(numDimDensity2 == 2)
            nDensities2 = 1;
        else
            nDensities2 = dimsDensity2[2];
    }
    
    // create outputs 
    outputs[0] = mxCreateDoubleMatrix(nDensities1,nDensities2, mxREAL);
    double* out = mxGetPr(outputs[0]);
    
    // auxiliary variables
    double **dist = (double**) mxMalloc(n*sizeof(double**));
    for(unsigned int i=0; i<n; i++)
        dist[i] = (double*) mxCalloc(n,sizeof(double));
    int *rowsol = (int*) mxMalloc(n*sizeof(int));
    int *colsol = (int*) mxMalloc(n*sizeof(int));
    double* u = (double*) mxMalloc(n*sizeof(double));
    double* v = (double*) mxMalloc(n*sizeof(double));
    
    // auxiliary variables with dummy initialization (this is necessary, see Eigen docs)
    Eigen::Map<Eigen::MatrixXd> densMat1(densities1,m,n), densMat2(densities1,m,n);

    // compute distances
    if(nInputs == 1){
        for(unsigned int i=0; i<nDensities1; i++){
            for(unsigned int j=i+1; j<nDensities1; j++){
                new (&densMat1) Eigen::Map<Eigen::MatrixXd>(densities1+i*m*n,m,n);
                new (&densMat2) Eigen::Map<Eigen::MatrixXd>(densities1+j*m*n,m,n);
                                
                out[i*nDensities1+j] = sqrt(DistMat(densMat1,densMat2,dist,rowsol,colsol,u,v)/n);
                out[j*nDensities1+i] = out[i*nDensities1+j];
            }
        }
    } else{
        for(unsigned int i=0; i<nDensities1; i++){
            for(unsigned int j=i+1; j<nDensities1; j++){
                new (&densMat1) Eigen::Map<Eigen::MatrixXd>(densities1+i*m*n,m,n);
                new (&densMat2) Eigen::Map<Eigen::MatrixXd>(densities2+j*m*n,m,n);
                
                out[j*nDensities1+i] = sqrt(DistMat(densMat1,densMat2,dist,rowsol,colsol,u,v)/n);
            }
        }
    }
    
    mxFree(dist); mxFree(rowsol); mxFree(colsol); mxFree(u); mxFree(v); 
    return;
}

/*------------------------------------------------------------------------------------------------------*/

double DistMat(Eigen::Map<Eigen::MatrixXd> &dens1, Eigen::Map<Eigen::MatrixXd> &dens2,double **dist, int *rowsol, int *colsol, double *u, double *v){
    
    unsigned int n = dens1.cols(), m = dens1.rows();
    
    Eigen::MatrixXd tmp;
    Eigen::VectorXd tmpVec;

    for(unsigned int i=0; i<n; i++){
        for(unsigned int j=0;j<n; j++){
            tmpVec = dens1.col(i) - dens2.col(j);
            dist[i][j] = tmpVec.squaredNorm();
        }
    }
    
    return lap(n,dist,rowsol,colsol,u,v);
}