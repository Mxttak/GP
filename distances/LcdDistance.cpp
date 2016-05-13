/*
    Implementation of the modified Cramer - von Mises distance for Dirac densities that
    was published in U.D. Hanebeck, "Optimal Reduction of Multivariate Dirac Mixture 
    Densities", at-Automatisierungstechnik, 2015.
    
    This implementation uses Eigen available at https://eigen.tuxfamily.org
    
    Build (on Windows): 
        mex LcdDistance.cpp -IEigen
    
        tested with Matlab R2015b on Win10 x64
        
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
    May 13, 2016
    
    Provided as is. No warranty, no commercial usage.
*/

#include <mex.h>
#include <cmath>
#include <Eigen>

/*------------------------------------------------------------------------------------------------------*/

double Gamma = 0.5772156649015328606;

/*------------------------------------------------------------------------------------------------------*/

double DiracDiracDistanceC(Eigen::Map<Eigen::MatrixXd> &densMat1, Eigen::Map<Eigen::VectorXd> &densMean1, Eigen::Map<Eigen::MatrixXd> &densMat2, Eigen::Map<Eigen::VectorXd> &densMean2, double *bmax);

inline double xlog(double x) {
	if (x == 0) {
		return 0;
	}	
	else {
		return (x)*log(x);
	}
}
    
/*------------------------------------------------------------------------------------------------------*/

void mexFunction(int nOnputs, mxArray* onputs[],			// output arguments
				int nInputs, const mxArray* inputs[]) {	// input arguments
    
    // get inputs
    // bmax    
    double* bmax = mxGetPr(inputs[0]);
    
    // first density
    double* densities1 = mxGetPr(inputs[1]);
        mwSize numDimDensity1 = mxGetNumberOfDimensions(inputs[1]);
        mwSize* dimsDensity1 = (mwSize*)mxGetDimensions(inputs[1]);
    
    // first mean
    double* means1 = mxGetPr(inputs[2]);
    
    // assign dimensions
    mwSize m = dimsDensity1[0], nSamples1 = dimsDensity1[1];
    
    mwSize nDensities1;
    if(numDimDensity1 == 2) 
        nDensities1 = 1;
    else 
        nDensities1 = dimsDensity1[2];
        
    mwSize nDensities2 = nDensities1;

    // second density and mean
    mwSize *dimsDensity2, *dimsMeans2, numDimDensity2;
    double *means2, *densities2;
    
    if(nInputs == 5){   
        densities2 = mxGetPr(inputs[3]);
            mwSize numDimDensity2 = mxGetNumberOfDimensions(inputs[3]);
            dimsDensity2 = (mwSize*) mxGetDimensions(inputs[3]);
        means2 = mxGetPr(inputs[4]);
                
        if(numDimDensity2 == 2)
            nDensities2 = 1;
        else
            nDensities2 = dimsDensity2[2];
    }
    
               
    // create outputs
    onputs[0] = mxCreateDoubleMatrix(nDensities1,nDensities2, mxREAL);
    double* out = mxGetPr(onputs[0]);
    
    // auxiliary variables with dummy initialization (this is necessary, see Eigen docs)
    Eigen::Map<Eigen::MatrixXd,Eigen::ColMajor> densMat1(densities1,m,nSamples1), densMat2(densities1,m,nSamples1);
    Eigen::Map<Eigen::VectorXd> densMean1(means1,m), densMean2(means1,m);

        
    // compute the distances
    if(nInputs == 5){
        // if two densities are provided
        for(mwSize i=0;i<nDensities1;i++){
            for(mwSize j=0;j<nDensities2;j++){                
                // map densities to matrices and means to vectors
                new (&densMat1) Eigen::Map<Eigen::MatrixXd>(densities1+i*m*nSamples1,m,nSamples1);
                new (&densMat2) Eigen::Map<Eigen::MatrixXd>(densities2+j*m*nSamples1,m,nSamples1);
                new (&densMean1) Eigen::Map<Eigen::VectorXd>(means1+i*m,m);
                new (&densMean2) Eigen::Map<Eigen::VectorXd>(means2+j*m,m);
                
                // compute distance
                out[i+j*nDensities1] = DiracDiracDistanceC(densMat1, densMean1, densMat2, densMean2, bmax);
            }
        }
    } else{
        // if one density is provided
        for(mwSize j=0;j<nDensities2;j++){
            for(mwSize i=j+1;i<nDensities1;i++){
                // map density to matrix and mean to vector
                new (&densMat1) Eigen::Map<Eigen::MatrixXd>(densities1+i*m*nSamples1,m,nSamples1);
                new (&densMat2) Eigen::Map<Eigen::MatrixXd>(densities1+j*m*nSamples1,m,nSamples1);
                new (&densMean1) Eigen::Map<Eigen::VectorXd>(means1+i*m,m);
                new (&densMean2) Eigen::Map<Eigen::VectorXd>(means1+j*m,m);
                    
                // compute distance                    
                out[j+i*nDensities1] = DiracDiracDistanceC(densMat1, densMean1, densMat2, densMean2, bmax);
                out[i+j*nDensities1] = out[j+i*nDensities1];
            }
        }
    }
    
    return;
}

/*------------------------------------------------------------------------------------------------------*/

double DiracDiracDistanceC(Eigen::Map<Eigen::MatrixXd,Eigen::ColMajor> &densMat1, Eigen::Map<Eigen::VectorXd> &densMean1, Eigen::Map<Eigen::MatrixXd,Eigen::ColMajor> &densMat2, Eigen::Map<Eigen::VectorXd> &densMean2, double *bmax){
    double Df = 0, Dfg = 0, Dg = 0;

    unsigned int numSamples1 = densMat1.cols();
    unsigned int m = densMat1.rows();
    
    Eigen::MatrixXd tmpMat(m,numSamples1);
    Eigen::VectorXd tmpVecDf,tmpVecDfg,tmpVecDg;//(m);
    
    for(unsigned int i=0; i<numSamples1; i++){
        for(unsigned int j=0; j<numSamples1; j++){
            
            tmpVecDf = densMat1.col(i)-densMat1.col(j);
            tmpVecDfg = densMat1.col(i)-densMat2.col(j);
            tmpVecDg = densMat2.col(i)-densMat2.col(j);
            
            Df = Df + xlog(tmpVecDf.squaredNorm());
            Dfg = Dfg + xlog(tmpVecDfg.squaredNorm());
            Dg = Dg + xlog(tmpVecDg.squaredNorm());
        }
    }
    
    tmpVecDf = densMean1 - densMean2;
    
    return  pow(M_PI, (double)m / 2)*((Df - 2 * Dfg + Dg) / (numSamples1*numSamples1) + 2 * (log(4*(*bmax)*(*bmax))-Gamma)*tmpVecDf.squaredNorm())/8;    
}