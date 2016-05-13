# GP
Reference implementations of (my) algorithms for Gaussian Process regression.

The following external code is used:
- Eigen from https://eigen.tuxfamily.org
- LAPJV: Linear Assignment Problem algorithm by R. Jonker and A. Volgenant
	(see "A Shortest Augmenting Path Algorithm for Dense and Sparse Linear Assignment Problems," Computing 38, 
	325-340, 1987)
	
Legacy implementation using BLAS are in distances_with_BLAS. The implementation using Eigen turned out to be faster. Thus, the BLAS-based implementations will no longer be maintained.