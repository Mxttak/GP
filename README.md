# GP
Reference implementations of (my) algorithms for Gaussian Process regression.

I used the self-contained built OpenBLAS implementation from https://bitbucket.org/carlkl/mingw-w64-for-python/downloads on Win10. To build the sources with this library, run, e.g.,

	mex DiracDiracDistance.cpp DiracDiracDistanceC.cpp libopenblaspy.dll.a