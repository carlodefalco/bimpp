../configure --prefix=${HOME}/bimpp \
	     CC=mpicc CXX=mpicxx F77=mpif90 \
	     CPPFLAGS="-DBIM_TIMING -UNDEBUG -fopenmp -std=c++11" \
	     CXXFLAGS="-g -O3" \
	     --with-blas-lapack="-L${mkOpenblasLib} -lopenblas" \
	     --with-lis-home=${mkLisPrefix} \
	     --with-mumps-home=${mkMumpsPrefix} \
	     --with-mumps-extra-libs="-L${mkMumpsLib} -L${mkScotchLib} -lptscotcherr" \
	     --with-octave-home=${mkOctavePrefix}/bin \
	     --with-p4est-home=${mkP4estPrefix}
