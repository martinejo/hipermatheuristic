CC=icc
COPT=-O3
MKLPATH=/opt/intel/composerxe/mkl/lib/intel64
MKLINCLUDE=/opt/intel/composerxe/mkl/include

SOS: Bender.c
	$(CC) $(COPT) -I$(MKLINCLUDE) -I$(AUXILIAR) -L$(MKLPATH) \
	Bender.c -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lmkl_lapack95_lp64 \
        -openmp -fPIC -I/home_nfs/mgonzes/CPLEX/cplex/include -L/home_nfs/mgonzes/CPLEX/cplex/lib/x86-64_linux/static_pic -lcplex -lm -lpthread -g -liomp5 -lpthread -lm -o Bender
		
		
		
		
