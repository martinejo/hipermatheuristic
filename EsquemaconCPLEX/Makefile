CC=icc
COPT=-O3
MKLPATH=/opt/intel/composerxe/mkl/lib/intel64
MKLINCLUDE=/opt/intel/composerxe/mkl/include

genetico: main.c chains_tools.c generador.c rngs.c eficiencia.c hiperheuristica.c genetico.c
	$(CC) $(COPT) -I$(MKLINCLUDE) -I$(AUXILIAR) -L$(MKLPATH) \
	main.c chains_tools.c generador.c rngs.c eficiencia.c hiperheuristica.c genetico.c -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lmkl_lapack95_lp64 \
        -fPIC -I/home_nfs/mgonzes/CPLEX/cplex/include -L/home_nfs/mgonzes/CPLEX/cplex/lib/x86-64_linux/static_pic -lcplex -lm -lpthread -g -liomp5 -lpthread -lm -o genetico