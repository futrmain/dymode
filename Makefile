###################################################
###                                             ###
###     Step 1. Set the compiler                ###
###                                             ###
###################################################

CC = /home/fromain/libs/openmpi-gcc/bin/mpiCC
CCOPTS =  -std=c++11 -Wfatal-errors -fpermissive
CCOPTS_TRAIL = -lgfortran


###################################################
###                                             ###
###     Step 2. Set external dependencies paths ###
###                                             ###
###################################################
MPI = 
HDF5 = /home/fromain/libs/gcc/hdf5
BOOST = /home/fromain/libs/gcc/boost


###################################################
###                                             ###
###     Step 3. Set up the ScaLAPACK            ###
###                                             ###
###################################################

###################################################
###	Case 1: Generic ScaLAPACK		###
###################################################
BLAS = /home/fromain/libs/gcc/BLAS
LAPACK = /home/fromain/libs/gcc/lapack-3.5.0
SCALAPACK = /home/fromain/libs/gcc/scalapack

PACKLIBS = -lscalapack -llapack -lblas

###################################################
###     Case 2: MKL (use MKL Link Line Advisor) ###
###################################################
MKLROOT = 
USE_THIS_LINK_LINE = 
COMPILER_OPTIONS = 				  # This is typically "-m64 -I$(MKLROOT)/include"




###################################################
###                                             ###
###     Text below should not be modified       ###
###                                             ###
###################################################

###################################################
###     Shipped libraries			###
###################################################
EIGEN = ./eigen
PEIGEN = ./peigen
TICTOC = ./tic-toc-profiler
YAMLHPP = ./yaml-hpp
TCLAP = ./tclap-1.2.1


###################################################
###     Link line	                        ###
###################################################
LINK = $(HDF5)/lib/libhdf5.a -lz


###################################################
###     Choose between MKL and regular ScaLAPACK###
###################################################
ifeq ($(MKLROOT),)
	PACKDIRS = -L$(BLAS) -L$(LAPACK) -L$(SCALAPACK)
	PACKLINK = $(PACKLIBS)
else
	PACKDIRS = $(COMPILER_OPTIONS)
	PACKLINK = $(USE_THIS_LINK_LINE)
endif


###################################################
###     Gather everything in single macros 	###
###################################################
INCLUDES = -I$(MPI)/include -I$(EIGEN) -I$(PEIGEN) -I$(BOOST) -I$(TICTOC) -I$(YAMLHPP) -I$(HDF5)/include -I$(TCLAP)/include 
LIBDIRS = -L$(MPI)/lib -L$(HDF5)/lib $(PACKDIRS)
LIBS = $(LINK) $(PACKLINK)


###################################################
###                                             ###
###     Build targets			        ###
###                                             ###
###################################################

debug:
	$(CC) $(CCOPTS) dymode/dymode.cpp $(INCLUDES) $(LIBDIRS) $(LIBS) $(CCOPTS_TRAIL) -o dymode.out

clean:
	rm dymode.out