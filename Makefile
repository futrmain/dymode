
MKL_ROOT = /pdc/vol/i-compilers/12.1.5/composer_xe_2011_sp1.11.339/mkl
MKL_LINK =  -Wl,--start-group  /pdc/vol/i-compilers/12.1.5/composer_xe_2011_sp1.11.339/mkl/lib/intel64/libmkl_intel_lp64.a /pdc/vol/i-compilers/12.1.5/composer_xe_2011_sp1.11.339/mkl/lib/intel64/libmkl_sequential.a /pdc/vol/i-compilers/12.1.5/composer_xe_2011_sp1.11.339/mkl/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm

CFS = /cfs/klemming/nobackup/f/fromain/peigen/
Yincl = yaml/include
Eincl = ../eigen-eigen-3afaf217ffc2

MKLROOT = $(MKL_ROOT)
PMKL_LINK = $(MKLROOT)/lib/intel64/libmkl_scalapack_lp64.a -Wl,--start-group  $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_intel_thread.a $(MKLROOT)/lib/intel64/libmkl_core.a $(MKLROOT)/lib/intel64/libmkl_blacs_intelmpi_lp64.a -Wl,--end-group -liomp5 -lpthread -lm
	
main:	
	CC main.cpp -O0 -o $(CFS)$@ -I$(Eincl) -I$(MKL_ROOT)/include $(MKL_LINK)

