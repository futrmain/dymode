SRC = aliasmanager.o \
     emitterutils.o \
	 ostream.o \
     simplekey.o \
 	 binary.o \
     exp.o \
	 parser.o \
     singledocparser.o \
	 conversion.o \
	 iterator.o \
	 regex.o \
	 stream.o \
	 directives.o   \
	 nodebuilder.o \
	 scanner.o  \
	 tag.o \
	 emitfromevents.o \
	 node.o  \
	 scanscalar.o \
	 emitter.o \
	 nodeownership.o \
	 scantag.o \
	 emitterstate.o \
	 null.o  \
	 scantoken.o

# This is for Lindgren	 
#MKL_ROOT = /pdc/vol/i-compilers/12.1.5/composer_xe_2011_sp1.11.339/mkl
#MKL_LINK =  -Wl,--start-group  /pdc/vol/i-compilers/12.1.5/composer_xe_2011_sp1.11.339/mkl/lib/intel64/libmkl_intel_lp64.a /pdc/vol/i-compilers/12.1.5/composer_xe_2011_sp1.11.339/mkl/lib/intel64/libmkl_sequential.a /pdc/vol/i-compilers/12.1.5/composer_xe_2011_sp1.11.339/mkl/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm

CFS = /pfs/nobackup/home/r/romainf/WengProject/oscill350/
Yincl = ../yaml-cpp/include
Ylib = ../yaml-cpp/lib
Eincl = ../eigen-eigen-3afaf217ffc2
PEincl = ../peigen
Bincl = /pfs/nobackup/home/r/romainf/boost/boost_1_53_0

H5incl = /pfs/nobackup/home/r/romainf/hdf5-1.8.7-parallel/gcc/include
H5lib = /pfs/nobackup/home/r/romainf/hdf5-1.8.7-parallel/gcc/lib/libhdf5.a
H5CC = /pfs/nobackup/home/r/romainf/hdf5-1.8.10-patch1/hdf5/bin/h5pcc

# This is for Abisko
MKLROOT = /lap/intel-mkl/10.3.3.174/mkl
#PMKL_LINK = $(MKLROOT)/lib/intel64/libmkl_scalapack_lp64.a -Wl,--start-group  $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_intel_thread.a $(MKLROOT)/lib/intel64/libmkl_core.a $(MKLROOT)/lib/intel64/libmkl_blacs_intelmpi_lp64.a -Wl,--end-group -liomp5 -lpthread -lm

MKL_LINK_IMPI = $(MKLROOT)/lib/intel64/libmkl_scalapack_lp64.a -Wl,--start-group  $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_sequential.a $(MKLROOT)/lib/intel64/libmkl_core.a $(MKLROOT)/lib/intel64/libmkl_blacs_intelmpi_lp64.a -Wl,--end-group -lpthread -lm 

MKL_LINK_GCC = $(MKLROOT)/lib/intel64/libmkl_scalapack_lp64.a -Wl,--start-group  $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_sequential.a $(MKLROOT)/lib/intel64/libmkl_core.a $(MKLROOT)/lib/intel64/libmkl_blacs_openmpi_lp64.a -Wl,--end-group -lpthread -lm

MKL_CFLAGS_GCC = -m64 -I$(MKLROOT)/include

LDFLAGS=$$MKL_LDFLAGS
CFLAGS=$$MKL_INCLUDE

	
main:	
	CC main.cpp -O0 -o $(CFS)$@ -I$(Eincl) -I$(PEincl) -I$(MKL_ROOT)/include $(MKL_LINK)
	
main2:	
	CC main2.cpp -O0 -o $(CFS)$@ -I$(Eincl) -I$(PEincl) -I$(MKL_ROOT)/include $(MKL_LINK) -wd2196
	
main2O3:	
	CC main2.cpp -O3 -o $(CFS)$@ -I$(Eincl) -I$(PEincl) -I$(MKL_ROOT)/include $(MKL_LINK) -wd2196
	
main3:	
	$(H5CC) main3.cpp -O0 -o $(CFS)$@ -I$(Bincl) -I$(Eincl) -I$(Yincl) $(addprefix $(Ylib)/gcc-4.6.3/,$(SRC)) -I$(PEincl) -I$(MKL_ROOT)/include $(MKL_LINK) -wd2196 -I$(H5incl) -L$(H5lib) -lhdf5 $(CFLAGS)
	
main5:	
	mpiCC main5.cpp -O0 -o $(CFS)/$@ -I$(Bincl) -I$(Eincl) -I$(Yincl)  $(Ylib)/libyamlcpp.a -I$(PEincl)   $(MKL_LINK_GCC) $(MKL_CFLAGS_GCC) -wd2196 -I$(H5incl) $(H5lib)
	
main5O3:	
	mpiCC main5.cpp -O3 -o $(CFS)/$@ -I$(Bincl) -I$(Eincl) -I$(Yincl)  $(Ylib)/libyamlcpp.a -I$(PEincl)   $(MKL_LINK_GCC) $(MKL_CFLAGS_GCC) -wd2196 -I$(H5incl) $(H5lib)

main3O3:	
	CC main3.cpp -O3 -o $(CFS)$@ -I$(Bincl) -I$(Eincl) -I$(Yincl) $(addprefix $(Ylib)/gcc-4.6.3/,$(SRC)) -I$(PEincl) -I$(MKL_ROOT)/include $(MKL_LINK) -wd2196
