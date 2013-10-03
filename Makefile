# Makefile for DDSCAT.7.3
vers=7.3.0

# upper-level targets:
#       calltarget
#	ddscat                
#	ddpostprocess
#       vtrconvert

#--------- do NOT alter the following definitions: -------------------------
MPI_f = mpi_subs.f90 \
	mpi_bcast_char.f90 mpi_bcast_cplx.f90 mpi_bcast_int.f90\
	mpi_bcast_int2.f90 mpi_bcast_real.f90
MPI_o = mpi_subs.o \
	mpi_bcast_char.o mpi_bcast_cplx.o mpi_bcast_int.o\
	mpi_bcast_int2.o mpi_bcast_real.o
MKL_f = cxfft3_mkl.f90 mkl_dfti.f90
MKL_o = cxfft3_mkl.o mkl_dfti.o
MKL_m = mkl_dfti.mod
#---------------------------------------------------------------------------

# Here we explain the different strings that need to be defined:

#                          PRECISION

# ddscat can be compiled to use either single (sp) or double precision (dp)
# For most applications, sp appears to work fine, and is recommended for
# normal use.
# However, if ddscat is converging only very slowly or not at all, you can
# experiment with switching from sp to dp to see whether roundoff errors
# are compromising performance of the conjugate gradient solver.

# ** Important** : if you change from sp to dp, or from dp to sp:
# you must type 
#   make clean
# before typing
#   make ddscat

# to use single precision, set
# PRECISION = sp

# to use double precision, set
# PRECISION = dp

#---------------------------------------------------------------------------

#                        OpenMP support

#    OpenMP (www.openmp.org) can be used to use multiple threads on 
#    shared-memory nodes (e.g., a dual quad-core cpus can support up to 
#    8 simultaneous threads).
# If OpenMP is not installed, leave OPENMP undefined:
# DOMP	    =
# OPENMP    =

#    If OpenMP is installed, and you would like to use it:
# DOMP	    = -Dopenmp
# OPENMP    = -openmp
# or, on some systems (e.g., gfortran)
# OPENMP    = -fopenmp

#-----------------------------------------------------------------------

#                        FFTMKL support

#    If the Intel MKL library is installed, then set
# CXFFTMKL.f = $(MKL_f)
# CXFFTMKL.o = $(MKL_o)
# MKLM       = $(MKL_m)

#    If the Intel MKL library is not installed on system, then set
# CXFFTMKL.f = cxfft3_mkl_fake.f90
# CXFFTMKL.o = cxfft3_mkl_fake.o
# MKLM       =

#-----------------------------------------------------------------------------
#                          MPI support 

#    Module DDSCAT.f90 require very minor editing to prepare it for
#    either non-MPI use or use on a system with MPI support.
#    Consult either the UserGuide (section 24.1) or comments within DDSCAT.f90
#    for instructions on which line needs to be enabled
#    and which line needs to be commented out.

#    The strings MPI.f and MPI.o need to be appropriately defined.

#    If MPI support is desired and the MPI library is installed, these are
# MPI.f	    = $(MPI_f)
# MPI.o	    = $(MPI_o)
# DMPI	    = -Dmpi
#
#    If MPI support is not required, these are
# MPI.f	    = mpi_fake.f90
# MPI.o	    = mpi_fake.o
# DMPI       =

#******************************************************************************
#
#                   Compiler and options
#
#    FC specifies the fortran 90 compiler
#    FFLAGS are compilation options
#    LFLAGS are flags for linking
#
#############################################################################
#                      Examples

# 1.  gfortran compiler
#     dp + no MKL + OpenMP + no MPI

# define the following:
PRECISION	= dp
CXFFTMKL.f	= cxfft3_mkl_fake.f90
CXFFTMKL.o	= cxfft3_mkl_fake.o
MKLM		=
DOMP		= -Dopenmp
OPENMP		= -fopenmp
MPI.f	    = mpi_fake.f90
MPI.o	    = mpi_fake.o
DMPI       =
FC		= gfortran
FFLAGS		= -O2 -march=native -ffloat-store
LFLAGS	 	=

# 1.b  gfortran compiler
#     dp + no MKL + no OpenMP + no MPI

# define the following:
#PRECISION	= sp
#CXFFTMKL.f	= cxfft3_mkl_fake.f90
#CXFFTMKL.o	= cxfft3_mkl_fake.o
#MKLM		=
#DOMP		=
#OPENMP		=
#MPI.f		= mpi_fake.f90
#MPI.o		= mpi_fake.o
#DMPI		=
#FC		= gfortran
#FFLAGS		= -O2 -march=native
#LFLAGS	 	=

#-------------------------------------------------------------------

# 2.  g95 compiler
#     sp + no MKL + no OpenMP + no MPI

# define the following:
#PRECISION	= sp
#CXFFTMKL.f	= cxfft3_mkl_fake.f90
#CXFFTMKL.o	= cxfft3_mkl_fake.o
#MKLM		=
#DOMP		=
#OPENMP		=
#MPI.f		= mpi_fake.f90
#MPI.o		= mpi_fake.o
#DMPI		=
#FC		= g95
#FFLAGS		= -O2
#LFLAGS	 	=

#----------------------------------------------------------------------
# 3.  ifort compiler
#     sp + no MKL + no OpenMP + no MPI 

# define the following:
#PRECISION	= sp
#CXFFTMKL.f	= cxfft3_mkl_fake.f90
#CXFFTMKL.o	= cxfft3_mkl_fake.o
#MKLM		=
#DOMP		=
#OPENMP		=
#MPI.f		= mpi_fake.f90
#MPI.o		= mpi_fake.o
#DMPI		=
#FC		= ifort
#FFLAGS		= -O2
#LFLAGS		=

#----------------------------------------------------------------------
# 3.b  ifort compiler
#     dp + no MKL + no OpenMP + no MPI,

# define the following:
#PRECISION	= dp
#CXFFTMKL.f	= cxfft3_mkl_fake.f90
#CXFFTMKL.o	= cxfft3_mkl_fake.o
#MKLM		=
#DOMP		=
#OPENMP		=
#MPI.f 		= mpi_fake.f90
#MPI.o		= mpi_fake.o
#DMPI		=
#FC		= ifort
#FFLAGS		= -O2 -C
#LFLAGS		=

#----------------------------------------------------------------------
# 5.  ifort compiler
#     sp + MKL + no OpenMP + no MPI,

# define the following:
#PRECISION	= sp
#CXFFTMKL.f	= $(MKL_f)
#CXFFTMKL.o	= $(MKL_o)
#MKLM		= $(MKL_m)
#DOMP		=
#OPENMP		=
#MPI.f		= mpi_fake.f90
#MPI.o		= mpi_fake.o
#DMPI		=
#FC		= ifort
#FFLAGS		= -O2
#LFLAGS	 	=

#----------------------------------------------------------------------
# 6.  ifort compiler
#     sp + no MKL + OpenMP + no MPI,

# define the following:
#PRECISION	= sp
#CXFFTMKL.f	= cxfft3_mkl_fake.f90
#CXFFTMKL.o	= cxfft3_mkl_fake.o
#MKLM		=
#DOMP		= -Dopenmp
#OPENMP		= -openmp
#MPI.f		= mpi_fake.f90
#MPI.o		= mpi_fake.o
#DMPI		=
#FC		= ifort
#FFLAGS		= -O2
#LFLAGS		=

#----------------------------------------------------------------------
# 7.  ifort compiler
#     sp + MKL + OpenMP + no MPI

#     following definitions work for artemis
# before compiling, type
#   module purge
#   module load intel-mkl

# define the following:
#PRECISION	= sp
#CXFFTMKL.f	= $(MKL_f)
#CXFFTMKL.o	= $(MKL_o)
#MKLM		= $(MKL_m)
#DOMP		= -Dopenmp
#OPENMP		= -openmp
#MPI.f		= mpi_fake.f90
#MPI.o		= mpi_fake.o
#DMPI		=
#FC		= ifort
#FFLAGS		= -O2
#LFLAGS		= -traceback -lmkl_em64t -lmkl_intel_thread -lmkl_core \
#		-lguide -lpthread -lmkl_intel_lp64

#----------------------------------------------------------------------
# 8.  ifort compiler (via mpif90)
#     sp + MKL + OpenMP + MPI

#     following definitions work for artemis
# before compiling, type
#   module purge
#   module load intel-mkl openmpi

# define the following:
#PRECISION	= sp
#CXFFTMKL.f	= $(MKL_f)
#CXFFTMKL.o	= $(MKL_o)
#MKLM		= $(MKL_m)
#DOMP		= -Dopenmp
#OPENMP		= -openmp
#MPI.f		= $(MPI_f)
#MPI.o		= $(MPI_o)
#DMPI		= -Dmpi
#FC		= mpif90
#FFLAGS		= -O2
#LFLAGS		= -traceback -lmkl_em64t -lmkl_intel_thread -lmkl_core \
#		-lguide -lpthread -lmkl_intel_lp64 -lmpi

#******************************************************************************
#
#                  End of option specifications.
#
#******************************************************************************

# general rule for compilation of most .o files:

%.o: %.f90 ddprecision.mod ddcommon_1.mod cgmodule.mod
	$(FC) -c $(FFLAGS) $(OPENMP) $< -o $@

# special cases:

DDSCAT.o: DDSCAT.f90 ddprecision.mod ddcommon_1.mod cgmodule.mod
	cpp -P -traditional-cpp $(DMPI) $(DOMP) DDSCAT.f90 DDSCAT_cpp.f90
	$(FC) -c $(FFLAGS) $(OPENMP) DDSCAT_cpp.f90 -o DDSCAT.o
	rm DDSCAT_cpp.f90

DDVTR.o: DDVTR.f90 ddprecision.mod vtr.mod
	$(FC) -c $(FFLAGS) DDVTR.f90 -o DDVTR.o

DDPOSTPROCESS.o: DDPOSTPROCESS.f90 ddprecision.mod readnf_bcom.mod \
	readnf_ecom.mod vtr.mod
	$(FC) -c $(FFLAGS) $(OPENMP) DDPOSTPROCESS.f90 \
	-o DDPOSTPROCESS.o

bself.o: bself.f90 ddprecision.mod
	cpp -P -traditional-cpp $(DOMP) bself.f90 bself_cpp.f90
	$(FC) -c $(FFLAGS) $(OPENMP) bself_cpp.f90 -o bself.o
	rm bself_cpp.f90

cgcommon.o: cgcommon.f90 ddprecision.mod
	cpp -P -traditional-cpp -D$(PRECISION) cgcommon.f90 cgcommon_cpp.f90
	$(FC) -c $(FFLAGS) $(OPENMP) cgcommon_cpp.f90 -o cgcommon.o
	rm cgcommon_cpp.f90

cxfft3_mkl.o: cxfft3_mkl.f90 ddprecision.mod mkl_dfti.mod
	$(FC) -c $(FFLAGS) $(OPENMP) cxfft3_mkl.f90 \
	-o cxfft3_mkl.o

eself.o: eself.f90 ddprecision.mod
	cpp -P -traditional-cpp $(DOMP) eself.f90 eself_cpp.f90
	$(FC) -c $(FFLAGS) $(OPENMP) eself_cpp.f90 -o eself.o
	rm eself_cpp.f90

scat.o: scat.f90 ddprecision.mod ddcommon_1.mod
	cpp -P -traditional-cpp $(DOMP) scat.f90 scat_cpp.f90
	$(FC) -c $(FFLAGS) $(OPENMP) scat_cpp.f90 -o scat.o
	rm scat_cpp.f90

readnf.o: readnf.f90 ddprecision.mod readnf_bcom.mod readnf_ecom.mod
	$(FC) -c $(FFLAGS) $(OPENMP) readnf.f90 \
	-o readnf.o

# dependencies for ddscat:
OBJS =  DDSCAT.o\
	alphadiag.o\
	blas.o\
	bself.o \
	ccgpack.o\
	cgcommon.o\
	cglib2.o\
	cglib3.o\
	cgsarkar2.o\
	cgsarkar3.o\
	cisi.o\
	copyit.o\
	cprod.o\
	cxfft3n.o\
	$(CXFFTMKL.o)\
	cxfftw_fake.o\
	ddcommon.o\
	dielec.o\
	divide.o\
	dsyevj3.o\
	dummy.o\
	errmsg.o\
	eself.o\
	evala.o\
	evale.o\
	besseli0.o \
	besseli1.o \
	besselk0.o \
	besselk1.o \
	evalq.o\
	extend.o\
	gasdev.o\
	getfml.o\
	getmueller.o\
        gpbicg.o\
	gpfa.o\
	interp.o\
	$(MPI.o)\
	namer.o\
	namer2.o\
        namid.o\
        nearfield.o\
	nuller.o\
	orient.o\
	p_lm.o\
	pbcscavec.o\
	pim.o\
	prinaxis.o\
        qmrpim2.o\
	ran3.o\
	reapar.o\
	reashp.o\
	reduce.o\
	restore.o\
	rot2.o\
	rotate.o\
	scat.o\
	scavec.o\
        tangcg.o\
	tar2el.o\
	tar2sp.o\
	tar3el.o\
	taranirec.o\
	tarblocks.o\
	tarcel.o\
	tarcyl.o\
	tarcylcap.o\
	tarell.o\
	target.o\
	targspher.o\
	tarhex.o\
	tarlyrslab.o\
	tarnas.o\
	tarnsp.o\
	tarpbxn.o\
	tarprsm.o\
	tarrctblk3.o\
	tarrctell.o\
	tarrec.o\
	tarrecrec.o\
	tarslbhol.o\
	tartet.o\
	timeit.o\
        unreduce.o\
	version.o\
	wrimsg.o\
	writebin.o\
	writefml.o\
	writepol.o\
	writesca.o\
	zbcg2wp.o

# dependencies for calltarget:

OBJS2 = CALLTARGET.o\
	ddcommon.o\
	dsyevj3.o\
	errmsg.o\
	gasdev.o\
	p_lm.o\
	prinaxis.o\
	ran3.o\
	reashp.o\
	sizer.o\
	tar2el.o\
	tar2sp.o\
	tar3el.o\
	taranirec.o\
	tarblocks.o\
	tarcel.o\
	tarcyl.o\
	tarcylcap.o\
	tarell.o\
	target.o\
	targspher.o\
	tarhex.o\
	tarlyrslab.o\
	tarnas.o\
	tarnsp.o\
	tarpbxn.o\
	tarprsm.o\
	tarrctblk3.o\
	tarrctell.o\
	tarrec.o\
	tarrecrec.o\
	tarslbhol.o\
	tartet.o\
	wrimsg.o

# dependencies for ddpostprocess:

OBJS3 =	DDPOSTPROCESS.o\
	readnf_bcom.o\
	readnf_ecom.o\
	readnf.o\
	vtr.o

# dependencies for vtrconvert:

OBJS4 = VTRCONVERT.o\
	vtr.o

all:	ddscat calltarget ddpostprocess vtrconvert

ddscat:	ddprecision.mod ddcommon_1.mod $(MKLM) $(OBJS)
	@echo 'LOADEDMODULES='$(LOADEDMODULES)
	@echo 'LOADEDMODULES_modshare='$(LOADEDMODULES_modshare)
	@echo 'LD_LIBRARY_PATH='$(LD_LIBRARY_PATH)
	@echo 'LD_LIBRARY_PATH_modshare='$(LD_LIBRARY_PATH_modshare)
	@echo 'MKLROOT='$(MKLROOT)
	@echo 'CPATH='$(CPATH)
	@echo 'FPATH='$(FPATH)
	@echo 'DYLD_LIBRARY_PATH='$(DYLD_LIBRARY_PATH)
	@echo 'INCLUDE='$(INCLUDE)
	@echo 'LIBRARY_PATH='$(LIBRARY_PATH)
	$(FC) -o ddscat \
	$(OBJS) $(LFLAGS) $(OPENMP)

calltarget: ddprecision.mod ddcommon_1.mod $(OBJS2)
	$(FC) -o calltarget \
	$(OBJS2) $(LFLAGS)

ddpostprocess: ddprecision.mod readnf_bcom.mod readnf_ecom.mod vtr.mod $(OBJS3)
	$(FC) -o ddpostprocess \
	$(OBJS3) $(LFLAGS) 

vtrconvert: ddprecision.mod vtr.o $(OBJS4)
	$(FC) -o vtrconvert \
	$(OBJS4) $(LFLAGS)

#--------------- modules ---------------------------------------------

cgmodule.mod: ddprecision.mod cgmodule.f90
	$(FC) -c cgmodule.f90 -o cgmodule.o

ddprecision.mod: ddprecision.f90
	cpp -P -traditional-cpp -D$(PRECISION) ddprecision.f90 \
	ddprecision_cpp.f90
	$(FC) -c $(FFLAGS) ddprecision_cpp.f90 -o ddprecision.o
	rm ddprecision_cpp.f90

ddcommon_1.mod:	ddprecision.mod ddcommon.f90
	$(FC) -c ddcommon.f90 -o ddcommon.o

mkl_dfti.mod: mkl_dfti.f90
	$(FC) -c mkl_dfti.f90 -o mkl_dfti.o

readnf_bcom.mod: readnf_bcom.f90
	$(FC) -c readnf_bcom.f90 -o readnf_bcom.o

readnf_ecom.mod: readnf_ecom.f90
	$(FC) -c readnf_ecom.f90 -o readnf_ecom.o

vtr.mod: vtr.f90
	$(FC) -c vtr.f90 -o vtr.o

#---------------------------------------------------------------------

clean:;	rm -f *.o make.out*  *.mod

veryclean: clean
	rm -f calltarget ddscat ddpostprocess
