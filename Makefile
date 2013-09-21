#!/bin/sh

# Makefile for DDSCAT.7.1.0 
vers=7.1.0

# upper-level targets:
#	ddscat                
#	ddfield
#	ddpol            

#--------- do NOT alter the following definitions: -------------------------
MPI_f   = mpi_subs.f90 \
	mpi_bcast_char.f90 mpi_bcast_cplx.f90 mpi_bcast_int.f90 \
	mpi_bcast_int2.f90 mpi_bcast_real.f90
MPI_o	= mpi_subs.o \
	mpi_bcast_char.o mpi_bcast_cplx.o mpi_bcast_int.o \
	mpi_bcast_int2.o mpi_bcast_real.o
MKL_f	= cxfft3_mkl.f90 mkl_dfti.f90
MKL_o	= cxfft3_mkl.o mkl_dfti.o
MKL_m	= mkl_dfti.mod
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
#     sp + no MKL + no OpenMP + no MPI

# define the following:
PRECISION	= dp
CXFFTMKL.f	= cxfft3_mkl_fake.f90
CXFFTMKL.o	= cxfft3_mkl_fake.o
MKLM		=
DOMP		=
OPENMP		=
MPI.f		= mpi_fake.f90
MPI.o		= mpi_fake.o
DMPI		=
FC		= gfortran
FFLAGS		= -O2
LFLAGS	 	=

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
# 3.  NAG f95 compiler
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
#FC		= f95
#FFLAGS		= -O2
#LFLAGS		=

#----------------------------------------------------------------------
# 4.  ifort compiler
#     sp + no MKL + no OpenMP + no MPI,

# define the following:
#PRECISION	= sp
#CXFFTMKL.f	= cxfft3_mkl_fake.f90
#CXFFTMKL.o	= cxfft3_mkl_fake.o
#MKLM		=
#DOMP		=
#OPENMP		=
#MPI.f 		= mpi_fake.f90
#MPI.o		= mpi_fake.o
#DMPI		=
#FC		= ifort
#FFLAGS		= -O2
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
#MPI.f		= $(MPI_f)
#MPI.o		= $(MPI_o)
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

%.o: %.f90 ddprecision.mod ddcommon_1.mod
	$(FC) -c $(FFLAGS) $(OPENMP) $< -o $@

# special cases:

ddscat.o: ddscat.f90 ddprecision.mod ddcommon_1.mod
	cpp -P -traditional-cpp $(DMPI) $(DOMP) ddscat.f90 \
	DDSCAT_cpp.f90
	$(FC) -c $(FFLAGS) $(OPENMP) DDSCAT_cpp.f90 -o ddscat.o
	rm DDSCAT_cpp.f90

ddfield.o: ddfield.f90 ddprecision.mod
	cpp -P -traditional-cpp $(DOMP) ddfield.f90 DDfield_cpp.f90
	$(FC) -c $(FFLAGS) $(OPENMP) DDfield_cpp.f90 -o ddfield.o
	rm DDfield_cpp.f90

cgcommon.o: cgcommon.f90 ddprecision.mod
	cpp -P -traditional-cpp -D$(PRECISION) cgcommon.f90 cgcommon_cpp.f90
	$(FC) -c $(FFLAGS) $(OPENMP) cgcommon_cpp.f90 -o cgcommon.o
	rm cgcommon_cpp.f90

eself.o: eself.f90 ddprecision.mod
	cpp -P -traditional-cpp $(DOMP) eself.f90 eself_cpp.f90
	$(FC) -c $(FFLAGS) $(OPENMP) eself_cpp.f90 -o eself.o
	rm eself_cpp.f90

scat.o: scat.f90 ddprecision.mod ddcommon_1.mod
	cpp -P -traditional-cpp $(DOMP) scat.f90 scat_cpp.f90
	$(FC) -c $(FFLAGS) $(OPENMP) scat_cpp.f90 -o scat.o
	rm scat_cpp.f90

cxfft3_mkl.o: cxfft3_mkl.f90 ddprecision.mod mkl_dfti.mod
	$(FC) -c $(FFLAGS) $(OPENMP) cxfft3_mkl.f90 \
	-o cxfft3_mkl.o

OBJS	= ddscat.o \
	alphadiag.o \
	blas.o \
	ccgpack.o \
	cgcommon.o \
	copyit.o \
	cprod.o \
	cxfft3n.o \
	$(CXFFTMKL.o) \
	cxfftw_fake.o \
	ddcommon.o \
	dielec.o \
	divide.o \
	dsyevj3.o \
	dummy.o \
	errmsg.o \
	eself.o \
	evala.o \
        evale.o \
	besseli0.o \
	besseli1.o \
	besselk0.o \
	besselk1.o \
	evalq.o \
	extend.o \
	gasdev.o \
	getfml.o \
	getmueller.o \
	gpfa.o \
	interp.o \
	$(MPI.o) \
	namer.o \
	namer2.o \
	namid.o \
	nuller.o \
	orient.o \
	p_lm.o \
	pbcscavec.o \
	pim.o \
	prinaxis.o \
	ran3.o \
	reapar.o \
	reashp.o \
	reduce.o \
	restore.o \
	rot2.o \
	rotate.o \
	scat.o \
	scavec.o \
	tar2el.o \
	tar2sp.o \
	tar3el.o \
	taranirec.o \
	tarblocks.o \
	tarcel.o \
	tarcyl.o \
	tarcylcap.o \
	tarell.o \
	target.o \
	targspher.o \
	tarhex.o \
	tarlyrslab.o \
	tarnas.o \
	tarnsp.o \
	tarpbxn.o \
	tarprsm.o \
	tarrctblk3.o \
	tarrec.o \
	tarrecrec.o \
	tarslbhol.o \
	tartet.o \
	timeit.o \
	version.o \
	wrimsg.o \
	writebin.o \
	writefml.o \
	writepol.o \
	writesca.o \
	zbcg2wp.o

# dependencies for DDfield:

OBJS2	= ddfield.o \
	readpol.o  \
        besseli0.o \
        besseli1.o \
        besselk0.o \
        besselk1.o

# dependencies for DDpol:

OBJS3	= ddpol.o \
	readpol.o

# dependencies for calltarget:

OBJS4	= calltarget.o\
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
	tarrec.o\
	tarrecrec.o\
	tarslbhol.o\
	tartet.o\
	wrimsg.o

all:	ddscat ddfield ddpol calltarget

ddscat:	ddprecision.mod ddcommon_1.mod $(MKLM) $(OBJS)
	@echo 'LOADEDMODULES='$(LOADEDMODULES)
	@echo 'LOADEDMODULES_modshare='$(LOADEDMODULES_modshare)
	@echo 'LD_LIBRARY_PATH='$(LD_LIBRARY_PATH)
	@echo 'LD_LIBRARY_PATH_modshare='$(LD_LIBRARY_PATH_modshare)
	$(FC) $(OBJS) $(LFLAGS) $(OPENMP) -o ddscat

ddfield: ddprecision.mod $(OBJS2)
	$(FC) $(OBJS2) $(OPENMP) $(LFLAGS) -o ddfield

ddpol:	$(OBJS3)
	$(FC) $(OBJS3) $(LFLAGS) -o ddpol

calltarget: ddprecision.mod ddcommon_1.mod $(OBJS4)
	$(FC) $(OBJS4) $(LFLAGS) -o calltarget

#--------------- modules ---------------------------------------------

ddprecision.mod: ddprecision.f90
	cpp -P -traditional-cpp -D$(PRECISION) ddprecision.f90 \
	ddprecision_cpp.f90
	$(FC) -c $(FFLAGS) ddprecision_cpp.f90 -o ddprecision.o
	rm ddprecision_cpp.f90

ddcommon_1.mod:	ddprecision.mod ddcommon.f90
	$(FC) -c ddcommon.f90

mkl_dfti.mod: mkl_dfti.f90
	$(FC) -c mkl_dfti.f90

#---------------------------------------------------------------------

clean:;	rm -f *.o make.out*  *.mod

veryclean: clean
	rm -f ddscat ddfield ddpol *~
