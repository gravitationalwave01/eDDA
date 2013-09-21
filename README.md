eDDA
====

Electron-Driven Discrete-Dipole Approximation (e-DDA): Extension of the optical DDA to include an electron beam source to simulate experiments in electron energy-loss spectroscopy (EELS).

PLEASE NOTE. This code is an extension of Draine and Flatau's DDSCAT 7.1. 
The following files have been unmodified from DDSCAT 7.1:

alphadiag.f90 
blas.f90 
calltarget.f90 
ccgpack.f90 
cgcommon.f90 
copyit.f90 
cprod.f90 
cxfft3_mkl.f90 
cxfft3_mkl_fake.f90 
cxfft3n.f90 
xfftw_fake.f90 
ddcommon.f90 
ddcommon_1.mod 
ddcommon_10.mod 
eddcommon_2.mod 
ddcommon_3.mod 
ddcommon_4.mod 
ddcommon_5.mod 
ddcommon_6.mod 
ddcommon_7.mod 
ddcommon_8.mod 
ddcommon_9.mod 
ddfield.f90 
ddpol.f90 
ddprecision.f90 
ddprecision.mod 
dielec.f90 
divide.f90 
dsyevj3.f90 
dummy.f90 
errmsg.f90 
eself.f90 
evala.f90 
evale.xbeam.f90 
extend.f90 
gasdev.f90 
getmueller.f90 
gpfa.f90 
interp.f90 
mpi_bcast_char.f90 
mpi_bcast_cplx.f90 
mpi_bcast_int.f90 
mpi_bcast_int2.f90 
mpi_bcast_real.f90 
mpi_fake.f90 
mpi_subs.f90 
namer.f90 
namer2.f90 
namid.f90 
nuller.f90 
orient.f90 
p_lm.f90 
pbcscavec.f90 
pim.f90 
prinaxis.f90 
ran3.f90 
readpol.f90 
reashp.f90 
reduce.f90 
restore.f90 
rot2.f90 
rotate.f90 
scat.f90 
scavec.f90 
sizer.f90 
tar2el.f90 
tar2sp.f90 
tar3el.f90 
taranirec.f90 
tarblocks.f90 
tarcel.f90 
tarcyl.f90 
tarcylcap.f90 
tarell.f90 
target.f90 
targspher.f90 
tarhex.f90 
tarlyrslab.f90 
tarnas.f90 
tarnsp.f90 
tarpbxn.f90 
tarprsm.f90 
tarrctblk3.f90 
tarrec.f90 
tarrecrec.f90 
tarslbhol.f90 
tartet.f90 
timeit.f90 
version.f90 
wrimsg.f90 
writebin.f90 
writefml.f90 
writepol.f90 
zbcg2wp.f90


The following files have either been modified or were written entirely by the Masiello Research Group 
at the University of Washington: http://faculty.washington.edu/masiello/Masiello_Group_Website/Home.html

besseli0.f90 
besseli1.f90 
besselk0.f90 
besselk1.f90 
ddscat.f90 
evale.f90 
evalq.f90 
getfml.f90 
Makfile 
reapar.f90 
writesca.f90
