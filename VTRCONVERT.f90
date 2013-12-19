    program vtrconvert
!!    USE DDPRECISION,ONLY : WP
! 
! Purpose:
! This code converts DDSCAT target shape file to VTK format
! 
! Calling sequence:
! vtrconvert input_ddscat_shape output_vtr_file
! where output_vtr_file is an output file name (without prefix)
!       input_ddscat_shape is file created by "calltarget"
! Thus typical sequance would be:
! calltarget < sphere40x40x40.shp" (to create target.out file - DDSCAT shape file)
! vtrtarget,"target.out","sphere40x40x40" (to do the conversion)
! where the input file "sphere40x40x40.shp"  could be
! ELLIPSOID
! 40 40 40
! 0 0 0
!
! To compile the VTRCONVER:
! "gfortran vtr.f90 vtrconvert.f90 -o vtrconvert"
!
! NOTES:
! You can use "paraview", "mayavi2" to visualize the data using VTK format
!
! History:
! Written by PJF 2011
! (PJF) Feb 2012 Added capability to output a1, a2 vectors to separate file
! Copyright (C) 2011 B.T. Draine and P.J. Flatau
! This code is covered by the GNU General Public License.
!***********************************************************************
     USE VTR
!!     USE DDPRECISION,ONLY : WP
     IMPLICIT NONE
!we can run in in single precision be defult
     INTEGER,PARAMETER :: WP=KIND(0.E0)

!(mayavi2, paraview) graphics 
      type(VTR_file_handle) :: fd


    integer mxnat, ioshp, ipad
    parameter (mxnat=100)
      INTEGER*2 :: ICOMP(MXNAT,3), ICOMP1,ICOMP2,ICOMP3
      INTEGER ::   IXYZ(MXNAT,3)
      integer:: jx, nat, JXX,IXYZ1,IXYZ2,IXYZ3
      integer:: mi1,mi2, mi3, mx1, mx2, mx3
      integer:: mi1ext,mi2ext, mi3ext, mx1ext, mx2ext, mx3ext

      REAL(WP) ::   A1(3), A2(3), DX(3), X0(3)
      REAL(WP) :: AX, AY, AZ
      character*80  cdesc 
!note that mayavi2 graphics fields have to be double precision
      real(kind=8), allocatable, dimension(:) :: x,y,z
      real(kind=8), allocatable, dimension(:,:,:) :: vtr8
      real(kind=8), dimension(1):: xvect, yvect, zvect
      real(kind=8), dimension(1,1,1):: u, v, w
      character*80 cddscat, cvtr
      data cvtr/'target'/
      data cddscat/'target.out'/

    ioshp=10

     if(iargc().eq.2) then 
         CALL GETARG(1 , cddscat)
         CALL GETARG(2 , cvtr)
         print*, ' VTRTARGET>> we will use input DDSCAT file shape format name: ', cddscat(1:20)
         print*, ' VTRTARGET>> we will output VTK name: ', cvtr(1:20)
     else
         print*, ' VTRTARGET>> you have to define 2 files on the command line output '
         print*, ' VTRTARGET>> vtrtarget  ddscat_file vtk_file (without extension)'
         print*, ' VTRTARGET>> for example '
         print*, ' VTRTARGET>> vtrtarget target.out  output'
         stop ' exiting now '
      endif


      OPEN(UNIT=IOSHP,FILE=cddscat,STATUS='UNKNOWN')
       read(ioshp,*) cdesc
       read(ioshp,*) nat
       read(ioshp,*) a1
       read(ioshp,*) a2
       read(ioshp,*) dx
       read(ioshp,*) x0
       read(ioshp,*) cdesc
 
       read(IOSHP,*)JXX,IXYZ1,IXYZ2,IXYZ3,ICOMP1,ICOMP2,ICOMP3
       mi1=ixyz1
       mi2=ixyz2
       mi3=ixyz3
       mx1=ixyz1
       mx2=ixyz2
       mx3=ixyz3

         DO JX=2,NAT
            read(IOSHP,*)JXX,IXYZ1,IXYZ2,IXYZ3,ICOMP1,ICOMP2,ICOMP3
!            print*, JXX,IXYZ1,IXYZ2,IXYZ3,ICOMP1,ICOMP2,ICOMP3
            if(ixyz1.lt.mi1) mi1=ixyz1
            if(ixyz2.lt.mi2) mi2=ixyz2
            if(ixyz3.lt.mi3) mi3=ixyz3

            if(ixyz1.gt.mx1) mx1=ixyz1
            if(ixyz2.gt.mx2) mx2=ixyz2
            if(ixyz3.gt.mx3) mx3=ixyz3

         ENDDO
         print*, mi1, mi2, mi3, mx1, mx2, mx3
         rewind(ioshp)
! allocate array and write
! I am padding to allocate a bit of space around the target so the graphics codes do not
! make objects "holow" (chek ipad=0 to see what I am talking about)
         ipad=5
         mi1ext=mi1-ipad
         mi2ext=mi2-ipad
         mi3ext=mi3-ipad

         mx1ext=mx1+ipad
         mx2ext=mx2+ipad
         mx3ext=mx3+ipad

         allocate(vtr8(mi1ext:mx1ext,mi2ext:mx2ext,mi3ext:mx3ext))
         allocate(x(mi1ext:mx1ext),y(mi2ext:mx2ext), z(mi3ext:mx3ext))

         read(ioshp,*) cdesc
         read(ioshp,*) nat
         read(ioshp,*) a1
         read(ioshp,*) a2
         read(ioshp,*) dx
         read(ioshp,*) x0
         read(ioshp,*) cdesc

         vtr8=0.
         DO JX=1,NAT
            read(IOSHP,FMT=*)JXX,IXYZ1,IXYZ2,IXYZ3,ICOMP1,ICOMP2,ICOMP3
            vtr8(ixyz1,ixyz2,ixyz3)=icomp1
         ENDDO
         CLOSE(UNIT=IOSHP)

        do jx=mi1ext,mx1ext
          x(jx)=jx
        enddo
        do jx=mi2ext,mx2ext
          y(jx)=jx
        enddo
        do jx=mi3ext,mx3ext
          z(jx)=jx
        enddo
!copy to kind=8 precision for graphics
        call VTR_open_file(PREFIX=cvtr, FD=fd)
        call VTR_write_mesh(FD=fd, X=x, Y=y, Z=z)
        call VTR_write_var(FD=fd, NAME="Composition", FIELD=vtr8)
        call VTR_close_file(FD=fd)
        call VTR_collect_file(fd)
!write vectors  a1, a2 to separate files
!copy to kind=8 precision for graphics
        xvect(1)=0.
        yvect(1)=0.
        zvect(1)=0.
        u(1,1,1)=a1(1)
        v(1,1,1)=a1(2)
        w(1,1,1)=a1(3)
        call VTR_open_file(PREFIX='a1a2', FD=fd)
        call VTR_write_mesh(FD=fd, X=xvect, Y=yvect, Z=zvect)
        u(1,1,1)=a1(1)
        v(1,1,1)=a1(2)
        w(1,1,1)=a1(3)
        call VTR_write_var(FD=fd, NAME='a1', VX=u,VY=v,VZ=w)
        u(1,1,1)=a2(1)
        v(1,1,1)=a2(2)
        w(1,1,1)=a2(3)
        call VTR_write_var(FD=fd, NAME='a2', VX=u,VY=v,VZ=w)
        call VTR_close_file(FD=fd)
        call VTR_collect_file(fd)
    stop
    end
