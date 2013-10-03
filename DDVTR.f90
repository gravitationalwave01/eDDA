      program ddvtr
!
! Program to read polarization and electric field data and output to MAYAVI2 
! or PARAVIEW (VTK/VTR dataformat).
! The code uses "vtr" module written by Jalel Chergui (LIMSI-CNRS)
! VTR is VTK rectilinear format to store XML-type format for rectilinear 
! datasets 
!
! Calling sequence
!  ./ddvtr onemul.bin ddfield
! onemul.bin is a binary polarization file (see "getfml.f90")
! ddfield_1.vtr - VTK ASCII file (for use with "mayavi2")
! ddfield.pvd - ASCII file for use with "paraview"
!
! History:
! Written by PJF June 2011
! Copyright (C) 2011
!               B.T. Draine and P.J. Flatau
! This code is covered by the GNU General Public License.
      USE VTR
      USE DDPRECISION,ONLY: WP
      IMPLICIT NONE
!mayavi2 graphics 
      type(VTR_file_handle) :: fd

!note that mayavi2 graphics fields have to be double precision
      real(kind=8), allocatable, dimension(:) :: x, y, z    ! mesh
      real(kind=8), allocatable, dimension(:,:,:) :: vtr8
      complex(wp),allocatable,dimension(:,:):: cxenat3
      complex(wp):: cxee(3)
      integer::ix,iy,iz

      complex(wp), allocatable, dimension(:):: cxpol,cxeinc, cxe,cxadia 
      integer*2, allocatable, dimension(:):: icomp

! 110809 (BTD) replace ixyz0 with ixyz
!      integer, allocatable, dimension(:,:):: ixyz0
       integer, allocatable, dimension(:,:):: ixyz
!***

      character*25 cflpar, cddfield, command(10)
      
      integer j,jx,jy,jz,nat3,nx,ny,nz, nxyz, nat0

      data cflpar/'onemul.bin'/ 
      data cddfield/'ddfield'/

      if(iargc().eq.1) then 
         CALL GETARG(1 , command(1))
         cflpar=command(1) 
      endif
      print*, ' DDvtr>>  use input binary file ', cflpar
      print*, ' DDvtr>>  use output vtr file ', cddfield

!      OPEN(UNIT=IOBIN,FILE=CFLE1,FORM='UNFORMATTED')
!      WRITE(IOBIN)NXYZ,NAT0,NAT3,NX,NY,NZ
!      WRITE(IOBIN)CXPOL1(1:3*NXYZ),CXEINC(1:3*NXYZ),CXESCA1(1:3*NXYZ)
!      WRITE(IOBIN)CXADIA(1:3*NXYZ),ICOMP(1:NXYZ)
!      CLOSE(IOBIN)

!this file was previously written by "getfml.f90"
      open(66,file=cflpar,form='unformatted') 
      read(66) nxyz, nat0, nat3, nx,ny,nz

      print*, ' nxyz, nat0,nat3, nx,ny,nz ', nxyz, nat0, nat3, nx,ny,nz

!!inconsistent below in nearfield.f90
!         WRITE(IOBIN)NXYZ,NAT3,NX,NY,NZ
!         WRITE(IOBIN)CXPOL2(1:3*NXYZ),CXEINC(1:3*NXYZ),CXESCA2(1:3*NXYZ), &
!                     CXADIA(1:3*NXYZ),ICOMP(1:NXYZ)
!         CLOSE(IOBIN)

      allocate(cxpol(1:3*nxyz),cxeinc(1:3*nxyz),cxe(1:3*nxyz),cxadia(1:3*nxyz), icomp(1:nxyz), ixyz(nxyz,3))
      read(66) cxpol(1:3*nxyz),cxeinc(1:3*nxyz),cxe(1:3*nxyz)
      read(66) cxadia(1:3*nxyz), icomp(1:nxyz)
      close(66)
!calculate total field
      cxe=cxeinc+cxe
!btd could skip storage and reading of ixyz and generate it here:
 
      do jz=1,nz
         do jy=1,ny
            do jx=1,nx
               j=jx+nx*(jy-1)+nx*ny*(jz-1)
               ixyz(j,1)=jx
               ixyz(j,2)=jy
               ixyz(j,3)=jz
            enddo
         enddo
      enddo


!supplementary arrays
      allocate(x(nx),y(ny),z(nz))
      allocate(cxenat3(nxyz,3))
      allocate(vtr8(nx,ny,nz))
!I wish we could reshape to the same array in Fortran90 (or pointer - it is possible in Fortran2003)        
      cxenat3=reshape(cxe,[nxyz,3])

!write VTK file for graphics
!define mesh for graphics
      call VTR_open_file(PREFIX=cddfield, FD=fd)

      do j=1, nxyz
         ix=ixyz(j,1)
         iy=ixyz(j,2)
         iz=ixyz(j,3)
         x(ix)=ix
         y(iy)=iy
         z(iz)=iz
      enddo
      call VTR_write_mesh(FD=fd, X=x, Y=y, Z=z)

!intensity
      do j=1,nxyz
         ix=ixyz(j,1)
         iy=ixyz(j,2)
         iz=ixyz(j,3)
         cxee(1:3)=cxenat3(j,1:3)
         vtr8(ix,iy,iz)=sqrt(sum(abs(cxee)**2))
      enddo
      call VTR_write_var(FD=fd, NAME="Intensity", FIELD=vtr8)

!composition
      do j=1,nxyz
         ix=ixyz(j,1)
         iy=ixyz(j,2)
         iz=ixyz(j,3)
         vtr8(ix,iy,iz)=icomp(j)
      enddo

      call VTR_write_var(FD=fd, NAME="Composition", FIELD=vtr8)
      call VTR_close_file(FD=fd)
      call VTR_collect_file(fd)  ! Produces "ddfield.vtr" and "ddfield.pvd" files
      print*, ' We are done writting vtr file '
      stop
      end program ddvtr
