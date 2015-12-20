!**************************************************************************!
!                                                                          !
!  Module:       ARTIFIC_DIF.F90                                           !
!                                                                          !
!  Programmer:   Julian M. Toumey                                          !
!                Madison, WI                                               !
!                                                                          !
!  Date:         Decemver 2015                                             !
!                                                                          !
!  Language:     FORTRAN90                                                 !
!                                                                          !
!  Description:  This code solves source-free pure convection in a 2-D     !
!                domain. The flow direction is 45* from the horizontal     !
!                in order to maximize the artificial diffusion term for    !
!                comparing discretization schemes.                         !
!                                                                          !
!                This program differs from ARTIFICDIF2D in that it treats  !
!                boundary conditions differently and allows either upwind  !
!                or central differencing for the convection term.          !
!                                                                          !
!**************************************************************************!
PROGRAM ARTIFIC_DIF
IMPLICIT NONE
!
integer nx,ny,ii,jj,kk,iter
real dx,dy,xmax,ymax
real, dimension(:), allocatable :: x,y
!
call read_input(xmax,ymax,nx,ny)
write(*,*)xmax,ymax,nx,ny
!
!   allocate memory for the coordinates of cell centers
!
allocate(x(nx),y(ny))
!
! setup grid
!
dx = xmax/float(nx)
do ii = 1,nx
   x(ii) = (ii - 0.5)*dx
end do
!
dy = ymax/float(ny)
do jj = 1,ny
   y(jj) = (jj - 0.5)*dy
end do

!
deallocate(x,y)
!
END
