!**************************************************************************!
!                                                                          !
!  Module:       ARTIFIC_DIF.F90                                           !
!                                                                          !
!  Programmer:   Julian M. Toumey                                          !
!                Madison, WI                                               !
!                                                                          !
!  Date:         December 2015                                             !
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
real, dimension(:), allocatable :: ap
real, dimension(:), allocatable :: an,as,aw,ae,Su,Sp
real, dimension(:), allocatable :: phi,phi_prev
real dx,dy,xmax,ymax
real rho,u,v
real Fx,Fy
real, dimension(:), allocatable :: x,y
integer np
!
call read_input(xmax,ymax,nx,ny,rho,u,v)
write(*,*)xmax,ymax,nx,ny,rho,u,v
!
!   allocate memory for the coordinates of cell centers
!
allocate(x(nx),y(ny))
!
! setup grid
!
np = nx*ny
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
!
!  Calculate fluxes 
!
Fx = rho * u
Fy = rho * v
!
!   calculate coefficients
!
!   allocate space for coefficent vectors
allocate(an(np),as(np),aw(np),ae(np),ap(np))
allocate(Su(np),Sp(np))
allocate(phi(np),phi_prev(np))
phi = 0.
phi_prev = 0.
!
call calc_fvm_coefficients(np,dx,dy,Fx,Fy,an,as,aw,ae,Su,Sp)
call set_boundary_condition(np,nx,ny,Fx,Fy,dx,dy,ap,an,as,aw,ae,Su,Sp)
!
call update_implicit(np,nx,ny,aw,ae,Su,phi_prev)
do ii = 1,np
   ap(ii) = ae(ii) + aw(ii) + an(ii) + as(ii) - Sp(ii)
end do
!
write(*,*)'|    aS    |    aW    |    aP    |    aE    |    aN    |   Su   |'
write(*,*)'================================================================='
do ii = 1,10
   write(6,301)as(ii),aw(ii),ap(ii),ae(ii),an(ii),Su(ii)
end do
write(*,*)'----//----'
do ii = 11,20
   write(6,301)as(ii),aw(ii),ap(ii),ae(ii),an(ii),Su(ii)
end do
write(*,*)'----//----'
do ii = 21,30
   write(6,301)as(ii),aw(ii),ap(ii),ae(ii),an(ii),Su(ii)
end do
write(*,*)'----//----'
do ii = 31,40
   write(6,301)as(ii),aw(ii),ap(ii),ae(ii),an(ii),Su(ii)
end do
!
call thomas(np,-as,ap,-an,Su,phi)
!end do
!
write(*,*)'Solution:'
do ii = 1,np
   write(*,*)phi(ii)
end do

!
!   deallocate data
!
deallocate(x,y)
deallocate(an,as,aw,ae)
deallocate(ap)
!
301 format(3x,f7.2,3x,f7.2,3x,f7.2,3x,f7.2,3x,f7.2,3x,f7.2)
END
