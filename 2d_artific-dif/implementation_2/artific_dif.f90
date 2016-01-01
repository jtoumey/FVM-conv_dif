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
integer :: u_bound,l_bound
real, dimension(:), allocatable :: Su_temp
integer Su_counter
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
!...Calculate fluxes 
!
Fx = rho * u
Fy = rho * v
!
!   calculate coefficients
!
!   allocate space for coefficent vectors
!
allocate(an(np),as(np),aw(np),ae(np),ap(np))
allocate(Su(np),Sp(np))
allocate(phi(np),phi_prev(np))
allocate(Su_temp(ny))
phi = 0.
phi_prev = 0.
!
call calc_fvm_coefficients(np,dx,dy,Fx,Fy,an,as,aw,ae,Su,Sp)
call set_boundary_condition(np,nx,ny,Fx,Fy,dx,dy,ap,an,as,aw,ae,Su,Sp)
call write_coeff_matrix(np,nx,ny,as,aw,ap,ae,an,Su,Sp)
!
!call update_implicit(np,nx,ny,aw,ae,Su,phi_prev)
!
do jj = 1,nx
   Su_counter = 1
   l_bound = (jj - 1)*ny + 1
   u_bound = l_bound + ny - 1
!   if (jj > ny .and. jj < np-ny) then
!   write(*,*)'West Expl: ',aw(jj)*phi_prev(kk+1-ny)
   do kk = l_bound,u_bound
      Su(kk) = Su(kk) + aw(kk)*phi_prev(kk+1-ny) + ae(kk)*phi_prev(kk+ny)
     ! Su_temp(Su_counter) = Su(kk) + aw(kk)*phi_prev(kk+1-ny) + ae(kk)*phi_prev(kk+ny)
     ! Su_counter = Su_counter + 1
!      write(*,*)'Su vector slice: ',Su(l_bound:u_bound)
   end do
 !  end if
   call thomas(ny,-as(l_bound:u_bound),ap(l_bound:u_bound),-an(l_bound:u_bound),Su(l_bound:u_bound),phi(l_bound:u_bound))
!k   call thomas(ny,-as(l_bound:u_bound),ap(l_bound:u_bound),-an(l_bound:u_bound),Su_temp,phi(l_bound:u_bound))
   phi_prev = phi 
 !  call update_implicit(np,nx,ny,aw,ae,Su,phi_prev)
!   write(*,*)'Su after update_implicit: ',Su(np/2)
end do
!--------------------------------------------------------------------------!
!
!...Write the results to a file
!
!--------------------------------------------------------------------------!
open(unit=7,file='phi_distr.dat',ACTION="write", STATUS="replace")
do ii = 1,nx
   do jj = 1,ny
      write(7,301)x(ii),y(jj),phi((ii-1)*ny+jj)
   end do
   write(7,*)
end do
!
!   deallocate data
!
deallocate(x,y)
deallocate(an,as,aw,ae)
deallocate(ap)
!
301 format(3x,f7.2,3x,f7.2,3x,f7.2,3x,f7.2,3x,f7.2,3x,f7.2)
401 format(3x,f12.5,3x,f12.5,3x,f12.5)
!
END
