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
!
real :: resid,tol,Frp
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
iter  = 0
resid = 1000.
tol   = 0.01
Frp   = 0.
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
!   Set the fvm coefficients for the entire domain, then overwrite
!   with the boundary conditions
!
call calc_fvm_coefficients(np,dx,dy,Fx,Fy,an,as,aw,ae,Su,Sp)
call set_boundary_condition(np,nx,ny,Fx,Fy,dx,dy,ap,an,as,aw,ae,Su,Sp)
!
!...Begin iteration to solve to within tolerance
!
do while (resid >= tol)
   !
   !...Begin W -> E sweep along each N-S line
   !
   do jj = 1,nx
      !
      ! Calculate the bounds for the current N-S line
      !
      l_bound = (jj - 1)*ny + 1
      u_bound = l_bound + ny - 1
      !
      ! Update Su with the explicit components from the W and E
      !
      call update_explicit(ny,np,l_bound,u_bound,Su_temp,Su,aw,ae,phi_prev)
      !
      !...Solve the tri-diagonal system for a given N-S line
      !
      call thomas(ny,-as(l_bound:u_bound),ap(l_bound:u_bound),-an(l_bound:u_bound),Su_temp,phi(l_bound:u_bound))
      !
      !...Save the solution for explicit treatment at the next N-S line
      !
      phi_prev = phi 
      !
   end do
   !
   !...Calculate residual
   !
   call calc_residual(np,nx,ny,as,aw,ap,ae,an,Su,phi_prev,resid,Frp)
   !
   if (iter < 1) then
      Frp = 1.
   end if
   resid = resid/Frp
   !
   !...Increment the iteration and print the iteration results to the screen
   !
   iter = iter + 1
   write(6,501)iter,resid
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
close(7)
call write_results_diagonal(np,nx,ny,y,phi)
!
!   Deallocate data
!
deallocate(x,y)
deallocate(an,as,aw,ae)
deallocate(ap)
!
301 format(3x,f7.2,3x,f7.2,3x,f7.2,3x,f7.2,3x,f7.2,3x,f7.2)
401 format(3x,f12.5,3x,f12.5,3x,f12.5)
501 format('*** Iteration: ',i6,3x,'Residual: ',f12.5,'   ***')
!
END
