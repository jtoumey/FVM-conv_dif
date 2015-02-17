!**************************************************************************!
!                                                                          !
!  Module:       CONVDIFCDS.F90                                            !
!                                                                          !
!  Programmer:   Julian M. Toumey                                          !
!                Madison, WI                                               !
!                                                                          !
!  Date:         February 2015                                             !
!                                                                          !
!  Language:     FORTRAN90                                                 !
!                                                                          !
!  Description:  This code solves source-free convection-diffusion of a    !
!                property "phi" in a velocity field "u". The method        !
!                follows example 5.1 in Versteeg and Malalasekera, 2nd ed. !
!                The boundary conditions are phi_0 = 1 and phi_L = 0.      !
!                We use the Central Differencing scheme for both the       !
!                convective and diffusive terms.                           !
!                                                                          !
!                Cases                                                     !
!                   (i)   u = 0.1 [m/s], n = 5                             !
!                   (ii)  u = 2.5 [m/s], n = 5                             !
!                   (iii) u = 2.5 [m/s], n = 20                            !
!                                                                          !
!**************************************************************************!
PROGRAM CONVDIFCDS
IMPLICIT NONE
!
integer n,ii,jj
parameter (n=20)
real dx,L,xmax,x(n)
real rho,u,gamma
real phiA,phiB,Dc,F
real aw,ae,Su,Sp,ap
real a(n),b(n),c(n),d(n),phi(n)
real phi_at(n)
!
!...System parameters
!   Density [kg/m^3], velocity [m/s], diffusion coeff [kg/m.s]
!
rho   = 1.
u     = 2.5
gamma = 0.1
phiA  = 1.
phiB  = 0.
!
!...Length [m], dx [m]
!   Diffusion conductance, convective mass flux per unit area
!
L  = 1.
dx = L/float(n)
Dc = gamma/dx
F  = rho*u
!
!...Set up the grid
!
do ii = 1,n
   x(ii) = (ii-0.5)*dx
end do
!
!...Set up system
!   Left boundary
!
aw =  0.
ae =  Dc - F/2.
Su =  (2.*Dc + F)*phiA
Sp = -(2.*Dc + F)
ap =  ae + aw - Sp
!
a(1) = -aw
b(1) =  ap
c(1) = -ae
d(1) =  Su
!
do ii = 2,n-1
   aw = Dc + F/2.
   ae = Dc - F/2.
   Su = 0.
   Sp = 0.
   ap = ae + aw - Sp
   !
   a(ii) = -aw
   b(ii) =  ap
   c(ii) = -ae
   d(ii) =  Su
end do
!...Right boundary
aw =  Dc + F/2.
ae =  0.
Su =  (2.*Dc - F)*phiB
Sp = -(2.*Dc - F)
ap =  ae + aw - Sp
!
a(n) = -aw
b(n) =  ap
c(n) = -ae
d(n) =  Su
!
!...Solve the linear system
!
call thomas(n,a,b,c,d,phi)
!
!...Calculate an analytical solution
!
call analyt_soln(n,x,rho,u,L,gamma,phi_at)
!
!...Write results
!
write(6,101)
open(unit=7,file='temp_distr.dat')
write(6,201)0.,phiA,(1. - exp(0.))/(exp(rho*u*L/gamma)-1) + 1.
write(7,201)0.,phiA,(1. - exp(0.))/(exp(rho*u*L/gamma)-1) + 1.
do jj = 1,n
   write(6,201)x(jj),phi(jj),phi_at(jj)
   write(7,201)x(jj),phi(jj),phi_at(jj)
end do
write(6,201)L,phiB,(exp(1.)-exp(L))/(exp(rho*u*L/gamma)-1)
write(7,201)L,phiB,(exp(1.)-exp(L))/(exp(rho*u*L/gamma)-1)
!
101 format(5x,'____x(j)___',5x,'__phi(j)___',5x,'__Analytic_')
201 format(3x,f12.5,3x,f12.5,3x,f12.5)
END

SUBROUTINE ANALYT_SOLN(n,x,rho,u,L,gamma,phi_at)
integer n,ii
real x(n),phi_at(n)
real rho,u,L,gamma
real d1
!
d1 = exp(rho*u*L/gamma) - 1.
do ii = 1,n
   phi_at(ii) = (1. - exp(rho*u*x(ii)/gamma))/d1 + 1.
end do
END SUBROUTINE ANALYT_SOLN