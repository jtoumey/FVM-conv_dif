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
!                                                                          !
!**************************************************************************!
PROGRAM CONVDIFCDS
IMPLICIT NONE
!
integer n,ii,
parameter (n=5)
real dx,L,xmax,x(n)
real rho,u,gamma
real phiA,phiB
!
!
!
L = 1.
rho = 1.
u = 0.1
gamma = 0.1
phiA = 1.
phiB = 0.


D = gamma/dx
F = rho*u


dx = L/float(n)
do ii = 1,n
   x(ii) = (ii-0.5)*dx
end do
!
!...Left boundary
!
aw =  0.
ae =  D - F/2.
Su =  (2.*D + F)*phiA
Sp = -(2.*D + F)
ap =  ae + aw - Sp
!
a(1) = -aw
b(1) =  ap
c(1) = -ae
d(1) =  Su
!
do ii = 2,n-1

end do
!...Right boundary

END