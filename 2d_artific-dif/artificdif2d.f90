!**************************************************************************!
!                                                                          !
!  Module:       ARTIFICDIF2D.F90                                          !
!                                                                          !
!  Programmer:   Julian M. Toumey                                          !
!                Madison, WI                                               !
!                                                                          !
!  Date:         February 2015                                             !
!                                                                          !
!  Language:     FORTRAN90                                                 !
!                                                                          !
!  Description:  This code solves source-free pure convection in a 2-D     !
!                domain. The flow direction is 45* from the horizontal     !
!                in order to maximize the artificial diffusion term for    !
!                comparting discretization schemes.                        !
!                                                                          !
!**************************************************************************!
PROGRAM ARTIFICDIF2D
IMPLICIT NONE
!
integer IL,JL,ii,jj,kk,iter
parameter (IL=100,JL=100)
real dx,dy,xmax,ymax,x(IL),y(JL)
real rho,u,v
real phiW,phiN,phiE,phiS,Fx,Fy
real aw,ae,as,an,Su,Sp,ap
real a(JL),b(JL),c(JL),d(JL)
real phisol(JL),resid
real, dimension(JL,IL) :: phi,phiprev
real t1,t2
!
!...Parameters for iteration
!
resid = 1000.
iter  = 0.
!
!...
!
phiW = 100.
phiN = 100.
phiE = 0.
phiS = 0.
!
!...set up the grid
!
xmax = 1
dx   = xmax/float(IL)
do ii = 1,IL
   x(ii) = (ii-0.5)*dx
end do
!
ymax = 1
dy   = ymax/float(JL)
do jj = 2,JL
   y(jj) = (jj-0.5)*dy
end do
!
!...physical properties
!   West heat flux [W/m^2], therm. conduct. [W/m.K], North fixed temp [*C]
!   plate thickness 1 cm, area [m^2], initial temperature 0 [*C] everywhere
!
rho = 1.
u   = 2.
v   = 2. 
!
Fx = rho*u
Fy = rho*v
phi = 0.
!--------------------------------------------------------------------------!
!
!...Begin outer WHILE iterative loop 
!
!--------------------------------------------------------------------------!
call cpu_time(t1)
do while (resid >= .00000001)
   !   save previous temperature distribution to compare errors
   phiprev = phi
   !-----------------------------------------------------------------------!
   !
   !...West Cells
   !
   !-----------------------------------------------------------------------!
   !
   !   SW Corner
   aw = 0.
   ae = 0.
   as = 0.
   an = 0.
   Sp = -(Fx*dy + Fy*dx)
   Su = Fx*dy*phiW + Fy*dx*phiS
   ap = aw + ae + as + an - Sp
   !
   a(1) = -as
   b(1) =  ap
   c(1) = -an
   d(1) =  Su
   !   West Interior cells
   do jj = 2,JL-1
      aw = 0.
      ae = 0.
      as = Fy*dx
      an = 0.
      Sp = -Fx*dy
      Su = Fx*dy*phiW
      ap = aw + ae + as + an - Sp
      !
      a(jj) = -as
      b(jj) =  ap
      c(jj) = -an
      d(jj) =  Su + ae*phi(jj,2)
   end do
   !   NW corner
   aw =  0.
   ae =  0.
   as =  Fy*dx
   an =  0.
   Sp = -Fx*dy
   Su =  Fx*dy*phiW + Fy*dx*phiN
   ap =  aw + ae + as + an - Sp
   !
   a(JL) = -as
   b(JL) =  ap
   c(JL) = -an
   d(JL) =  Su + ae*phi(JL,2)
   !
   !...solve N-S system with the TDMA
   !
   call thomas(JL,a,b,c,d,phisol)
   !...Store temperature solution
   phi(:,1) = phisol
   !-----------------------------------------------------------------------!
   !
   !...Interior Cells, marching W to E
   !
   !-----------------------------------------------------------------------!
   do ii = 2,IL-1
      !   South boundary
      aw = Fx*dy
      ae = 0.
      as = 0.
      an = 0.
      Sp = -Fy*dx
      Su = 0.
      ap = aw + ae + as + an - Sp
      !
      a(1) = -as
      b(1) =  ap
      c(1) = -an
      d(1) =  Su + ae*phi(1,ii+1) + aw*phi(1,ii-1)
      !
      do jj = 2,JL-1
         aw = Fx*dy
         ae = 0.
         as = Fy*dx
         an = 0.
         Sp = 0.
         Su = 0.
         ap = aw + ae + as + an - Sp
         !
         a(jj) = -as
         b(jj) =  ap
         c(jj) = -an
         d(jj) =  Su + ae*phi(jj,ii+1) + aw*phi(jj,ii-1)
      end do
      !...North boundary
      aw =  Fx*dy
      ae =  0.
      as =  Fy*dx
      an =  0.
      Sp =  0.
      Su =  Fy*dx*phiN
      ap = aw + ae + as + an - Sp
      !
      a(JL) = -as
      b(JL) =  ap
      c(JL) = -an
      d(JL) =  Su + ae*phi(jj,ii+1) + aw*phi(jj,ii-1)
      !
      !...solve N-S system with the TDMA
      !
      call thomas(JL,a,b,c,d,phisol)
      !...Store temperature solution
      phi(:,ii) = phisol
   end do
   !-----------------------------------------------------------------------!
   !
   !...East Cells
   !
   !-----------------------------------------------------------------------!
   !
   !   SE Corner
   aw = Fx*dy
   ae = 0.
   as = 0.
   an = 0.
   Sp = -Fy*dx
   Su = Fy*dx*phiS + Fx*dy*phiE
   ap = aw + ae + as + an - Sp
   !
   a(1) = -as
   b(1) =  ap
   c(1) = -an
   d(1) =  Su + aw*phi(1,IL-1)
   !
   do jj = 2,JL-1
      aw = Fx*dy
      ae = 0.
      as = Fy*dx
      an = 0.
      Sp = 0.
      Su = Fx*dy*phiE
      ap = aw + ae + as + an - Sp
      !
      a(jj) = -as
      b(jj) =  ap
      c(jj) = -an
      d(jj) =  Su + aw*phi(jj,IL-1)
   end do
   !...NE corner
   aw =  Fx*dy
   ae =  0.
   as =  Fy*dx
   an =  0.
   Sp =  0.
   Su =  Fx*dy*phiE + Fy*dx*phiN
   ap =  aw + ae + as + an - Sp
   !
   a(JL) = -as
   b(JL) =  ap
   c(JL) = -an
   d(JL) =  Su + aw*phi(JL,IL-1)
   !
   !...solve system with the TDMA
   !
   call thomas(JL,a,b,c,d,phi)
   !...Store N-S temperature solution
   phi(:,IL) = phisol
   !
   !...Recompute the error, increment the iteration
   !
   resid = maxval(abs(phi - phiprev))
   iter  = iter + 1
   write(6,401)iter,resid
end do
call cpu_time(t2)
!--------------------------------------------------------------------------!
!
!...Write the results to a file
!
!--------------------------------------------------------------------------!
open(unit=7,file='plate_temp.dat',ACTION="write", STATUS="replace")
do jj = 1,JL
   write(7,'(1000f12.5)') (phi(jj,ii),ii=1,IL)
end do
write(6,201)t2 - t1
201 format(3x,f12.5)
401 format(3x,'*** Iteration : ',i8,3x,'Residual :',f12.5,'  ***')
END