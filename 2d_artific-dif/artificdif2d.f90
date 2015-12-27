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
!                comparing discretization schemes.                         !
!                                                                          !
!**************************************************************************!
PROGRAM ARTIFICDIF2D
IMPLICIT NONE
!
integer IL,JL,ii,jj,kk,iter
parameter (IL=10,JL=10)
real dx,dy,xmax,ymax,x(IL),y(JL)
real rho,u,v
real phiW,phiN,phiE,phiS,Fx,Fy
real aw,ae,as,an,Su,Sp,ap
real a(JL),b(JL),c(JL),d(JL)
real phisol(JL),resid
real, dimension(JL,IL) :: phi,phiprev
real t1,t2
real Frp
!
!...Parameters for iteration
!
resid = 1000.
iter  = 0.
Frp   = 0.
!
!...Boundary conditions -- temperature [*C]
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
do jj = 1,JL
   y(jj) = (jj-0.5)*dy
end do
!
!...physical properties
!   density [kg/m^3], velocity [m/s], convective mass flux per unit area F
!   plate thickness 1 cm, area [m^2], initial temperature 0 [*C] everywhere
!
rho = 1.
u   = 2.
v   = 2.
!
Fx  = rho*u
Fy  = rho*v
phi = 0.
!--------------------------------------------------------------------------!
!
!...Begin outer WHILE iterative loop 
!
!--------------------------------------------------------------------------!
call cpu_time(t1)
do while (resid >= .001)
   !   save previous phi distribution to compare errors
   phiprev = phi
   !************************************************************************
   !
   !   Set the residual to zero to begin the summation. 
   !
   !************************************************************************
   resid = 0.
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
   write(*,*)'SW aP: ',ap
   !
   a(1) = -as
   b(1) =  ap
   c(1) = -an
   d(1) =  Su
   !
   resid = resid + abs(ae*phi(1,2) + an*phi(2,1) + Su - ap*phi(1,1))
   Frp = Frp + abs(ap*phi(1,1))
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
      !
      resid = resid + abs(ae*phi(jj,2) + an*phi(jj+1,1) + as*phi(jj-1,1) + Su - ap*phi(jj,1))
      Frp = Frp + abs(ap*phi(jj,1))
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
   resid = resid + abs(ae*phi(JL,2) + as*phi(JL-1,1) + Su - ap*phi(JL,1))
   Frp = Frp + abs(ap*phi(JL,1))
   !
   !...solve N-S system with the TDMA
   !
   call thomas(JL,a,b,c,d,phisol)
   !...Store phi solution
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
      resid = resid + abs(aw*phi(1,ii-1) + ae*phi(1,ii+1) + an*phi(2,ii) + Su - ap*phi(1,ii))
      Frp = Frp + abs(ap*phi(1,ii))
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
         !
         resid = resid + abs(aw*phi(jj,ii-1) + ae*phi(jj,ii+1) + an*phi(jj+1,ii) + as*phi(jj-1,ii) + Su - ap*phi(jj,ii))
         Frp = Frp + abs(ap*phi(jj,ii))
         !
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
      resid = resid + abs(aw*phi(JL,ii-1) + ae*phi(JL,ii+1) + as*phi(JL-1,ii) + Su - ap*phi(JL,ii))
      Frp = Frp + abs(ap*phi(JL,ii))
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
   resid = resid + abs(aw*phi(1,IL-1) + an*phi(2,IL) + Su - ap*phi(JL,IL))
   Frp = Frp + abs(ap*phi(JL,IL))
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
      !
      resid = resid + abs(aw*phi(jj,IL-1) + an*phi(jj+1,IL) + as*phi(jj-1,IL) + Su - ap*phi(jj,IL))
      Frp = Frp + abs(ap*phi(jj,IL))
      !
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
   resid = resid + abs(aw*phi(JL,IL-1) + as*phi(JL-1,IL) + Su - ap*phi(JL,IL))
   Frp = Frp + abs(ap*phi(JL,IL))
   !
   !...solve system with the TDMA
   !
   call thomas(JL,a,b,c,d,phisol)
   !...Store N-S temperature solution
   phi(:,IL) = phisol
   !
   !...Recompute the error, increment the iteration
   !
   if (iter == 0) then
      Frp = 1.
   end if
   resid = resid/Frp
   !
   iter  = iter + 1
   write(6,401)iter,resid
end do
call cpu_time(t2)
!--------------------------------------------------------------------------!
!
!...Write the results to a file
!
!--------------------------------------------------------------------------!
open(unit=7,file='phi_distr.dat',ACTION="write", STATUS="replace")
do ii = 1,IL
   do jj = 1,JL
      write(7,301)x(ii),y(jj),phi(jj,ii)
   end do
   write(7,*)
end do
write(6,201)t2 - t1
201 format(3x,f12.5)
301 format(3x,f12.5,3x,f12.5,3x,f12.5)
401 format(3x,'*** Iteration : ',i8,3x,'Residual :',f12.5,'  ***')
END
