SUBROUTINE CALC_FVM_COEFFICIENTS(np,nx,ny,dx,dy,Fx,Fy)
use dynamic_coeff
!
IMPLICIT NONE
!
! variables passed in
integer, intent(inout) :: np, nx, ny
real, intent(inout) :: dx, dy, Fx, Fy
!
! variables used only in this subroutine
integer ii
!
do ii = 1,np
   an(ii + nx)


end do
do ii = 1,np
   aw(ii) = Fx*dy
   an(ii) = Fy*dx
   ! 
   
end do
!
END SUBROUTINE CALC_FVM_COEFFICIENTS
