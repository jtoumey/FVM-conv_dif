SUBROUTINE CALC_FVM_COEFFICIENTS(np,dx,dy,Fx,Fy,an,as,aw,ae,Su,Sp)
!
IMPLICIT NONE
!
!   variables passed in
integer, intent(inout) :: np
real, intent(inout) :: dx, dy, Fx, Fy
real, dimension(:), intent(inout) :: an(np),as(np),aw(np),ae(np)
real, dimension(:), intent(inout) :: Su(np),Sp(np)
!
!   variables used only in this subroutine
integer ii
!
!...caculate neighbor coefficients using upwind scheme
!
do ii = 1,np
   !
   ae(ii) = 0.
   aw(ii) = Fx*dy
   an(ii) = 0.
   as(ii) = Fy*dx
   Sp(ii) = 0.
   Su(ii) = 0.
   ! 
end do
!
END SUBROUTINE CALC_FVM_COEFFICIENTS
