SUBROUTINE SET_BOUNDARY_CONDITION(np,an,as,aw,ae,Su,Sp)
!
IMPLICIT NONE
!
!   variables passed in
integer, intent(inout) :: np
real, dimension(:), intent(inout) :: an(np),as(np),aw(np),ae(np)
real, dimension(:), intent(inout) :: Su(np),Sp(np)
!
!   variables used only in this subroutine
integer ii
!

!
!   West boundary
do ii = 1,



END SUBROUTINE SET_BOUNDARY_CONDITION
