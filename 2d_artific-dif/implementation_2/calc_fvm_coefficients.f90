SUBROUTINE CALC_FVM_COEFFICIENTS(np,an,as,aw,ae)
!
IMPLICIT NONE
!
! variables passed in
integer, intent(out) :: np
real, intent(inout) :: an(:),as(:),aw(:),ae(:)
!
! variables used only in this subroutine
integer ii
do ii = 1,np
   an(ii) = ii*2 
end do
!
END SUBROUTINE CALC_FVM_COEFFICIENTS
