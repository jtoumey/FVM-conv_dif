SUBROUTINE CREATE_ONED_ARRAYS(np,nx,ny)
use dynamic_coeff
!
IMPLICIT NONE
!
!   variables passed in
integer, intent(inout) :: np,nx,ny
!
!   variables used only in this subroutine
integer ii,jj
!
!   Start outer iteration. This iteration goes from W to E and creates a 1D array for each
!   N-S line of coefficients. The Thomas algorithm can solve each N-S line directly in 
!   O(n) time. Note that W and E coefficients are treated explicity (calculated from \phi
!   at the previous iteration).
!
do ii = 1,nx
   

end do




END SUBROUTINE CREATE_ONED_ARRAYS
