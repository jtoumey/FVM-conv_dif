SUBROUTINE UPDATE_IMPLICIT(np,nx,ny,aw,ae,Su,phi_prev)
!
implicit none
!
!   variables passed in
integer, intent(inout) :: np,nx,ny
real, dimension(:), intent(inout) :: aw(np),ae(np),Su(np),phi_prev(np)
!
!   variables used only in this subroutine
integer ii
!
!...caculate neighbor coefficients using upwind scheme
!
do ii = ny,np-ny
   Su(ii) = Su(ii) + aw(ii-ny)*phi_prev(ii-ny) + ae(ii+ny)*phi_prev(ii+ny)
end do

!
END SUBROUTINE UPDATE_IMPLICIT
