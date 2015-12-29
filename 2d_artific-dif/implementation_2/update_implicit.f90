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
!...calculate neighbor coefficients using upwind scheme
!
do ii = ny,np-ny
   write(*,301)Su(ii),aw(ii+1-ny),phi_prev(ii+1-ny),ae(ii+ny),phi_prev(ii+ny)
   Su(ii) = Su(ii) + aw(ii)*phi_prev(ii+1-ny) + ae(ii)*phi_prev(ii+ny)
end do
!
301 format(3x,f7.2,3x,f7.2,3x,f7.2,3x,f7.2,3x,f7.2,3x)
!
END SUBROUTINE UPDATE_IMPLICIT
