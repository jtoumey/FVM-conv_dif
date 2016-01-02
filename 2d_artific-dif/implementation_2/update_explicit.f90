SUBROUTINE UPDATE_EXPLICIT(ny,aw,ae,Su,phi_prev_w,phi_prev_e,Su_temp,jj,nx)
!
implicit none
!
!   variables passed in
integer, intent(inout) :: ny,nx,jj
real, dimension(:), intent(inout) :: aw(ny),ae(ny),phi_prev_w(ny),phi_prev_e(ny)
real, dimension(:), intent(inout) :: Su(ny),Su_temp(ny)
!
!   variables used only in this subroutine
integer ii
integer :: check_east,check_west
!
if (jj == 1) then
   check_west = 0
   check_east = 1
else if (jj == nx) then
   check_west = 1
   check_east = 0
else 
   check_west = 1
   check_east = 1
end if
!
!...calculate neighbor coefficients using upwind scheme
!
Su_temp = Su + aw*phi_prev_w*check_west + ae*phi_prev_e*check_east
!
END SUBROUTINE UPDATE_EXPLICIT
