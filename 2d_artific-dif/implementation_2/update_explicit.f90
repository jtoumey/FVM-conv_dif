SUBROUTINE UPDATE_EXPLICIT(ny,np,l_bound,u_bound,Su_temp,Su,aw,ae,phi_prev)
!
implicit none
!
!   variables passed in
integer, intent(inout) :: ny,np,l_bound,u_bound
real, dimension(:), intent(inout) :: aw(np),ae(np),Su(np),phi_prev(np)
real, dimension(:), intent(inout) :: Su_temp(ny)
!
!
if (u_bound <= ny) then
   !   West-most row: discard West boundary condition contribution
   !   Calculate neighbor coefficients using upwind scheme
   !
   Su_temp = Su(l_bound:u_bound) + ae(l_bound:u_bound)*phi_prev(l_bound+ny:u_bound+ny)
   !
else if (u_bound == np) then
   !   East-most row: discard East boundary condition contribution
   Su_temp = Su(l_bound:u_bound) + aw(l_bound:u_bound)*phi_prev(l_bound-ny:u_bound-ny)
   !
else 
   !   Consider explicit contributions from the West and the East
   !
   Su_temp = Su(l_bound:u_bound) + aw(l_bound:u_bound)*phi_prev(l_bound-ny:u_bound-ny) + &
   ae(l_bound:u_bound)*phi_prev(l_bound+ny:u_bound+ny)
   !
end if
!
END SUBROUTINE UPDATE_EXPLICIT
