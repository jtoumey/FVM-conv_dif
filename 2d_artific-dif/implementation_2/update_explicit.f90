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
integer :: check_index_east,check_index_west
!
!...calculate neighbor coefficients using upwind scheme
!
do ii = 1,np
   ! this is the index of the South East cell 
   ! in the SE cell, there is no solution to the East to treat explicitly
   check_index_east = np - ny + 1

   ! this is the index of the South cell one East of the West boundary
   ! in the column of the West boundary, there is no solution to the West 
   ! to treat explicitly
   check_index_west = ny + 1

   !
   if (ii < check_index_west) then
      Su(ii) = Su(ii) + ae(ii)*phi_prev(ii+ny)
   else if (ii < check_index_east) then
      Su(ii) = Su(ii) + aw(ii)*phi_prev(ii+1-ny) + ae(ii)*phi_prev(ii+ny)
   else
      Su(ii) = Su(ii) + aw(ii)*phi_prev(ii+1-ny) 
   end if
end do
!
END SUBROUTINE UPDATE_IMPLICIT
