SUBROUTINE WRITE_RESULTS_DIAGONAL(np,ny,y,phi)
!
implicit none
!
! variables passed in
integer :: np,ny
real, dimension(np), intent(inout) :: phi
real, dimension(ny), intent(inout) :: y
!
! variables used only in this subroutine
integer :: diag_index,grid_index
!
diag_index = ny
grid_index = 1
!
open(unit=8,file='phi_diagonal.dat',ACTION="write", STATUS="replace")
!
do while (diag_index > 0)
   write(8,201)y(grid_index),phi(diag_index)
   diag_index = diag_index - 1
   grid_index = grid_index + 1
end do
!
!   Close the output file
close(8)
!
201 format(3x,f7.2,3x,f7.2)
!
END SUBROUTINE WRITE_RESULTS_DIAGONAL
