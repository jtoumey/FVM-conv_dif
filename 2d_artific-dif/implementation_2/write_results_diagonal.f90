SUBROUTINE WRITE_RESULTS_DIAGONAL(np,nx,ny,y,phi)
!
implicit none
!
! variables passed in
integer :: np,nx,ny
real, dimension(np), intent(inout) :: phi
real, dimension(ny), intent(inout) :: y
!
! variables used only in this subroutine
integer :: diag_index,ii
!
open(unit=8,file='phi_diagonal.dat',ACTION="write", STATUS="replace")
!
do ii = 1,nx
   diag_index = ii*ny - (ii-1)
   write(8,201)y(ii),phi(diag_index)
end do
!
!   Close the output file
close(8)
!
201 format(3x,f7.2,3x,f7.2)
!
END SUBROUTINE WRITE_RESULTS_DIAGONAL
