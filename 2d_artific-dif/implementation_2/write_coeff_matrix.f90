SUBROUTINE WRITE_COEFF_MATRIX(np,nx,ny,as,aw,ap,ae,an,Su,Sp)
!
IMPLICIT NONE
!
!   variables passed in
integer, intent(inout) :: np,nx,ny
real, dimension(:), intent(inout) :: as(np),aw(np),ap(np),ae(np),an(np)
real, dimension(:), intent(inout) :: Su(np),Sp(np)
!
!   variables used only in this subroutine
integer :: ii,jj
integer :: coeff_index

open(unit=11,file='coeff_matrix.dat',ACTION="write", STATUS="replace")

!
do jj = 1,nx
   write(11,*)'|    aS    |    aW    |    aP    |    aE    |    aN    |   Su   |   Sp   |'
   write(11,*)'=========================================================================='
   do ii = 1,ny
      coeff_index = (jj - 1)*ny + ii
      write(11,301)as(coeff_index),aw(coeff_index),ap(coeff_index),ae(coeff_index),an(coeff_index),Su(coeff_index),Sp(coeff_index)
   end do
end do
!
301 format(2x,f7.2,4x,f7.2,4x,f7.2,4x,f7.2,4x,f7.2,4x,f7.2,4x,f7.2)
!
END SUBROUTINE WRITE_COEFF_MATRIX
