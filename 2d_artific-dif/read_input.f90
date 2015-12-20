SUBROUTINE READ_INPUT
IMPLICIT NONE
!
character(len=25) :: grid_file_name

grid_file_name = "grid.dat"
write(*,*)grid_file_name

open(unit=2,file = grid_file_name) 
write(*,*)'READING INPUT FILE...'
END SUBROUTINE READ_INPUT
