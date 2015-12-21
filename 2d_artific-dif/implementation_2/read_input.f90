SUBROUTINE READ_INPUT(xmax,ymax,nx,ny,rho,u,v)
IMPLICIT NONE
!
! variables passed in, to be modified
real, intent(out) :: xmax, ymax
integer, intent(out) :: nx, ny
real, intent(out) :: rho, u, v
!
! variables used only in this subroutine
character(len=25 ) :: grid_file_name
character(len=100) :: line_buffer
integer iblnk,io_status
!
io_status = 1
grid_file_name = "input.dat"
write(*,*)grid_file_name

open(unit=2,file = grid_file_name) 
write(*,*)'READING INPUT FILE...'

do while (io_status .GT. 0)
   ! read each line, save to variable
   read(2,*,IOSTAT = io_status)line_buffer, xmax
   read(2,*)line_buffer, ymax
   read(2,*)line_buffer, nx 
   read(2,*)line_buffer, ny 
   read(2,*)line_buffer, rho
   read(2,*)line_buffer, u 
   read(2,*)line_buffer, v 
   !
end do
write(*,*)'SUCCESSFULLY READ THE INPUT FILE.'
! Figure this out later
!do
!   read(2,*,IOSTAT = io_status)line_buffer

!   if (io_status == 0) then 

!      write(*,*)xmax
!   else if (io_status < 0) then
!      write(*,*)'END OF FILE.'
!      exit
!   end if
!end do 
!read(2,*)ymax
!read(2,*)nx
!read(2,*)ny

close(2)
!
END SUBROUTINE READ_INPUT
