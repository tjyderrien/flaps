! An example of reading a simple control with labels.
!
! Jason Blevins <jrblevin@sdf.lonestar.org>
! Durham, May 6, 2008

subroutine control_file
  implicit none

  ! Input related variables
  character(len=100) :: buffer, label
  integer :: pos
  integer, parameter :: fh = 15
  integer :: ios = 0
  integer :: line = 0
  real(8) :: pi = acos(-1d0)

  ! Control file variables
  real(8) :: phiMie0
  integer(8) :: PolarizationSource

  open(fh, file='input_Mie.txt')

  ! ios is negative if an end of record condition is encountered or if
  ! an endfile condition was detected.  It is positive if an error was
  ! detected.  ios is zero otherwise.

  do while (ios == 0)
     read(fh, '(A)', iostat=ios) buffer
     if (ios == 0) then
        line = line + 1

        ! Find the first instance of whitespace.  Split label and data.
        pos = scan(buffer, ' 	')
        label = buffer(1:pos)
        buffer = buffer(pos+1:)

        select case (label)
        case ('phiMie0')
           read(buffer, *, iostat=ios) phiMie0
           print *, 'Read phiMie0: ', phiMie0, 'deg -> rad'
           phiMie0=phiMie0*pi/180d0
        case ('PolarizationSource')
           read(buffer, *, iostat=ios) PolarizationSource
           print *, 'Read PolarizationSource TM: ', PolarizationSource
        case default
           print *, 'Skipping invalid label at line', line
        end select
     end if
  end do

end subroutine control_file
