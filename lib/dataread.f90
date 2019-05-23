module dataread
 implicit none

contains
 subroutine write_data_phantom() bind(c)

  print*,'hello from Fortran'

 end subroutine write_data_phantom

end module dataread