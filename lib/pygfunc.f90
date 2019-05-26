module write_data_phantom_interface
use iso_c_binding, only: c_double, c_int  !magic?
use phantread_module, only: write_data_phantom
implicit none
contains

!! pass MY python things here, through this, and build daniel's things out of them

subroutine c_gfunc(ngas,mgas,x,y,z,h,u,msink) bind(c)!!x, n, m, a, b, c 

	!! ACCEPTS THINGS FROM PYTHON
	!! NOT COMPLETELY WRITTEN YET


    ! real(c_double), intent(in) :: x
    ! integer(c_int), intent(in) ::  n, m
    ! real(c_double), dimension(n), intent(in) :: a
    ! real(c_double), dimension(m), intent(in) :: b
    ! real(c_double), dimension(n, m), intent(out) :: c

    !integer(c_int), intent(in) :: nsink, ntypes, npartoftype, massoftype
    integer(c_int), intent(in) :: ngas
    integer(c_int), intent(in) :: mgas(ngas)
    real(c_double), intent(in) :: x(ngas), y(ngas), z(ngas), h(ngas), u(ngas), msink(ngas)

    print *, "Starting"
    print *, ngas
    print *, mgas
    print *, x
    print *, y
    print *, z
    print *, h
    print *, u
    print *, msink

    !n = ngas

!    real, allocatable :: dat(:,:), datsink(:,:)

!    nsink = 1
!    allocate(dat(8,n),datsink(6,nsink))
    !! BUILD the dat, datsink, label, and label_sink data structures in the form compatible with phantom reader


!    dat(1,:) = x
!    dat(2,:) = y
!    dat(3,:) = z
!    dat(4,:) = h
!    dat(5,:) = 0. !vx(:)
!    dat(6,:) = 0. !vy(:)
!    dat(7,:) = 0. !vz(:)
!    dat(8,:) = u

    !label = (/'x','y','z','h','u'/)

!    datsink(1,1) = xsink
!    datsink(2,1) = ysink
!    datsink(3,1) = zsink
!    datsink(4,1) = msink
!    datsink(5,1) = rsink
!    datsink(6,1) = hsoftsink

!    label_sink = (/'x','y','z','m','h','hsoft'/)

!    ntypes = 1
    !npartoftype(1) = ngas
    !massoftype(1) = mgas

!    call write_data_phantom(n,dat,label,datsink,label_sink)
    print *, "Hello?"

!    deallocate(dat,datsink)

end subroutine
end module
