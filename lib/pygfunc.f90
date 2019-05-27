module write_data_phantom_interface
use iso_c_binding, only: c_double, c_int  !magic?
!use phantread_module, only: write_data_phantom
use write_data_phantom_DPmod, only: write_sphdata_phantom
implicit none
contains

!! pass MY python things here, through this, and build daniel's things out of them

subroutine c_gfunc(ngas,mgas,x,y,z,h,u,msink) bind(c) 
	!! ACCEPTS THINGS FROM PYTHON
	!! NOT COMPLETELY WRITTEN YET

    integer(c_int), intent(in) :: ngas 
    real(c_double), dimension(1), intent(in) :: mgas(ngas), x(ngas), y(ngas), z(ngas), h(ngas), u(ngas)
    real(c_double), intent(in) ::  msink


	! integer, parameter :: int8 = selected_int_kind(10)
	! integer, parameter :: sing_prec = c_float
	! integer, parameter :: doub_prec = c_double
	! character(len=10), parameter, public :: formatname='phantom'
	! ! integer, parameter :: lentag = 16

	!  integer(c_int), intent(in)          :: ndim,ntotal,ntypes,ncolumns
	!  integer(c_int), intent(in)          :: npartoftype(:)
	!  real(c_double), intent(in)             :: time,gamma
	!  real(c_double), intent(in)             :: dat(ntotal,ncolumns)
	!  real(c_double), intent(in)             :: masstype(:)
	!  real(c_double), intent(in)  :: udist,umass,utime,umagfd
	!  character(len=*), intent(in) :: labeltype(ntypes),label_dat(ncolumns)
	!  integer(c_int),          intent(in) :: ix(3),ivx,ih,iBfirst,ipmass,iutherm
	!  character(len=*), intent(in) :: filename


    print *, "Starting (loc 5)"
    print *, "ngas:", ngas
    print *, "mgas:", mgas(1)
    print *, "x_array", x(1)
    print *, "y_array",y(1)
    print *, "z_array",z(1)
    print *, "h_array",h(1)
    print *, "u_array",u(1)
    print *, "msink", msink!(1)

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


	! call write_sphdata_phantom(time,gamma,dat,ndim,ntotal,ntypes,npartoftype, &
 !                                 masstype,ncolumns,udist,umass,utime,umagfd,labeltype,&
 !                                 label_dat,ix,ih,ivx,iBfirst,ipmass,iutherm,filename)
	!(n,dat,label,datsink,label_sink)
    print *, "Hello? (loc 6)\n\n"

!    deallocate(dat,datsink)

end subroutine
end module
