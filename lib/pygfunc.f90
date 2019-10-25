module write_data_phantom_interface
 use iso_c_binding,      only:c_double,c_float,c_int,c_char  !magic?
 use write_data_phantom, only:write_sphdata_phantom
 implicit none

contains

!! pass MY python things here, through this, and build daniel's things out of them

subroutine c_gfunc(ngas,mgas,x,y,z,h,u,msink,hsink_soft) bind(c)

 !! ACCEPTS THINGS FROM PYTHON
 integer(c_int), intent(in) :: ngas
 real(c_float), intent(in) :: mgas(ngas), x(ngas), y(ngas), z(ngas), h(ngas), u(ngas)
 real(c_float), intent(in) :: msink, hsink_soft

 integer, parameter :: ncolumns = 9 ! number of quantities to write

 integer :: i,ndim,ntotal,ntypes,nsink
 integer :: npartoftype(5),npart
 real(c_float) :: time
 real(c_float) :: gamma

 real(c_float), allocatable  :: dat(:,:)
 real(c_float) :: masstype(5)
 character(len=16) :: label(ncolumns)

 real(c_double) :: udist,umass,utime,umagfd
 character(len=*), parameter :: labeltype(5) = (/'gas ','sink','    ','    ','    '/)

 integer :: ix(3),ivx,ih,iBfirst,ipmass,iutherm
 character(len=120) :: filename

 time = 0.0
 gamma = 5.0/3.0

 print *, "Starting (loc 5)"
 print *, "ngas:", ngas
 print *, "mgas:", mgas(1)
 print *, "x_array", x(1)
 print *, "y_array",y(1)
 print *, "z_array",z(1)
 print *, "h_array",h(1)
 print *, "u_array",u(1)
 print *, "msink", msink!(1)

 umass = 1.98d33
 udist = 6.95d10
 utime = sqrt(udist**3/(6.672041d-8*umass))
 umagfd = 0.

 nsink = 1
 npart = ngas+nsink
 allocate(dat(npart,ncolumns))

 ! BUILD the dat, datsink, label, and label_sink data structures in the form compatible with phantom reader
 ix(1:3) = (/1,2,3/)
 ipmass = 4
 ih = 5
 ivx = 6
 iBfirst = 0
 iutherm = 9

 dat = 0.
 do i=1,ngas
    dat(i,1) = x(i)/udist
    dat(i,2) = y(i)/udist
    dat(i,3) = z(i)/udist
    dat(i,4) = mgas(i)/umass
    dat(i,5) = h(i)/udist
    dat(i,6) = 0. !vx(:)
    dat(i,7) = 0. !vy(:)
    dat(i,8) = 0. !vz(:)
    dat(i,9) = u(i)/(udist/utime)**2
 enddo

 ! sink particle properties
 do i=ngas+1,ngas+nsink
    dat(i,4) = msink/umass
    dat(i,5) = 0.
 enddo

 ntypes = 2
 npartoftype(:) = 0
 masstype(:) = 0.
 npartoftype(1) = ngas
 masstype(1) = mgas(1)
 npartoftype(2) = nsink
 masstype(2) = msink

 label = (/'x ','y ','z ','m ','h ','vx','vy','vz','u '/)
 ndim = 3

 filename = 'star_00000'

 call write_sphdata_phantom(time,gamma,dat,ndim,npart,ntypes,npartoftype, &
                            masstype,ncolumns,udist,umass,utime,umagfd,labeltype,&
                            label,ix,ih,ivx,iBfirst,ipmass,iutherm,filename, hsink_soft)

 deallocate(dat)

end subroutine

end module
