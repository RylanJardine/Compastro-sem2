program assign
  ! Import modules for use
  ! use set
  use set2
  use inout
  use density
  use integrator
  ! use ghost
  ! use equation
  implicit none
  ! set parameters and arrays of length nx
  integer, parameter :: nx=1500
  real :: x(nx), v(nx), m(nx), rho(nx),  p(nx), cs(nx),h(nx),a(nx), ek,po,mt,du(nx),u(nx)
  ! real :: u(nx)
  ! define the max and min of your grid and number of particles, n
  real :: xmax,xmin,t,dt,tprint,dtout,b
  integer :: n,ifile,z,ng,i
  real,parameter :: pi=4.*atan(1.)




  !DO not ever use my global integers z and y


  print*,'Select Standing Wave (1), Isothermal Shock tube (2) or Adiabatic Shock Tube Problem (3)'
  read*,z
  if (z==1) then
    ! call setup(rho,nx,x,v,xmin,dx,m,cs,n,h,xmax,ng,a,p)
  else if ((z==2) .or. (z==3)) then
    call setup_shock(rho,nx,x,v,xmin,m,cs,n,h,xmax,a,p,u,z)
  else
    print*, 'Please Select (1),(2) or (3)'
  endif

  ! print*, 'Add viscosity? (yes/no)'
  ! read*,b
  ! if ((b==yes) .or. (b==y)) then
  !
  ! else if ((b==no) .or. b==(n)) then
  !
  ! else
  !   print*, 'Default without viscosity'
  !
  ! endif

  ! call to setup initial conditions at time t=0

  ! call setup(rho,nx,x,v,xmin,dx,m,cs,n,h,xmax,ng,a,p)
  open(2,file='try.dat',status='replace',action='write')
  write(2,*)'# x rho'
  do i=1,n
    write(2,*)x(i),rho(i)
  enddo
  close(2)



  call derivs(cs,rho,p,n,a,nx,x,m,h,dx,v,ng,du,u)
  ! open(2,file='try.dat',status='replace',action='write')
  ! write(2,*)'# x rho'
  ! do i=1,n+ng
  !   write(2,*)x(i),rho(i)
  ! enddo
  ! close(2)
  ! read*,
  open(1,file='kin.dat',status='replace',action='write')
  po=sum(m(1:n)*a(1:n))
  ek=sum(0.5*m(1:n)*v(1:n)**2)
  mt=sum(m(1:n))
  t=0.
  write(1,*)'#t ek po mt'
  write(1,*)t,ek,po,mt


  ifile=0
  dtout=0.01
  tprint=ifile*dtout
  print*,tprint
  call output(n,x,v,h,nx,rho,m,p,cs,a,t,ifile,ng,u,du)

  call tim(x,v,a,nx,dt,cs,rho,p,n,m,h,dx,ng,du,u)

  ! print*,'b',rho

  do while (t<0.3)
    t=t+dt
    if (t>tprint) then
      ifile=ifile+1
      tprint=ifile*dtout
      print*,tprint
      call output(n,x,v,h,nx,rho,m,p,cs,a,t,ifile,ng,u,du)
    end if
    call tim(x,v,a,nx,dt,cs,rho,p,n,m,h,dx,ng,du,u)
    po=sum(m(1:n)*a(1:n))
    mt=sum(m(1:n))
    ek=sum(0.5*m(1:n)*v(1:n)**2)
    if (ek/ek .ne. 1) then
      call output(n,x,v,h,nx,rho,m,p,cs,a,t,ifile,ng,u,du)
      read*,
    endif
    write(1,*)t,ek,po,mt
    ! print*,a(1)

  enddo
  close(1)

  print*,'hello world'


end program
