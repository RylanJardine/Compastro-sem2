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
  real :: x(nx), v(nx), m(nx), rho(nx),  p(nx), cs(nx),h(nx),a(nx), ek,po,mt
  ! real :: u(nx)
  ! define the max and min of your grid and number of particles, n
  real :: dx,xmax,xmin,t,dt,tprint,dtout
  integer :: n,ifile,j,ng
  real,parameter :: pi=4.*atan(1.)


  print*,'Select Standing Wave (1) or Shock Tube Problem (2)'
  read*,j
  if (j==1) then
    ! nx=200
    ! call setup(rho,nx,x,v,xmin,dx,m,cs,n,h,xmax,ng,a,p)
  else
    call setup_shock(rho,nx,x,v,xmin,m,cs,n,h,xmax)
  endif
  ! call to setup initial conditions at time t=0

  ! call setup(rho,nx,x,v,xmin,dx,m,cs,n,h,xmax,ng,a,p)
  call derivs(cs,rho,p,n,a,nx,x,m,h,dx,v,ng)
  open(1,file='kin.dat',status='replace',action='write')
  po=sum(m(1:n)*v(1:n))
  ek=sum(0.5*m(1:n)*v(1:n)**2)
  mt=sum(m(1:n))
  t=0.
  write(1,*)'#t ek po mt'
  write(1,*)t,ek,po,mt


  ! print*,'a',rho

  ifile=0
  dtout=0.01
  tprint=ifile*dtout
  print*,tprint
  call output(n,x,v,h,nx,rho,m,p,cs,a,t,ifile,ek,ng)

  call tim(x,v,a,nx,dt,cs,rho,p,n,m,h,dx,ek,ng)

  ! print*,'b',rho

  do while (t<0.1)
    t=t+dt
    if (t>tprint) then
      ifile=ifile+1
      tprint=ifile*dtout
      print*,tprint
      call output(n,x,v,h,nx,rho,m,p,cs,a,t,ifile,ek,ng)
    end if
    call tim(x,v,a,nx,dt,cs,rho,p,n,m,h,dx,ek,ng)
    po=sum(m(1:n)*v(1:n))
    mt=sum(m(1:n))
    write(1,*)t,ek,po,mt

  enddo
  close(1)

  print*,'hello world'


end program
