program assign
  ! Import modules for use
  use set
  use inout
  use density
  use integrator
  ! use ghost
  ! use equation
  implicit none
  ! set parameters and arrays of length nx
  integer, parameter :: nx=120
  real :: x(nx), v(nx), m(nx), rho(nx),  p(nx), cs(nx),h(nx),a(nx)
  ! real :: u(nx)
  ! define the max and min of your grid and number of particles, n
  real :: dx,xmax,xmin,t,dt,tprint,dtout
  integer :: n,ifile
  real,parameter :: pi=4.*atan(1.)


  ! call to setup initial conditions at time t=0
  call setup(rho,nx,x,v,xmin,dx,m,cs,n,h,xmax)

  ! call to set up ghost particles on either side of grid
  ! call set_ghosts(rho,nx,x,v,xmin,dx,m,cs,n,h,xmax)
  ! call to define density using kernel definition
  ! call get_density(m,x,rho,nx,n,h)

  ! Calculate pressure/sound speed
  ! call equation_of_state(cs,rho,p,n)

  ! call get_accel(cs,rho,p,n,a,nx,x,m,h)

  call derivs(cs,rho,p,n,a,nx,x,m,h,dx,v)
  ! print*,'a',rho
  t=0.
  ifile=0
  dtout=0.05
  tprint=ifile*dtout
  call output(n,x,v,h,nx,rho,m,p,cs,a,t,ifile)

  call tim(x,v,a,nx,dt,cs,rho,p,n,m,h,dx)
  ! print*,'b',rho

  do while (t<5.)
    t=t+dt
    if (t>tprint) then
      ifile=ifile+1
      tprint=ifile*dtout
      print*,tprint
      call output(n,x,v,h,nx,rho,m,p,cs,a,t,ifile)
    end if
    call tim(x,v,a,nx,dt,cs,rho,p,n,m,h,dx)
    ! read*,
  enddo

  ! write output to files



  print*,'hello world'


end program
