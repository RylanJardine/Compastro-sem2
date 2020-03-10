program assign
  ! Import modules for use
  use set
  use inout
  use density
  ! use ghost
  ! use equation
  implicit none
  ! set parameters and arrays of length nx
  integer, parameter :: nx=106
  real :: x(nx), v(nx), m(nx), rho(nx), u(nx), p(nx), cs(nx),h(nx),a(nx)
  ! define the max and min of your grid and number of particles, n
  real :: dx,xmax,xmin
  integer :: n

  ! call to setup initial conditions at time t=0
  call setup(rho,nx,x,v,xmin,dx,m,cs,n,h,xmax)
  ! call to set up ghost particles on either side of grid
  ! call set_ghosts(rho,nx,x,v,xmin,dx,m,cs,n,h,xmax)
  ! call to define density using kernel definition
  call get_density(m,x,rho,nx,n,h)
  ! Calculate pressure/sound speed
  call equation_of_state(cs,rho,p,n)

  ! write output to files
  call output(n,x,v,h,nx,rho,m,p,cs)


  print*,'hello world'


end program
