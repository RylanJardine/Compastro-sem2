program assign
  use set
  use inout
  implicit none
  integer, parameter :: nx=106
  real :: x(nx), v(nx), m(nx), rho(nx), u(nx), p(nx), cs(nx),h(nx)
  real :: dx
  real,parameter :: xmax=1.,xmin=0.,rho0=1.
  integer :: n

  rho(:)=rho0
  cs(:)=1.



  call setup(rho,nx,x,v,xmin,dx,m,rho0,cs,n,h,xmax)


  call output(n,x,v,h,nx)


  print*,'hello world'

  print*,v,x

    print*,n



end program
