module integrator
  use physics
  implicit none

contains
  subroutine tim(x,v,a,nx,dt,cs,rho,p,n,m,h,ng,du,u)
    integer,intent(in) :: nx,n
    real, intent(inout) :: x(nx),v(nx),a(nx),h(nx),rho(nx),m(nx),p(nx),cs(nx),du(nx),u(nx)
    real :: dt,a0(nx),du0(nx)
    integer,intent(inout) :: ng
    !set timestep with stability condition
    dt=0.2*minval(h(1:n)/cs(1:n))
    !define initial du and a
    a0=a
    du0=du

    !update x and v and u half-steps
    x(1:n)=x(1:n)+dt*v(1:n)+0.5*(dt)**2*a0(1:n)

    v(1:n)=v(1:n)+dt*a(1:n)
    u(1:n)=u(1:n)+dt*du(1:n)
    !update a and du
    call derivs(cs,rho,p,n,a,nx,x,m,h,v,ng,du,u)

    !update final step of v and u
    v(1:n)=v(1:n)+0.5*dt*(a(1:n)-a0(1:n))
    u(1:n)=u(1:n)+0.5*dt*(du(1:n)-du0(1:n))


  end subroutine
end module
