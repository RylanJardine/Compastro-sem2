module integrator
  use density
  implicit none

contains
  subroutine tim(x,v,a,t,nx,dt,cs,rho,p,n,m,h,dx)
    real, intent(inout) :: x(nx),v(nx),a(nx),h(nx),rho(nx),m(nx),p(nx),cs(nx),dx
    integer,intent(in) :: nx,n
    real, intent(out) :: t
    real :: vs(nx),ab(nx),dt,a0(nx)
    dt=0.001
    a0=a
    x(:n)=x(:n)+dt*v(:n)+0.5*(dt)**2*a(:n)

    vs(:n)=v(:n)+dt*a(:n)
    call derivs(cs,rho,p,n,a,nx,x,m,h,dx,v)
    print*,
    ! v(:)=vs(:)+0.5*dt(a(:)-ab(:))
    ! v(:n)=vs(:n)+0.5*dt*(a(:n)-a0(:n))
    v(:n)=v(:n)+dt*a(:n)+0.5*dt*(a(:n)-a0(:n))
    print*,(a(1)-a0(1))
    print*,v(1),a(1),a0(1),dt
    a=a0
    ! v(1)=0.
    ! v(n)=0.


  end subroutine
end module
