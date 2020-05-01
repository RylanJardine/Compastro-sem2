module integrator
  implicit none
  ! set global G and M
  real,parameter :: G=1.,m=1.


contains
  subroutine tim(x,v,a,dt)
    real,intent(inout) :: x(2),v(2)
    real, intent(out) :: a(2)
    real,intent(in) :: dt
    real :: r,a0(2)

        ! Leapfrog integrator
        a0=a

        !update x and v  half-steps
        x=x+dt*v+0.5*(dt)**2*a0
        v=v+dt*a

        !update a
        r=sqrt(x(1)**2+x(2)**2)
        a(:)=-G*M*x(:)/r**3

        !update final step of v and u
        v=v+0.5*dt*(a-a0)

  end subroutine


  subroutine rung(x,v,a,dt)
    real,intent(inout) :: x(2),v(2)
    real, intent(out) :: a(2)
    real,intent(in) :: dt
    real :: r
    real :: k0(2,2),k1(2,2),k2(2,2),k3(2,2)


    ! Fourth Order Runge Kutta

    !Define r and calculate intermediate steps
    r=sqrt(x(1)**2+x(2)**2)
    k0=dt*f1(x,v)
    k1=dt*f1(x+0.5*k0(:,1),v+0.5*k0(:,2))
    k2=dt*f1(x+0.5*k1(:,1),v+0.5*k1(:,2))
    k3=dt*f1(x+k2(:,1),v+k2(:,2))

    ! calculate new position, velocity and accleration
    x=x+1./6.*(k0(:,1)+2.*k1(:,1)+2.*k2(:,1)+k3(:,1))
    v=v+1./6.*(k0(:,2)+2.*k1(:,2)+2.*k2(:,2)+k3(:,2))

    a(:)=-G*M*x(:)/r**3


  end subroutine


    function f1(x,v)
      real,intent(in) :: x(2),v(2)
      real, dimension(2,2) :: f1
      real ::r

      ! Calculate function for intermediate steps and update r
      f1(:,1)=v(:)
      r=sqrt(x(1)**2+x(2)**2)
      f1(:,2)=-G*m*x(:)/r**3

    end function


end module
