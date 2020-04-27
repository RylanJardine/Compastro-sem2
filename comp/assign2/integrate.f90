module integrator
  implicit none
  real,parameter :: e=0.7,G=1.,m=1.


contains
  subroutine tim(x,v,a,dt)
    real,intent(inout) :: x(2),v(2)
    real, intent(out) :: a(2)
    real,intent(in) :: dt
    real :: r,a0(2)

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
    ! real :: k0(2),k1(2),k2(2),k3(2),l0(2),l1(2),l2(2),l3(2)


    ! r=sqrt(x(1)**2+x(2)**2)
    !a=-gmx/r**3
    !set x1=x, x2=v, thus
    ! dx1/dt=f1(t,x1,x2)=x2=v
    !dx2/dt=dv/dt=-GMx1/r**3
    r=sqrt(x(1)**2+x(2)**2)



    ! k0=dt*f1(x,v)
    ! l0=dt*f2(x,v)
    !
    ! k1=dt*f1(x+0.5*k0,v+0.5*l0)
    ! l1=dt*f2(x+0.5*k0,v+0.5*l0)
    !
    ! k2=dt*f1(x+0.5*k1,v+0.5*l1)
    ! l2=dt*f2(x+0.5*k1,v+0.5*l1)
    !
    ! k3=dt*f1(x+k2,v+l2)
    ! l3=dt*f2(x+k2,v+l2)
    !
    !
    ! x=x+1./6.*(k0+2.*k1+2.*k2+k3)
    ! v=v+1./6.*(l0+2.*l1+2.*l2+l3)
    !
    ! a(1)=-G*M*x(1)/r**3
    ! a(2)=-G*M*x(2)/r**3

    r=sqrt(x(1)**2+x(2)**2)
    k0=dt*f1(x,v)


    k1=dt*f1(x+0.5*k0(:,1),v+0.5*k0(:,2))


    k2=dt*f1(x+0.5*k1(:,1),v+0.5*k1(:,2))


    k3=dt*f1(x+k2(:,1),v+k2(:,2))



    x=x+1./6.*(k0(:,1)+2.*k1(:,1)+2.*k2(:,1)+k3(:,1))
    v=v+1./6.*(k0(:,2)+2.*k1(:,2)+2.*k2(:,2)+k3(:,2))

    a(:)=-G*M*x(:)/r**3


  end subroutine


    function f1(x,v)
      real,intent(in) :: x(2),v(2)
      real, dimension(2,2) :: f1
      real ::r


      ! f1(:)=v(:)

      f1(:,1)=v(:)
      r=sqrt(x(1)**2+x(2)**2)
      f1(:,2)=-G*m*x(:)/r**3

    end function

  !
  !   function f2(x,v)
  !     real,intent(in) :: x(2),v(2)
  !     ! real,intent(in) :: r
  !     real, dimension(2) :: f2
  !     real :: r
  !
  !     r=sqrt(x(1)**2+x(2)**2)
  !     f2(:)=-G*m*x(:)/r**3
  !
  ! end function

end module
