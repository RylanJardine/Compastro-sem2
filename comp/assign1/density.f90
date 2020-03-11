module density
  implicit none
  public :: dw


contains
  subroutine get_density(m,x,rho,nx,n,h)
    ! set arrays as inputs and outputs of subroutine
    real,intent(in) :: m(nx), x(nx),h(nx)
    real,intent(out) :: rho(nx)
    ! define place holder density, kernel as well as specific smoothing length and particle position
    real :: w(nx),hin,xin
    ! define integers
    integer :: i,j
    integer,intent(in) :: nx,n


    !call kernel function and to assign density to each particle
    do i=1,nx
      ! call kernel
      call kern(h(i),nx,x,x(i),w,n)
      ! sum over all elements to find density per particle
      rho(i)=sum(m*w)
    enddo

  end subroutine



  subroutine kern(hin,nx,x,xin,w,n)
    ! specific smoothing length and specific and general positions
    real,intent(in) :: hin,x(nx),xin
    ! kernel
    real,intent(out) :: w(nx)
    ! set integers
    integer :: i
    integer,intent(in) :: nx,n
    ! define additional parameters
    real :: sig3,q(nx)

    ! set normalising factor for 1D
    sig3=2./(3.*hin)

    ! create array of kernel values
    do i=1,nx
      ! define parameter q
      q(i)=abs(xin-x(i))/hin


      ! set kernel hybrid function
      if (0. .LE. q(i) .and. q(i) <= 1.) then
        w(i)=sig3*(1.-3./2.*q(i)**2.*(1-q(i)/2.))
      elseif (1. < q(i) .and. q(i) <= 2.) then
        w(i)=sig3/4.*(2.-q(i))**3.
      else
        w(i)=0.
      endif

    enddo

  end subroutine


  subroutine equation_of_state(cs,rho,p,nx,n)
    ! setup density, pressure and sound speed arrays
    integer,intent(in) :: nx,n
    real,intent(in) :: rho(nx)
    real,intent(out) :: p(nx), cs(nx)

    ! calculate isothermal pressure

    p(:nx)=cs(:nx)**2*rho(:nx)



  end subroutine

  subroutine get_accel(cs,rho,p,n,a,nx,x,m,h)
    integer,intent(in) :: n,nx
    real,intent(in) :: rho(nx),m(nx)
    real,intent(inout) :: cs(nx), p(nx),x(nx),h(nx)
    real, intent(out) :: a(nx)
    integer :: i,j
    real :: qa,qb,ap(nx)

    qa=0.
    qb=0.
    ! call equation_of_state(cs,rho,p,nx,n)

    do i=1,nx
      ! do j=1,n
        ! call kernel
        ! call kern(h(i),nx,x,x(i),w,n)
        ! compute the indivdual terms of the density sum
        do j=1,nx
          ap(j)=m(j)*((p(i)+qa)/rho(i)**2*dw(x(i),x(j),h(i),j)+(p(j)+qa)/rho(j)**2*dw(x(i),x(j),h(j),j))
        enddo
        print*,'bob'
        ! sum over all elements to find density per particle
        a(i)=-sum(ap)
      enddo

  end subroutine


  real function dw(xin,xon, hin,i)
    real,intent(in) :: xin,xon, hin
    integer,intent(in) :: i
    real :: sig3,q,dq
    ! integer :: i

    sig3=2./(3.*hin)

    ! create array of kernel values
    ! do i=1,nx
      ! define parameter q
      q=abs(xin-xon)/hin
      dq=(xin-xon)/(hin*abs(xin-xon))

      if (xin-xon == 0.) then
        dq=0.
      endif


      ! print*,q,dq
      ! set kernel hybrid function
      if (0. .LE. q .and. q <= 1.) then
        dw=sig3*(-3.*q*dq+9./4.*q**2*dq)
      elseif (1. < q .and. q <= 2.) then
        dw=3.*sig3/4.*(2.-q)**2.*(-dq)
      else
        dw=0.
      endif
      if (dw .ne. 0.) then
        print*,dw,i
      endif
    ! enddo


  end function


  subroutine derivs(cs,rho,p,n,a,nx,x,m,h)
    real, intent(in) ::  m(nx)
    real, intent(inout) :: rho(nx), cs(nx),p(nx), x(nx),h(nx)
    real,intent(out) :: a(nx)
    integer,intent(in) :: n,nx

    call get_density(m,x,rho,nx,n,h)

    call equation_of_state(cs,rho,p,nx,n)

    call get_accel(cs,rho,p,n,a,nx,x,m,h)

  end subroutine

end module
