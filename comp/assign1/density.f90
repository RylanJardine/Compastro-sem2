module density
  implicit none



contains
  subroutine get_density(m,x,rho,nx,n,h)
    ! set arrays as inputs and outputs of subroutine
    real,intent(in) :: m(nx), x(nx),h(nx)
    real,intent(out) :: rho(nx)
    ! define place holder density, kernel as well as specific smoothing length and particle position
    real :: rhop(nx),w(nx),hin,xin
    ! define integers
    integer :: i,j
    integer,intent(in) :: nx,n


    !call kernel function and to assign density to each particle
    do i=1,nx
      ! call kernel
      call kern(h(i),nx,x,x(i),w,n)
      ! compute the indivdual terms of the density sum
      do j=1,nx
        rhop(j)=m(j)*w(j)
      enddo
      ! sum over all elements to find density per particle
      rho(i)=sum(rhop)
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

  subroutine equation_of_state(cs,rho,p,n)
    ! setup density, pressure and sound speed arrays
    integer,intent(in) :: n
    real,intent(in) :: rho(n)
    real,intent(out) :: p(n), cs(n)

    ! calculate isothermal pressure
    p=cs**2*rho


  end subroutine

  ! subroutine get_accel(cs,rho,p,n,a)
  !   integer,intent(in) :: cs, rho, p, n
  !   real,intent(in) :: rho(n), p(n), cs(n)
  !   real, intent(out) :: a(n)
  !   integer :: i,j
  !
  !   do i=1,n
  !     do j=1,n
  !       ! call kernel
  !       call kern(h(i),nx,x,x(i),w,n)
  !       ! compute the indivdual terms of the density sum
  !       do j=1,nx
  !         rhop(j)=m(j)*w(j)
  !       enddo
  !       ! sum over all elements to find density per particle
  !       rho(i)=sum(rhop)
  !     enddo
  !
  !
  !   enddo
  !
  !
  !   call equation_of_state(cs,rho,p,n)
  !
  ! end subroutine

end module
