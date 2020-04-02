module set2
  implicit none


contains
  subroutine setup_shock(rho,nx,x,v,xmin,m,cs,n,h,xmax)

        real, intent(out) :: xmax,xmin
        integer,intent(inout) :: n
        integer,intent(in) :: nx
        ! integer, intent(in) :: nx

        real,intent(out) :: x(:),m(:),v(:),h(:),rho(:),cs(:)
        ! real, intent (inout) :: a(n)
        integer :: l,k,i
        real,parameter :: pi=4.*atan(1.)
        real :: dx, dx2,x2,x3
        ! n=100
        ! print*,'greg'
        xmin=-0.5
        xmax=0.5



        dx=0.001
        dx2=0.01
        ! l=(xmax-xmin)/2/dx
        ! k=(xmax-xmin)/2/dx2
        ! n=l+m

        ! do i=2,l
        !   x(i)=x(i-1)+dx
        !   rho(i)=1.
        !
        ! enddo

        l=1
        x2=xmin
        x(1)=xmin
        rho(1)=1.
        m(1)=dx
        do while (x2<0-dx)

          l=l+1
          x2=x2+dx
          x(l)=x2
          rho(l)=1.
          ! h(l)=1.2*dx
          m(l)=dx

          ! print*,x2
        enddo
        h(1:l)=1.2*dx


        k=1
        x3=0.
        x(l+1)=0.
        rho(l+1)=0.1
        h(l+1)=1.2*dx2
        m(l+1)=dx2
        do while (x3<0.5)
          x3=x3+dx2
          k=k+1
          x(l+k)=x3
          rho(l+k)=0.1

          ! print*,x2
        enddo
        n=l+k
        cs(1:n)=1.
        h(l+1:n)=1.2*dx2
        v(1:n)=cs(1:n)*10.**(-4)*sin(2*pi*x(1:n))
        v=0.
        m(l+1:n)=dx2




        ! do i=l+1,n
        !   x(i)=x(i-1)+dx2
        !   rho(i)=0.1
        ! enddo

        ! do i=1,n
        !   if (x(i)<0) then
        !     rho(i)=1.
        !   else if (x(i)>=0) then
        !     rho(i)=0.1
        !   endif
        ! enddo



        !
        ! dx=(xmax-xmin)/(n-1)
        !
        ! h=1.2*dx
        ! do i=2,n
        !
        !   x(i)=x(i-1)+dx
        !   v(i)=cs(i)*10.**(-4)*sin(2*pi*x(i))
        !
        ! enddo
        !
        ! m=rho*dx

        ! call set_ghosts(rho,nx,x,v,dx,m,cs,n,h)

      end subroutine

      subroutine set_ghosts(rho,nx,x,v,dx,m,cs,n,h,a,p,ng)
        ! real, intent(in) :: dx
        integer,intent(in) :: n,nx
        integer,intent(out) ::ng
        real,intent(inout) :: rho(nx),cs(nx),a(nx),p(nx)
        real,intent(out) :: x(nx),m(nx),v(nx),h(nx)
        real,parameter :: xmax=1.,xmin=0.
        real :: l

        integer :: i
        ! real,parameter :: pi=4.*atan(1.)
        real :: dx,dx2
        dx=0.001
        dx2=0.01
        ! ng=(nx-n)/2
        ! h(n:)=1.2*dx
        ! cs(n:)=1.
        l=xmax-xmin
        ng=12

        ! ng=(nx-n)/2
        ! h(n:)=1.2*dx
        ! cs(n:)=1.

        !
        ! ng=0
        ! do i=1,n
        !   if (x(i)+2*h(i)>xmax) then
        !     ng=ng+1
        !     x(n+ng)=x(i)-(xmax-xmin)
        !   elseif (x(i)-2*h(i)<xmin) then
        !     ng=ng+1
        !     x(n+ng)=x(i)+(xmax-xmin)
        !   endif
        ! enddo
        !
        ! print*,x


        ! the ghost points which lead immediately after the particles e.g. 101, 102 ...
        do i=n+1,n+ng/2
          x(i)=x(i-1)+dx
          ! x(i)=x(i-n+1)+l
          v(i)=v(i-n+1)
          rho(i)=rho(i-n+1)
          m(i)=m(i-n+1)
          cs(i)=cs(i-n+1)
          h(i)=h(i-n+1)
          p(i)=p(i-n+1)
          a(i)=a(i-n+1)
        enddo


        ! the ghost particles which exist before the first particles e.g. -1, -2...
        do i=1,ng/2
          ! x(i+n+ng)=-dx*i
          ! x(i+n+ng/2)=x(n-i)-l
          x(i+n+ng)=x(1)-dx2*i
          v(i+n+ng/2)=v(n-i)
          rho(i+n+ng/2)=rho(n-i)
          m(i+n+ng/2)=m(n-i)
          cs(i+n+ng/2)=cs(n-i)
          h(i+n+ng/2)=h(n-i)
          p(i+n+ng/2)=p(n-i)
          a(i+n+ng/2)=a(n-i)
        enddo
        ! v(n:ng)=0.
        ! print*,h(1:n)


      end subroutine

      ! subroutine set_ghosts(rho,nx,x,v,dx,m,cs,n,h,a,p,ng)
      !   real, intent(in) :: dx
      !   integer,intent(in) :: n,nx
      !   integer,intent(out) ::ng
      !   real,intent(inout) :: rho(nx),cs(nx),a(nx),p(nx)
      !   real,intent(out) :: x(nx),m(nx),v(nx),h(nx)
      !   real,parameter :: xmax=1.,xmin=0.
      !   real :: l
      !
      !   integer :: i
      !   ! real,parameter :: pi=4.*atan(1.)
      !
      !   ! ng=(nx-n)/2
      !   ! h(n:)=1.2*dx
      !   ! cs(n:)=1.
      !   l=xmax-xmin
      !   ng=551
      !
      !   ! ng=(nx-n)/2
      !   ! h(n:)=1.2*dx
      !   ! cs(n:)=1.
      !
      !   !
      !   ! ng=0
      !   ! do i=1,n
      !   !   if (x(i)+2*h(i)>xmax) then
      !   !     ng=ng+1
      !   !     x(n+ng)=x(i)-(xmax-xmin)
      !   !   elseif (x(i)-2*h(i)<xmin) then
      !   !     ng=ng+1
      !   !     x(n+ng)=x(i)+(xmax-xmin)
      !   !   endif
      !   ! enddo
      !   !
      !   ! print*,x
      !
      !
      !   ! the ghost points which lead immediately after the particles e.g. 101, 102 ...
      !   do i=n+1,n+ng/2
      !     ! x(i)=x(i-1)+dx
      !     x(i)=x(i-n+1)+l
      !     v(i)=v(i-n+1)
      !     rho(i)=rho(i-n+1)
      !     m(i)=m(i-n+1)
      !     cs(i)=cs(i-n+1)
      !     h(i)=h(i-n+1)
      !     p(i)=p(i-n+1)
      !     a(i)=a(i-n+1)
      !   enddo
      !
      !
      !   ! the ghost particles which exist before the first particles e.g. -1, -2...
      !   do i=1,ng/2
      !     ! x(i+n+ng)=-dx*i
      !     x(i+n+ng/2)=x(n-i)-l
      !     ! x(i+n+ng)=x(1)-dx*i
      !     v(i+n+ng/2)=v(n-i)
      !     rho(i+n+ng/2)=rho(n-i)
      !     m(i+n+ng/2)=m(n-i)
      !     cs(i+n+ng/2)=cs(n-i)
      !     h(i+n+ng/2)=h(n-i)
      !     p(i+n+ng/2)=p(n-i)
      !     a(i+n+ng/2)=a(n-i)
      !   enddo
      !   ! print*,h(1:n)
      !
      !
      ! end subroutine
      !


end module
