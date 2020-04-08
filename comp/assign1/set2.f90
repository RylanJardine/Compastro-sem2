module set2
  use param
  implicit none
  ! real :: dx,dx2
  ! real :: gamma,xmax,xmin
  ! integer :: y


contains
  subroutine setup_shock(rho,x,v,m,cs,n,h,a,p,u,z)

        integer,intent(inout) :: n
        integer,intent(in) :: z
        real,intent(out) :: x(:),m(:),v(:),h(:),rho(:),cs(:),a(:),p(:),u(:)
        integer :: l,i,g,n_1,n_2
        ! real :: xmid,xmid2,xmin2
        !
        !
        !
        !
        !   xmin=-1.5
        !   xmid=-1.
        !   xmid2=0.
        !   xmax=0.5
        !   xmin2=-0.5
          n_1=nint(abs(xmax)/dx)
          n_2=nint(abs(xmin2)/dx2)


          do i=1,n_2
            x(i)=xmin+(i-0.5)*dx2
            v(i) = 0.
            m(i) = dx
            h(i) = 1.2*dx2
            rho(i) = 0.125
            a(i) = 0.
            cs(i) = 1.
            p(i) =0.1
          enddo

          l=n_2
          do i=1,2*n_1
            x(l+i)=xmid2-(i-0.5)*dx
              v(l+i) = 0.
              m(l+i) = dx
              h(l+i) = 1.2*dx
              cs(l+i) = 1.
              a(l+i) = 0.
              p(l+i) = 1.
              rho(l+i) = 1.
          enddo

          g=2*n_1+n_2
          do i=1,n_2
            x(g+i)=xmax-(i-0.5)*dx2
            v(g+i) = 0.
            m(g+i) = dx
            h(g+i) = 1.2*dx2
            rho(g+i) = 0.125
            a(g+i) = 0.
            cs(g+i) = 1.
            p(g+i) =0.1
          enddo
          n=2*n_1+2*n_2
          u(1:n)=1.
          if (z==3) then
            u(1:n)=p(1:n)/((gamma-1)*rho(1:n))
          endif


      end subroutine

      subroutine set_ghosts(rho,nx,x,v,m,cs,n,h,a,p,ng,u,du)
        integer,intent(in) :: n,nx
        integer,intent(out) ::ng
        real,intent(inout) :: rho(nx),cs(nx),a(nx),p(nx),u(nx),du(nx)
        real,intent(out) :: x(nx),m(nx),v(nx),h(nx)
        real :: l
        integer :: i,j


        l=xmax-xmin
        ng=0
        do i=1,n
          if (x(i)+2.*h(i)>xmax) then
            ng=ng+1
            j=n+ng
            x(j)=x(i)-(xmax-xmin)
            v(j)=v(i)
            h(j)=h(i)
            m(j)=m(i)
            p(j)=p(i)
            rho(j)=rho(i)
            a(j)=a(i)
            cs(j)=cs(i)
            u(j)=u(i)
            du(j)=du(i)
          elseif (x(i)-2.*h(i)<xmin) then
            ng=ng+1
            j=n+ng
            x(j)=x(i)+(xmax-xmin)
            v(j)=v(i)
            h(j)=h(i)
            m(j)=m(i)
            p(j)=p(i)
            rho(j)=rho(i)
            a(j)=a(i)
            cs(j)=cs(i)
            u(j)=u(i)
            du(j)=du(i)
          endif
        enddo

      end subroutine


      subroutine setup(rho,x,v,m,cs,n,h,u)

        integer,intent(inout) :: n
        ! integer,intent(in) :: z
        real,intent(out) :: x(:),m(:),v(:),h(:),rho(:),cs(:),u(:)
        integer :: i
        real,parameter :: pi=4.*atan(1.)
        n=100
        ! y=z
        ! gamma=1.
        !
        ! xmin=0.
        ! xmax=1.
        cs(1:n)=1.
        rho(1:n)=1.


        dx=(xmax-xmin)/(n)

        h(1:n)=1.2*dx
        do i=1,n

          x(i)=xmin+(i-0.5)*dx
          v(i)=cs(i)*10.**(-4)*sin(2*pi*x(i))

        enddo
        u(1:n)=1.
        m(1:n)=rho(1:n)*dx


      end subroutine



end module
