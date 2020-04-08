module set2
  !use param for appropriate grid sizing
  use param
  implicit none



contains
  subroutine setup_shock(rho,x,v,m,cs,n,h,a,p,u,z)
    !deine variable types
        integer,intent(inout) :: n
        integer,intent(in) :: z
        real,intent(out) :: x(:),m(:),v(:),h(:),rho(:),cs(:),a(:),p(:),u(:)
        integer :: l,i,g,n_1,n_2

        !define number of particles in each region of shock tube
          n_1=nint(abs(xmax)/dx)
          n_2=nint(abs(xmin2)/dx2)

          !loop from x= -1.5 to -1.
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

          !loop from x = -1 to 0
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


          !loop from x=0.5 to 0
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
          !define number of particles and internal energy
          n=2*n_1+2*n_2
          u(1:n)=1.
          if (z==3) then
            u(1:n)=p(1:n)/((gamma-1)*rho(1:n))
          endif


      end subroutine

      subroutine set_ghosts(rho,nx,x,v,m,cs,n,h,a,p,ng,u,du)
        !define variable types
        integer,intent(in) :: n,nx
        integer,intent(out) ::ng
        real,intent(inout) :: rho(nx),cs(nx),a(nx),p(nx),u(nx),du(nx)
        real,intent(out) :: x(nx),m(nx),v(nx),h(nx)
        real :: l
        integer :: i,j


        l=xmax-xmin
        ng=0
        ! loop to create ghosts on LHS boundary
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
            ! loop to create ghosts on RHS boundary
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
        real,intent(out) :: x(:),m(:),v(:),h(:),rho(:),cs(:),u(:)
        integer :: i
        real,parameter :: pi=4.*atan(1.)

        ! set number of particles, and spacing
        n=100
        dx=(xmax-xmin)/(n)

        !set other parameters
        cs(1:n)=1.
        rho(1:n)=1.
        h(1:n)=1.2*dx
        !loop over posiiton starting half step from boundary
        do i=1,n

          x(i)=xmin+(i-0.5)*dx
          v(i)=cs(i)*10.**(-4)*sin(2*pi*x(i))

        enddo
        !additional parameters
        u(1:n)=1.
        m(1:n)=rho(1:n)*dx


      end subroutine



end module
