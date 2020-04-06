module set2
  implicit none


contains
  subroutine setup_shock(rho,nx,x,v,xmin,m,cs,n,h,xmax,a,p)

        real, intent(out) :: xmax,xmin
        integer,intent(inout) :: n
        integer,intent(in) :: nx
        real,intent(out) :: x(:),m(:),v(:),h(:),rho(:),cs(:),a(:),p(:)
        integer :: l,k,i,g
        real,parameter :: pi=4.*atan(1.)
        real :: dx, dx2,x2,x3,x1,xmid,xmid2,len
        !
        ! xmin=-1.5
        ! xmax=0.5
        !
        !
        ! dx=0.001
        ! dx2=0.01
        !
        ! x1=xmin+dx
        ! ! x(1)=xmin
        !
        !
        ! g=0
        !
        ! ! x(l+1)=0.
        ! ! rho(1)=0.1
        !
        ! do while (x1<-1.)
        !
        !   g=g+1
        !   x(g)=x1
        !   rho(g)=0.1
        !   x1=x1+dx2
        ! enddo
        ! h(1:g)=1.2*dx2
        !
        ! x2=-1.
        ! ! rho(g+1)=1.
        ! l=0
        ! ! x(g+1)=x2
        ! do while (x2<0)
        !
        !   l=l+1
        !
        !   x(l+g)=x2
        !   rho(l+g)=1.
        !   x2=x2+dx
        ! enddo
        ! h(g+1:l)=1.2*dx
        !
        !
        ! k=0
        ! x3=0.
        ! ! x(l+1+g)=0.
        ! ! rho(l+1+g)=0.1
        ! ! h(l+1+g)=1.2*dx2
        ! do while (x3<0.5+dx2)
        !
        !   k=k+1
        !   x(l+k+g)=x3
        !   rho(l+k+g)=0.1
        !   x3=x3+dx2
        ! enddo
        ! n=l+k+g
        ! cs(1:n)=1.
        ! h(l+1:n)=1.2*dx2
        ! ! v(1:n)=cs(1:n)*10.**(-4)*sin(2*pi*x(1:n))
        ! v=0.
        ! m(1:n)=rho(1)*dx2


        ! h(1:n)=maxval(1.2*m(1:n)/rho(1:n))


          ! integer ::  n_l, n_r
          ! real :: rho_l, rho_r, dx_l, dx_r, b
          ! real :: middle
          !
          ! rho_l = 1.
          ! rho_r = 0.1
          ! dx_l = 0.001
          ! dx_r = 0.01
          ! xmin = -0.5
          ! xmax = 0.5
          ! middle = 0
          ! n_l = nint(abs(xmax)/dx_l)  !500
          ! n_r = nint(abs(xmin)/dx_r)  !50
          ! n = 2*n_l + 2*n_r
          ! b= middle
          !
          !
          !
          !
          ! do i = 1,n_l
          !   !for left hand side start at -0.5 and go to 0 (500)
          !   b = xmin + (i-0.5)*dx_l
          !   x(i) = b
          !   v(i) = 0.
          !   m(i) = dx_l
          !   h(i) = 1.2*dx_l
          !   cs(i) = 1.
          !   a(i) = 0.
          !   p(i) = 1.
          !   rho(i) = 1.
          ! enddo
          ! b = middle
          ! do i = n_l+1, n_l+n_r
          !   !for right hand side start at 0 and go to 0.5 (50)
          !   b = middle + ((i - n_l)-0.5)*dx_r
          !   x(i) = b
          !   v(i) = 0.
          !   m(i) = dx_l
          !   h(i) = 1.2*dx_r
          !   rho(i) = 0.1
          !   a(i) = 0.
          !   cs(i) = 1.
          !   p(i) =0.1
          !
          ! enddo
          ! !now so we can simplify the ghost part
          ! !lets reflect what we have done above around -0.5 (500)
          ! do i = n_l+n_r+1, 2*n_l+n_r
          !   !set the high density part from -0.5 to -1
          !   b = x(1) - ((i - (n_l+n_r)))*dx_l
          !   x(i) = b
          !   v(i) = 0.
          !   m(i) = dx_l
          !   h(i) = 1.2*dx_l
          !   rho(i) = 1.
          !   a(i) = 0.
          !   cs(i) = 1.
          !   p(i) =1.
          ! enddo
          ! do i = 2*n_l+n_r+1, 2*n_l+2*n_r
          !   !set the low density part from -1 to -1.5 (50)
          !   b = x(2*n_l+n_r) - ((i - (2*n_l+n_r)))*dx_r
          !   x(i) = b
          !   v(i) = 0.
          !   m(i) = dx_l
          !   h(i) = 1.2*dx_r
          !   rho(i) = 0.1
          !   a(i) = 0.
          !   cs(i) = 1.
          !   p(i) =0.1
          ! enddo

          n=1100
          xmin=-1.5
          xmid=-1.
          xmid2=0.
          xmax=0.5
          len=xmax-xmin
          dx=0.001
          dx2=0.01

          do i=1,50
            x(i)=xmin+(i-0.5)*dx2
            v(i) = 0.
            m(i) = dx
            h(i) = 1.2*dx2
            rho(i) = 0.1
            a(i) = 0.
            cs(i) = 1.
            p(i) =0.1
          enddo

          l=50
          do i=1,1000
            x(l+i)=xmid+(i-0.5)*dx
              v(l+i) = 0.
              m(l+i) = dx
              h(l+i) = 1.2*dx
              cs(l+i) = 1.
              a(l+i) = 0.
              p(l+i) = 1.
              rho(l+i) = 1.
          enddo

          g=1050
          do i=1,50
            x(g+i)=xmid2+(i-0.5)*dx2
            v(g+i) = 0.
            m(g+i) = dx
            h(g+i) = 1.2*dx2
            rho(g+i) = 0.1
            a(g+i) = 0.
            cs(g+i) = 1.
            p(g+i) =0.1
          enddo

          ! print*,x(1:n)
          ! read*,



      end subroutine

      subroutine set_ghosts(rho,nx,x,v,dx,m,cs,n,h,a,p,ng)
        ! real, intent(in) :: dx
        integer,intent(in) :: n,nx
        integer,intent(out) ::ng
        real,intent(inout) :: rho(nx),cs(nx),a(nx),p(nx)
        real,intent(out) :: x(nx),m(nx),v(nx),h(nx)
        real,parameter :: xmax=0.5,xmin=-1.5
        real :: l

        integer :: i,j

        real :: dx,dx2
        dx=0.001
        dx2=0.01
        !
        l=xmax-xmin
        ng=20
        !
        !
        ! ! the ghost points which lead immediately after the particles e.g. 101, 102 ...
        ! do i=n+1,n+ng/2
        !   ! x(i)=x(i-1)+dx2
        !   x(i)=x(i-n+1)+l
        !   v(i)=v(i-n+1)
        !   rho(i)=rho(i-n+1)
        !   m(i)=m(i-n+1)
        !   cs(i)=cs(i-n+1)
        !   h(i)=h(i-n+1)
        !   p(i)=p(i-n+1)
        !   a(i)=a(i-n+1)
        ! enddo
        !
        !
        ! ! the ghost particles which exist before the first particles e.g. -1, -2...
        ! do i=1,ng/2
        !   ! x(i+n+ng)=-dx*i
        !   x(i+n+ng/2)=x(n-i)-l
        !   ! x(i+n+ng)=x(1)-dx2*i
        !   v(i+n+ng/2)=v(n-i)
        !   rho(i+n+ng/2)=rho(n-i)
        !   m(i+n+ng/2)=m(n-i)
        !   cs(i+n+ng/2)=cs(n-i)
        !   h(i+n+ng/2)=h(n-i)
        !   p(i+n+ng/2)=p(n-i)
        !   a(i+n+ng/2)=a(n-i)
        ! enddo
        ! print*,x(n:n+ng)

        !
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
          endif
        enddo

        ! do i=n,ng
        !   v(i)=0.
        !   h(i)=h(i)
        !   m(i)=m(i)
        !   p(i)=p(i)
        !   rho(i)=0.1
        !   a(i)=a(i)
        !   cs(i)=cs(i)
        ! enddo

      end subroutine



end module