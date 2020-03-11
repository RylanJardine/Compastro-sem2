module inout
  implicit none


contains
  subroutine output(n,x,v,h,nx,rho,m,p,cs,a)
    ! define the input parameters for printing
  real,intent(in) :: x(nx),v(nx),h(nx),rho(nx),m(nx),p(nx),cs(nx),a(nx)
  integer,intent(in) :: n,nx
  ! setup parameters for file writing
  real :: time
  integer :: i,iunit
  character(len=100) :: filename

  write(filename,"(a,i5.5)") 'snap_'
  open(newunit=iunit,file=filename,status='replace')
  ! write in column headers and values
  write(iunit,*) '# x v h rho m p cs a'
  write(iunit,*) time
  do i=1,n

    write(iunit,*) x(i),v(i),h(i),rho(i),m(i),p(i),cs(i),a(i)
  enddo
  close(iunit)


  end subroutine
end module
