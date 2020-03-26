module inout
  implicit none


contains
  subroutine output(n,x,v,h,nx,rho,m,p,cs,a,t,ifile,ek)
    ! define the input parameters for printing
  integer,intent(in) :: n,nx,ifile
  real,intent(in) :: x(nx),v(nx),h(nx),rho(nx),m(nx),p(nx),cs(nx),a(nx),t,ek

  ! setup parameters for file writing
  integer :: i,iunit
  character(len=100) :: filename

  write(filename,"(a,i5.5)") 'snap_',ifile
  open(newunit=iunit,file=filename,status='replace')

  print "(a,f8.3)", ' writing '//trim(filename)// ' t =',t
  ! write in column headers and values
  write(iunit,*) '# x v h rho m p cs a t ek'
  write(iunit,*) t
  do i=1,nx

    write(iunit,*) x(i),v(i),h(i),rho(i),m(i),p(i),cs(i),a(i), t, ek
  enddo
  close(iunit)



  end subroutine
end module
