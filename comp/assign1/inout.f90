module inout
  implicit none


contains
  subroutine output(n,x,v,h,nx,rho,m,p,cs,a,t,ifile,ng,u,du)
    ! define the input parameters for printing
  integer,intent(in) :: n,nx,ifile,ng
  real,intent(in) :: x(nx),v(nx),h(nx),rho(nx),m(nx),p(nx),cs(nx),a(nx),t,u(nx),du(nx)

  ! setup parameters for file writing
  integer :: i,iunit
  character(len=100) :: filename


  write(filename,"(a,i5.5)") 'snap_',ifile
  ! print*,filename, iunit, 'bob'
  open(newunit=iunit,file=filename,status='replace')

  print "(a,f8.3)", ' writing '//trim(filename)// ' t =',t
  ! write in column headers and values
  write(iunit,*) '# x v_x h rho m p cs a u du'
  write(iunit,*) t
  do i=1,n+ng

    write(iunit,*) x(i),v(i),h(i),rho(i),m(i),p(i),cs(i),a(i),u(i),du(i)
  enddo
  close(iunit)



  end subroutine
end module
