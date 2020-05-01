module inout2
  implicit none


contains
  subroutine output(x,v,a,t,l,u,l-dl,ifile)
    ! define the input parameters for printing
  integer,intent(in) :: ifile
  real,intent(in) :: x,v,a,t,l,u,l-dl)

  ! setup parameters for file writing
  integer :: i,iunit
  character(len=100) :: filename

  !write name of output files
  write(filename,"(a,i5.5)") 'snap_',ifile
  !open output files
  open(newunit=iunit,file=filename,status='replace')

  print "(a,f8.3)", ' writing '//trim(filename)// ' t =',t
  ! write in column headers and values
  write(iunit,*) '# x y v_x v_y a_x a_y t l e dl'
  write(iunit,*) t
  do i=1,n+ng
    !write in variables to file
    write(iunit,*) x(1),x(2),v(1),v(2),a(1),a(2),t,l,u,l-dl
  enddo
  close(iunit)



  end subroutine
end module
