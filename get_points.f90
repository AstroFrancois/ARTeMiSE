program create_points
implicit none
character(30) :: my_input_file
character(30) :: my_output_file
integer :: ndata, nvar
integer :: i,j,ios
real(8), dimension(:,:), allocatable :: dataMat

my_input_file = 'data.dat'
!my_input_file = 'file.data'
my_input_file = trim(my_input_file)

call read_data(my_input_file,ndata, nvar)

open (unit=30, file=trim(my_output_file), status='old',iostat=ios)
do i=1, ndata
!call get_points(r,phi,X,Y)
enddo
close(30)

end program create_points


subroutine read_data(my_input_file, ndata, nvar)
character(30) :: my_input_file
integer:: ndata, nvar
character(1000) :: string, string_aux
integer:: ios,filetty,filelog
integer:: rownum, error, nchar, nchar_aux
integer:: alloc_status
integer:: i,j,serror
integer:: separator
integer:: columnnum
real(8), dimension(:,:), allocatable :: dataMat
real(8) :: par1, par2



   open (unit=10, file=trim(my_input_file), status='old', iostat=ios)
    if( ios /= 0 ) then
       write(filetty,'(a,a)') 'error opening file:', my_input_file, '"'
       write(filelog,'(a,a)') 'error opening file:', my_input_file, '"'
       stop
    endif
      rownum=0
      j = 0
      do
         read(10,'(A180)' , iostat = error) string
         string = adjustl(string)
         nchar = len( trim(adjustl(string)) )
         string_aux = string
         separator = scan(string_aux,' ')
         nchar_aux = len( trim(adjustl(string_aux)) )
         do 
            if (separator > 0) then
            string_aux = trim(adjustl(string_aux(separator:nchar_aux)))
            end if
            j = j +1
            columnnum = j
            separator = scan(string_aux,' ')
            !print*, 'separator out', separator, string_aux
            nchar_aux = len( trim(adjustl(string_aux)) )
            if (separator == 1) exit
         enddo
         j = 0
         !print*, 'columnnum', columnnum

          if (error .ne. 0) exit
          if ( ( string(1:1) /= '!' ).and.  &
                   ( string(1:1) /= '#' ).and.  &
                   ( nchar > 0 ) &
           ) then
            rownum = rownum + 1
          endif

      enddo
      rewind(10)
      ndata  = rownum
   
   close(10)
   !print*, 'ndata', ndata


   allocate(dataMat(ndata,2), stat= alloc_status )
       if ( alloc_status /= 0 ) stop "Memory allocation error for: dataMat"

   !open(unit=20,file=trim(my_input_file))
     open (unit=10, file=trim(my_input_file), status='old', iostat=ios)
        do i=1,ndata+1
         read(10,'(A180)' , iostat = error) string
         string = adjustl(string)
         nchar = len( trim(adjustl(string)) )
          if (error .ne. 0) exit
          if ( ( string(1:1) /= '!' ).and.  &
                   ( string(1:1) /= '#' ).and.  &
                   ( nchar > 0 ) &
           ) then


           read(string,* , iostat = serror) par1, par2
           dataMat(i,1) = par1
           dataMat(i,2) = par2
           !print*, dataMat(i,1), dataMat(i,2)
           !write(20,*) xvec(i), yvec(i), densvec(i)
           endif
        enddo
     close(10)
   !close(20)
end subroutine read_data

subroutine get_points(r,phi,X,Y)
real(8):: r, phi
real(8):: X, Y
!
X = r*cos(phi)
Y = r*sin(phi)
end subroutine get_points

