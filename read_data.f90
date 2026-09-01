module read_data
implicit none
!character(30) :: my_input_file
!character(30) :: my_output_file
!integer :: ndata, nvar
!integer :: i,j,ios
real(8), dimension(:,:), allocatable :: FARGOdataMat

!my_input_file = 'data.dat'
!my_input_file = 'file.data'
!my_input_file = trim(my_input_file)

!call read_data(my_input_file,ndata, nvar)

!open (unit=30, file=trim(my_output_file), status='old',iostat=ios)
!do i=1, ndata
!call get_points(r,phi,X,Y)
!enddo
!close(30)

contains


subroutine read_data_from_FARGO(my_input_file, no_data, no_variables, block_size)
character(30) :: my_input_file
integer:: no_data, no_variables, block_size
character(1000) :: string, string_aux, string_empty
integer:: ios,filetty,filelog
integer:: rownum, error, nchar, nchar_aux
integer:: alloc_status
integer:: i,j,serror
integer:: separator, separator_empty
integer:: columnnum
!real(8), dimension(:,:), allocatable :: FARGOdataMat
real(8) :: par1, par2
integer :: bz, nblock
integer :: k, knew, kold
real(8), dimension(:,:), allocatable :: FARGOdataMatAux
character(30) :: my_aux_file
integer :: header

my_aux_file = 'aux_data.dat'
my_aux_file = trim(my_aux_file)

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
         string = trim(adjustl(string))
         nchar = len( trim(adjustl(string)) )
         string_aux = string
         string_empty = string
            !print*, 'string',string, string_aux
         separator = scan(string_aux,' ')
         nchar_aux = len( trim(adjustl(string_aux)) )
         do 
            if (separator > 0) then
            string_aux = trim(adjustl(string_aux(separator:nchar_aux)))
            end if
            j = j +1
            columnnum = j
            separator = scan(string_aux,' ')
            !print*, 'j',  j
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
          separator_empty= scan(string_empty,'    ')
          if (separator_empty==1) then
            bz = rownum 
!            print*, 'bz', bz
          endif

      enddo
      rewind(10)
      no_data  = rownum
      no_variables = columnnum
      block_size = no_data - bz
   
   close(10)
   !print*, 'ndata', no_data, 'nvar', no_variables
      if (columnnum == 1) then
          print*, 'Insuficent data, I can start'
          STOP
      endif


   allocate(FARGOdataMat(no_data,no_variables), stat= alloc_status )
       if ( alloc_status /= 0 ) stop "Memory allocation error for: FARGOdataMat"
   allocate(FARGOdataMatAux(no_data+block_size-1,no_variables+1), stat= alloc_status )
       if ( alloc_status /= 0 ) stop "Memory allocation error for: FARGOdataMatAux"

 open (unit=40, file=trim(my_aux_file))
     open (unit=10, file=trim(my_input_file), status='old', iostat=ios)
        do i=1,no_data+block_size!+1
         read(10,'(A180)' , iostat = error) string
         string = adjustl(string)
         nchar = len( trim(adjustl(string)) )
          if (error .ne. 0) exit
          if ( ( string(1:1) /= '!' ).and.  &
                   ( string(1:1) /= '#' ).and.  &
                   ( string(1:1) /= ' ' ).and.  &
                   ( nchar > 0 ) &
           ) then
           
           !print*, i, nchar, trim(string)
           read(string,* , iostat = serror) par1, par2
           write(40,*) par1 , par2
           FARGOdataMatAux(i,1) = par1
           FARGOdataMatAux(i,2) = par2
           endif
           
          if ( ( string(1:1) == '!' ).or.  &
                   ( string(1:1) == '#' ).or.  &
                   ( string(1:1) == ' ' )) &
          then
              header = 1 ! Means the input FARGO file has a header
              FARGOdataMatAux(i,no_variables + 1) = header
              
          endif   
        enddo
     close(10)
close(40)

nblock = no_data/block_size

kold = 1
do i=1,nblock
   knew = i*block_size
   do k = kold, knew
      if (FARGOdataMatAux(1,no_variables + 1 )== 1) then
         j = k + (i - 1) + 1
      else
         j = k + (i - 1) 
      endif
      FARGOdataMat(k,1) = FARGOdataMatAux(j,1) 
      FARGOdataMat(k,2) = FARGOdataMatAux(j,2) 
!      print*, 'nblock',nblock, 'kold',kold, 'knew',knew, 'k',k,'j', j
   enddo
   kold = knew + 1
enddo
   
do i=1, no_data
  ! print*, 'sal',FARGOdataMat(i,1), FARGOdataMat(i,2), i
enddo
end subroutine read_data_from_FARGO

subroutine get_points(r,phi,X,Y)
real(8):: r, phi
real(8):: X, Y
real(8):: my_pi
!
my_pi = 3.14159265358979 
phi = phi*my_pi/180.0
X = r*cos(phi)
Y = r*sin(phi)
end subroutine get_points

end module read_data
