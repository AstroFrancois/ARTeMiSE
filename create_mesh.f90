module create_mesh
use intersections
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Written by: Francisco Rendon
!On August 25, 2015
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
implicit none
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!GLOBAL VARIABLES to use in the main program (This variables 
!need not to be defined there in) and inside subroutines 
!(here).
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!type declaration
type simulation_data
   integer :: id(2)
   real(8) :: X,Y,RHO
   real(8) :: VAR_PHI
   real(8) :: PHI
   real(8) :: vertex(4,2) 
   real(8) :: PHI_min
   real(8) :: PHI_max
   real(8) :: R_min
   real(8) :: R_max
   real(8) :: z_min
   real(8) :: z_max
   real(8), dimension(:,:), allocatable :: IMat !This matrix saves the
                                                   !interesction of the 
                                                   !cell with the rays 
   integer :: rays_per_cell
end type simulation_data
!!!!!!!!!!!!!!!!!!!!!!!!!
   
!declaring array of mesh
type(simulation_data), dimension(:), allocatable :: mesh 
character(30) :: welcome




contains
subroutine get_mesh(data_mat,vec_data)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!The input is the matrix containing all data
!x1,y1,z1,den1,...,vz1
!x2,y2,z2,den2,...,vz2
!. ,. , ., .  ,..., .
!xn,yn,zn,denn,...,vzn
real(8), dimension(:,:), allocatable :: data_mat !the cartesian meshh
real(8), dimension(10):: vec_data ! This vector contains: ndata,nvar_data
                                  ! block_size,Delta_phi,rmin,Delta_r, etc.
integer :: block_size !this is the size of the block
real(8) :: Delta_phi
integer :: row_id_max !this is the size of all rows from the ordered file
integer :: row_id
integer :: col_id
integer :: i,j,k
integer :: ndata, nvar
integer :: N_rays
integer :: i_id, j_id
integer :: alloc_status
real(8) :: Delta_r
real(8) :: Delta_Z
real(8) :: cell_Zsize
real(8) :: rmin !from simulation
integer :: N_rays_per_cell
integer :: kmax
real(8) :: Z_max
real(8) :: Delta_Z_ray
real(8), dimension(4) :: Zcomp_vertex
character(30) :: my_output_file

!!!!!!!!!!!!!!!!!!!!!!!
!Asigning parameter data
!!!!!!!!!!!!!!!!!!!!!!

my_output_file = trim('vertex.dat')

ndata = vec_data(1)
nvar = vec_data(2)
block_size = vec_data(3)
Delta_phi = vec_data(4)
rmin = vec_data(5)
Delta_r = vec_data(6)
N_rays = int(vec_data(7))
!print*, 'NRays', N_rays
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Allocating memory for the mesh (type) variable
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
allocate(mesh(ndata), stat= alloc_status )
   if ( alloc_status /= 0 ) stop "Memory allocation error for: mesh in subroutine get_mesh"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
k = 1
open(unit=20,file=my_output_file)
do k=1, ndata
   call get_mesh_id(k,ndata,block_size,i_id,j_id)
      mesh(k)%id(1) = i_id
      mesh(k)%id(2) = j_id
      mesh(k)%X = data_mat(k,1)
      mesh(k)%Y = data_mat(k,2)
      !print*, 'X=',mesh(k)%X, 'Y=',mesh(k)%Y
      mesh(k)%RHO = data_mat(k,3)
      mesh(k)%VAR_PHI = get_VAR_PHI(mesh(k)%X,mesh(k)%Y)
      mesh(k)%PHI = get_REAL_PHI(mesh(k)%VAR_PHI,Delta_phi)
      mesh(k)%PHI_min = mesh(k)%PHI - 0.5*Delta_phi
      mesh(k)%PHI_max = mesh(k)%PHI + 0.5*Delta_phi
      mesh(k)%R_min = rmin + (j_id - 1)*Delta_r
      mesh(k)%R_max = rmin + j_id*Delta_r
      mesh(k)%vertex(1,1) = get_Xcomp(mesh(k)%R_max,mesh(k)%PHI_max)
      mesh(k)%vertex(1,2) = get_Ycomp(mesh(k)%R_max,mesh(k)%PHI_max)
      mesh(k)%vertex(2,1) = get_Xcomp(mesh(k)%R_min,mesh(k)%PHI_max)
      mesh(k)%vertex(2,2) = get_Ycomp(mesh(k)%R_min,mesh(k)%PHI_max)
      mesh(k)%vertex(3,1) = get_Xcomp(mesh(k)%R_min,mesh(k)%PHI_min)
      mesh(k)%vertex(3,2) = get_Ycomp(mesh(k)%R_min,mesh(k)%PHI_min)
      mesh(k)%vertex(4,1) = get_Xcomp(mesh(k)%R_max,mesh(k)%PHI_min)
      mesh(k)%vertex(4,2) = get_Ycomp(mesh(k)%R_max,mesh(k)%PHI_min)
      write(20,*) mesh(k)%vertex(1,1), mesh(k)%vertex(1,2), k, mesh(k)%id(1), mesh(k)%id(2)
      write(20,*) mesh(k)%vertex(2,1), mesh(k)%vertex(2,2), k, mesh(k)%id(1), mesh(k)%id(2) 
      write(20,*) mesh(k)%vertex(3,1), mesh(k)%vertex(3,2), k, mesh(k)%id(1), mesh(k)%id(2)
      write(20,*) mesh(k)%vertex(4,1), mesh(k)%vertex(4,2), k, mesh(k)%id(1), mesh(k)%id(2)
      write(20,*) k
      do j=1,4
      Zcomp_vertex(j) = mesh(k)%vertex(j,2)
      enddo
      mesh(k)%z_min = minval(Zcomp_vertex)
      mesh(k)%z_max = maxval(Zcomp_vertex)
      !do i=1,4 
      !      print*,  mesh(k)%vertex(i,1), mesh(k)%vertex(i,2)  
      !enddo
      !cell_Zsize = get_cell_Zsize(mesh(k)%z_min,mesh(k)%z_max)
      !N_rays_per_cell = get_number_of_rays(Delta_z,cell_Zsize)
enddo
close(20)
!deallocate(mesh)

kmax = ndata - block_size + 1!id of the first cell of the last column 
!call get_mesh_id(kmax,ndata,block_size,i_kmax,j_kmax)
Z_max = mesh(kmax)%vertex(1,2)
!print*, 'Z_max',Z_max
Delta_Z_ray = Z_max/N_rays

do k=1, ndata
     cell_Zsize = get_cell_Zsize(mesh(k)%z_min,mesh(k)%z_max)
     N_rays_per_cell = get_number_of_rays(Delta_Z_ray,cell_Zsize)
     !print*, 'k',k, cell_Zsize, N_rays_per_cell
     mesh(k)%rays_per_cell = N_rays_per_cell
     allocate(mesh(k)%IMat(N_rays_per_cell,4), stat= alloc_status )
        if ( alloc_status /= 0 ) stop "Memory allocation error for: Interesction Matrix in subroutine get_mesh"
     
enddo
return
end subroutine get_mesh

function get_VAR_PHI(x,y) result(VAR_PHI)
real(8) :: x, y
real(8) :: VAR_PHI

VAR_PHI = atan(y/x)
end function get_VAR_PHI

function get_Xcomp(r,theta) result(X)
real(8) :: r, theta
real(8) :: X
X = r*cos(theta)
end function get_Xcomp

function get_Ycomp(r,theta) result(Y)
real(8) :: r, theta
real(8) :: Y
Y = r*sin(theta)
end function get_Ycomp

function get_REAL_PHI(VAL_PHI,Delta_phi) result(PHI)
real(8) :: VAL_PHI
real(8) :: Delta_phi
real(8) :: a
integer :: a_integer
integer :: a_hat
real(8) :: PHI

a = 2.0*VAL_PHI/Delta_phi
!print*, 'Dphi',Delta_phi, 'a',a
a_integer = nint(a)
a_hat = real(a_integer)
PHI = 0.5*a_hat*Delta_phi

end function get_REAL_PHI

subroutine get_data_mat_id(i,j,k)
!!!!!!!!!!!!!!!!!!!!
!input parameter
integer :: i,j
!!!!!!!!!!!!!!!!!!!!
!output parameters
integer :: k
!!!!!!!!!!!!!!!!!!!!
integer :: block_size

k = (j -1)*block_size + i

end subroutine get_data_mat_id

subroutine get_mesh_id(k,ndata,block_size,i,j)
!!!!!!!!!!!!!!!!!!!!
!input parameter
integer :: k
integer :: ndata
integer :: block_size
!!!!!!!!!!!!!!!!!!!!
!output parameters
integer :: i,j
!!!!!!!!!!!!!!!!!!!!
integer :: i_max, j_max, k_max
integer :: j_integer
real(8) :: j_real
!block_size = 5
k_max = ndata
i_max = block_size
j_max = k_max/i_max

j_real = real(k)/real(block_size)
j_integer = k/block_size
if (j_real>real(j_integer) ) then
   j = j_integer + 1
else
   j = j_integer
endif
!j = k/(i_max + 1) + 1
i = k -(j - 1)*block_size
!print*, 'imax',i_max, 'jmax',j_max, 'k',k, 'kmax', k_max, 'i',i, 'j',j
!return

!print*, 'k=',k,'i=',i, 'j=',j

end subroutine get_mesh_id


subroutine get_cell_data(i,j,cell_info)
!!!!!!!!!!!!!!!!!!!!!!!!
!Variables Deffinition
!!!!!!!!!!!!!!!!!!!!!!!!
integer :: i,j
integer :: iter_n
real(8), dimension(3) :: cell_info 
real(8) :: X,Y,RHO

X = 1.0
Y = 2.0
RHO = 3.0

cell_info(1) = X
cell_info(2) = Y
cell_info(3) = RHO

!print*, X, Y, RHO
end subroutine get_cell_data

subroutine hello()
welcome = 'ARTeMiSE vous accueille!'
welcome = trim(welcome)
print*, welcome
end subroutine 




end module create_mesh
