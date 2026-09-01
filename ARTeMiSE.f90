program ARTeMiSE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Written by: Francisco Rendon
!On August 27, 2015
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use create_mesh
use intersections
use read_data
implicit none
!!!!!!!!!!!!!!!!!!!!!!!!
!Variables Deffinition
!!!!!!!!!!!!!!!!!!!!!!!!
real(8) :: N
integer :: i,j,k
integer :: id_i,id_j
real(8) :: X,Y,RHO
integer :: alloc_status
integer :: row_id_max, col_id
integer :: block_size
integer :: ndata, nvar_data
real(8) :: Delta_phi
real(8) :: rmin
real(8) :: Delta_r
real(8) :: Z_ray
real(8) :: Delta_Z_ray
real(8) :: Z_max
integer :: N_rays
integer :: kmax
integer :: i_kmax, j_kmax
real(8) :: vertex_aux(4,2)
real(8), dimension(:,:), allocatable :: data_mat !the cartesian meshh
!real(8), dimension(:), allocatable :: vec !the cartesian meshh
real(8), dimension(10):: vec_data ! This vector contains: ndata,nvar_data
                                  ! block_size,Delta_phi,rmin,Delta_r, etc.
real(8) :: Ivec1(2), Ivec2(2)
character(50):: my_input_file
character(50):: my_output_file
real(8) :: X_p,Y_p
real(8) :: parameter, my_pi = 3.14159265358979
integer :: i_id, j_id


call hello()
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Reading FARGO file
my_input_file = 'data.dat'
!my_input_file = 'file.data'
my_input_file = trim(my_input_file)
call read_data_from_FARGO(my_input_file, ndata, nvar_data, block_size)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Input parameters
!print*, ndata, nvar_data, block_size
Delta_phi = 9*my_pi/180.0 !from simulation
rmin = 1.0 !from simulation
Delta_r = 0.25 !size of the cell, from simulation
N_rays = 40
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
my_output_file = trim('intersections.dat')

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Asigning the initial parameters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
vec_data(1) = ndata
vec_data(2) = nvar_data
vec_data(3) = block_size
vec_data(4) = Delta_phi
vec_data(5) = rmin
vec_data(6) = Delta_r
vec_data(7) = real(N_rays)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Allocating memory to create data_mat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
allocate(data_mat(ndata,nvar_data), stat= alloc_status )
   if ( alloc_status /= 0 ) stop "Memory allocation error for: data_mat"

do i=1, ndata
   call get_points(FARGOdataMat(i,1),FARGOdataMat(i,2),X_p,Y_p)
   data_mat(i,1) = X_p
   data_mat(i,2) = Y_p
   print*, X_p, Y_p, i
enddo


call get_mesh(data_mat,vec_data)

kmax = ndata - block_size + 1!id of the first cell of the last column 
Z_max = mesh(kmax)%vertex(1,2)
Delta_Z_ray = Z_max/N_rays

open(unit=10,file=my_output_file)
Z_ray = Z_max !- Delta_Z_ray
do k=1,ndata 
   j = 1
   do 
       if ((Z_ray.ge.mesh(k)%z_min).and.(Z_ray.le.mesh(k)%z_max)) then
          do i=1,4
             vertex_aux(i,1) =  mesh(k)%vertex(i,1)
             vertex_aux(i,2) =  mesh(k)%vertex(i,2)
          enddo
          call get_intersections(Z_ray,vertex_aux,Ivec1,Ivec2)
          call get_mesh_id(k,ndata,block_size,i_id,j_id)

          mesh(k)%IMat(j,1) = Ivec1(1)
          mesh(k)%IMat(j,2) = Ivec1(2)
          mesh(k)%IMat(j,3) = Ivec2(1)
          mesh(k)%IMat(j,4) = Ivec2(2)
          write(10,*) mesh(k)%IMat(j,1),mesh(k)%IMat(j,2),mesh(k)%IMat(j,3),mesh(k)%IMat(j,4), k, Z_ray
          j = j + 1
       endif
       !print*, 'k', k, 'j', j, Z_ray, Delta_Z_ray

       if (Z_ray<mesh(k)%z_min)  exit
       Z_ray = Z_ray - Delta_Z_ray
   enddo
   Z_ray = Z_max - Delta_Z_ray
   Z_ray = Z_max! - Delta_Z_ray
enddo
close(10)
deallocate(data_mat)
deallocate(mesh)

end program ARTeMiSE

