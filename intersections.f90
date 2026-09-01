module intersections
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Written by: Francisco Rendon
!On August 25, 2015
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
implicit none


contains 

subroutine get_intersections(Z_ray,vertex,Ivec1,Ivec2)
real(8) :: Z_ray
real(8) :: vertex(4,2)
real(8) :: Ivec1(2)
real(8) :: Ivec2(2)
real(8) :: Pvec1(2)
real(8) :: Pvec2(2)
real(8) :: phi
real(8) :: rC1, rC2
real(8) :: slope
real(8) :: Z_V1, Z_V2, Z_V3, Z_V4

Z_V1 = vertex(1,2)
Z_V2 = vertex(2,2)
Z_V3 = vertex(3,2)
Z_V4 = vertex(4,2)

!print*, 'Z_ray in subroutine', Z_ray 
if (Z_V2.ge.Z_V4) then

   if ((Z_ray.ge.Z_V2).and.(Z_ray.le.Z_V1)) then !intersection with L1 and C2
       !print*, 'intersection with line 1 and circunference 2' 
       !intersection with line L1 
       Pvec1(1) = vertex(1,1)
       Pvec1(2) = vertex(1,2)
       Pvec2(1) = vertex(2,1)
       Pvec2(2) = vertex(2,2)
       slope = get_slope(Pvec1,Pvec2)
       Ivec1 = get_intersection_wLINE(Z_ray,slope)
       !intersection with circunference C2
       rC2 = sqrt(vertex(1,1)**2.0 + vertex(1,2)**2.0)
       Ivec2 = get_intersection_wCIRCUNFERENCE(Z_ray,rC2)
   endif

   if ((Z_ray.ge.Z_V4).and.(Z_ray<Z_V2)) then !intersection with C1 and C2
       !print*, 'intersection with circunfernce 1 and circunference 2' 
       !intersection with circunference C1
       rC1 = sqrt(vertex(2,1)**2.0 + vertex(2,2)**2.0)
       Ivec1 = get_intersection_wCIRCUNFERENCE(Z_ray,rC1)
       !intersection with circunference C2
       rC2 = sqrt(vertex(1,1)**2.0 + vertex(1,2)**2.0)
       Ivec2 = get_intersection_wCIRCUNFERENCE(Z_ray,rC2)
   endif

   if ((Z_ray.ge.Z_V3).and.(Z_ray<Z_V4)) then !intersection with C1 and L2
       !intersection with circunference C1
       rC1 = sqrt(vertex(2,1)**2.0 + vertex(2,2)**2.0) 
       Ivec1 = get_intersection_wCIRCUNFERENCE(Z_ray,rC1)
       !print*, 'intersection with circunference 1 and line 2' 
       !intersection with line L2
       Pvec1(1) = vertex(3,1)
       Pvec1(2) = vertex(3,2)
       Pvec2(1) = vertex(4,1)
       Pvec2(2) = vertex(4,2)
       slope = get_slope(Pvec1,Pvec2)
       Ivec2 = get_intersection_wLINE(Z_ray,slope)
   endif

elseif(Z_V4>Z_V2) then
   if ((Z_ray.ge.Z_V4).and.(Z_ray.le.Z_V1)) then !intersection with L1 and C2
       !print*, 'intersection with line 1 and circunference 2' 
       !intersection with line L1 
       Pvec1(1) = vertex(1,1)
       Pvec1(2) = vertex(1,2)
       Pvec2(1) = vertex(2,1)
       Pvec2(2) = vertex(2,2)
       slope = get_slope(Pvec1,Pvec2)
       Ivec1 = get_intersection_wLINE(Z_ray,slope)
       !intersection with circunference C2
       rC2 = sqrt(vertex(1,1)**2.0 + vertex(1,2)**2.0)
       Ivec2 = get_intersection_wCIRCUNFERENCE(Z_ray,rC2)
   endif

   if ((Z_ray.ge.Z_V2).and.(Z_ray<Z_V4)) then !intersection with L1 and L2
       !print*, 'intersection with circunfernce 1 and circunference 2' 
       !intersection with circunference L1
       Pvec1(1) = vertex(1,1)
       Pvec1(2) = vertex(1,2)
       Pvec2(1) = vertex(2,1)
       Pvec2(2) = vertex(2,2)
       slope = get_slope(Pvec1,Pvec2)
       Ivec1 = get_intersection_wLINE(Z_ray,slope)
       !intersection with circunference L2
       Pvec1(1) = vertex(3,1)
       Pvec1(2) = vertex(3,2)
       Pvec2(1) = vertex(4,1)
       Pvec2(2) = vertex(4,2)
       slope = get_slope(Pvec1,Pvec2)
       Ivec2 = get_intersection_wLINE(Z_ray,slope)
   endif

   if ((Z_ray.ge.Z_V3).and.(Z_ray<Z_V2)) then !intersection with C1 and L2
       !intersection with circunference C1
       rC1 = sqrt(vertex(2,1)**2.0 + vertex(2,2)**2.0) 
       Ivec1 = get_intersection_wCIRCUNFERENCE(Z_ray,rC1)
       !print*, 'intersection with circunference 1 and line 2' 
       !intersection with line L2
       Pvec1(1) = vertex(3,1)
       Pvec1(2) = vertex(3,2)
       Pvec2(1) = vertex(4,1)
       Pvec2(2) = vertex(4,2)
       slope = get_slope(Pvec1,Pvec2)
       Ivec2 = get_intersection_wLINE(Z_ray,slope)
   endif
endif

!print*, Z_ray, Ivec1(1), Ivec1(2), Ivec2(1), Ivec2(2)
!print*, 'in subroutine', Ivec1(1), Ivec1(2), Ivec2(1), Ivec2(2)
!print*, vertex(1,1), vertex(1,2)
!print*, vertex(2,1), vertex(2,2)
!print*, vertex(3,1), vertex(3,2)
!print*, vertex(4,1), vertex(4,2)

return
end subroutine get_intersections


function get_cell_Zsize(zmin,zmax) result(Cell_Zsize)
real(8) :: zmin
real(8) :: zmax
real(8) :: Cell_Zsize

Cell_Zsize = zmax - zmin
end function get_cell_Zsize

function get_number_of_rays(Delta_z,Cell_Zsize) result(N_rays)
real(8) :: Delta_z
real(8) :: Cell_Zsize
real(8) :: N_rays_real
integer :: N_rays

N_rays_real = Cell_Zsize/Delta_z
N_rays = nint(N_rays_real)
if (N_rays_real>real(N_rays)) then
N_rays = N_rays + 1
endif
!print*, 'Cell_Z_size', Cell_Zsize, 'Dz', Delta_z, 'NReal',N_rays_real, 'N',N_rays

end function get_number_of_rays

function get_slope_by_angle(phi) result(slope)
real(8)  :: phi
real(8)  :: slope

slope = atan(phi)
end function get_slope_by_angle

function get_slope(P1,P2) result(slope)
real(8)  :: P1(2), P2(2)
real(8)  :: slope

slope = P2(2) - P1(2)
slope = slope/(P2(1) - P1(1))
end function get_slope

function get_intersection_wLINE(Y,slope) result(Ivec)
real(8) :: Y
real(8) :: slope
real(8) :: Ivec(2)

Ivec(1) = Y/slope
Ivec(2) = Y
end function get_intersection_wLINE

function get_intersection_wCIRCUNFERENCE(Y,radius) result(Ivec)
real(8) :: Y
real(8) :: radius
real(8) :: Ivec(2)

Ivec(1) = sqrt(radius**2.0 - Y**2.0)
Ivec(2) = Y
end function get_intersection_wCIRCUNFERENCE

end module intersections
