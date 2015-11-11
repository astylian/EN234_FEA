!     Subroutines for basic 3D linear elastic elements 
!
!==========================SUBROUTINE el_linelast_3dbasic ==============================
subroutine el_linelast_3dbasic(lmn, element_identifier, n_nodes, node_property_list, &           ! Input variables
    n_properties, element_properties, element_coords, length_coord_array, &                      ! Input variables
    dof_increment, dof_total, length_dof_array, &                                                ! Input variables
    n_state_variables, initial_state_variables, &                                                ! Input variables
    updated_state_variables,element_stiffness,element_residual, fail)      ! Output variables                          ! Output variables
    use Types
    use ParamIO
    !  use Globals, only: TIME,DTIME  For a time dependent problem uncomment this line to access the time increment and total time
    use Mesh, only : node
    use Element_Utilities, only : N => shape_functions_3D
    use Element_Utilities, only : dNdxi => shape_function_derivatives_3D
    use Element_Utilities, only:  dNdx => shape_function_spatial_derivatives_3D
    use Element_Utilities, only : xi => integrationpoints_3D, w => integrationweights_3D
    use Element_Utilities, only : dxdxi => jacobian_3D
    use Element_Utilities, only : initialize_integration_points
    use Element_Utilities, only : calculate_shapefunctions
    use Element_Utilities, only : invert_small
    use Element_Utilities, only : dNbardx => vol_avg_shape_function_derivatives_3D
    implicit none

    integer, intent( in )         :: lmn                                                    ! Element number
    integer, intent( in )         :: element_identifier                                     ! Flag identifying element type (specified in .in file)
    integer, intent( in )         :: n_nodes                                                ! # nodes on the element
    integer, intent( in )         :: n_properties                                           ! # properties for the element
    integer, intent( in )         :: length_coord_array                                     ! Total # coords
    integer, intent( in )         :: length_dof_array                                       ! Total # DOF
    integer, intent( in )         :: n_state_variables                                      ! # state variables for the element

    type (node), intent( in )     :: node_property_list(n_nodes)                  ! Data structure describing storage for nodal variables - see below
    !  type node
    !      sequence
    !      integer :: flag                          ! Integer identifier
    !      integer :: coord_index                   ! Index of first coordinate in coordinate array
    !      integer :: n_coords                      ! Total no. coordinates for the node
    !      integer :: dof_index                     ! Index of first DOF in dof array
    !      integer :: n_dof                         ! Total no. of DOF for node
    !   end type node
    !   Access these using node_property_list(k)%n_coords eg to find the number of coords for the kth node on the element

    real( prec ), intent( in )    :: element_coords(length_coord_array)                     ! Coordinates, stored as x1,(x2),(x3) for each node in turn
    real( prec ), intent( in )    :: dof_increment(length_dof_array)                        ! DOF increment, stored as du1,du2,du3,du4... for each node in turn
    real( prec ), intent( in )    :: dof_total(length_dof_array)                            ! accumulated DOF, same storage as for increment

    real( prec ), intent( in )    :: element_properties(n_properties)                       ! Element or material properties, stored in order listed in input file
    real( prec ), intent( in )    :: initial_state_variables(n_state_variables)             ! Element state variables.  Defined in this routine
  
    logical, intent( out )        :: fail                                                   ! Set to .true. to force a timestep cutback
    real( prec ), intent( inout ) :: updated_state_variables(n_state_variables)             ! State variables at end of time step
    real( prec ), intent( out )   :: element_stiffness(length_dof_array,length_dof_array)   ! Element stiffness (ROW,COLUMN)
    real( prec ), intent( out )   :: element_residual(length_dof_array)                     ! Element residual force (ROW)
          

    ! Local Variables
    integer      :: n_points,kint,i

    real (prec)  ::  strain(6), dstrain(6)             ! Strain vector contains [e11, e22, e33, 2e12, 2e13, 2e23]
    real (prec)  ::  stress(6)                         ! Stress vector contains [s11, s22, s33, s12, s13, s23]
    real (prec)  ::  D(6,6)                            ! stress = D*(strain+dstrain)  (NOTE FACTOR OF 2 in shear strain)
    real (prec)  ::  B(6,length_dof_array)             ! strain = B*(dof_total+dof_increment)
    real (prec)  ::  dxidx(3,3), determinant           ! Jacobian inverse and determinant
    real (prec)  ::  x(3,length_coord_array/3)         ! Re-shaped coordinate array x(i,a) is ith coord of ath node
    real (prec)  ::  E, xnu, D44, D11, D12             ! Material properties
    real (prec)  ::  Difference_matrix(6,3*n_nodes)    ! Empty matrix to be filled
    real (prec)  ::  integralargumentdNbardx
    real (prec)  ::  integralargumentelvol
    real (prec)  ::  elvol
    real (prec)  ::  F(3,3)
    real (prec)  ::  Iden(6)
    real (prec)  ::  u(length_dof_array)
    real (prec)  ::  umat(3,n_nodes)
    real (prec)  ::  integrand_r
    real (prec)  ::  r
    real (prec)  ::  Finv(3,3)
    real (prec)  ::  J
    real (prec)  ::  dNdy(n_nodes,3)
    real (prec)  ::  Bhyper(3,3)
    real (prec)  ::  Bhypervector(6)
    real (prec)  ::  Bhyperinverse(3,3)
    real (prec)  ::  useless
    real (prec)  ::  Bhypervectorinverse(6)
    real (prec)  ::  dsubmat(6,6)
    real (prec)  ::  G(6,9)
    real (prec)  ::  Bstar(9,6)
    real (prec)  ::  outerib(6,6)
    real (prec)  ::  outerii(6,6)
    real (prec)  ::  outerbb(6,6)
    real (prec)  ::  Pvec(length_dof_array)
    real (prec)  ::  Pmat(length_dof_array,length_dof_array)
    real (prec)  ::  Svec(length_dof_array)
    real (prec)  ::  Smat(length_dof_array,length_dof_array)
    real (prec)  ::  S(3,length_dof_array/3)
    real (prec)  ::  Sigma(length_dof_array,length_dof_array)



    !     Subroutine to compute element stiffness matrix and residual force vector for 3D linear elastic elements
    !     El props are:

    !     element_properties(1)         Young's modulus
    !     element_properties(2)         Poisson's ratio

    fail = .false.
    
    x = reshape(element_coords,(/3,length_coord_array/3/))

    if (n_nodes == 4) n_points = 1
    if (n_nodes == 10) n_points = 4
    if (n_nodes == 8) n_points = 8
    if (n_nodes == 20) n_points = 27

    call initialize_integration_points(n_points, n_nodes, xi, w)

    element_residual = 0.d0
    element_stiffness = 0.d0
	
    D = 0.d0
    E = element_properties(1)
    xnu = element_properties(2)
        d44 = 0.5D0*E/(1+xnu)
    d11 = (1.D0-xnu)*E/( (1+xnu)*(1-2.D0*xnu) )
    d12 = xnu*E/( (1+xnu)*(1-2.D0*xnu) )
    D(1:3,1:3) = d12
    D(1,1) = d11
    D(2,2) = d11
    D(3,3) = d11
    D(4,4) = d44
    D(5,5) = d44
    D(6,6) = d44
  
    elvol = 0.d0
    dNbardx = 0.d0

    !     --  Loop over integration points
    do kint = 1, n_points
        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        integralargumentelvol = determinant*w(kint)

        dNbardx = dNbardx + dNdx*w(kint)*determinant
        elvol = elvol + integralargumentelvol
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end do
    dNbardx = dNbardx/elvol


   !     --  Loop over integration points
    do kint = 1, n_points
        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)
!        B = 0.d0
!        B(1,1:3*n_nodes-2:3) = dNdx(1:n_nodes,1)
!        B(2,2:3*n_nodes-1:3) = dNdx(1:n_nodes,2)
!        B(3,3:3*n_nodes:3)   = dNdx(1:n_nodes,3)
!        B(4,1:3*n_nodes-2:3) = dNdx(1:n_nodes,2)
!        B(4,2:3*n_nodes-1:3) = dNdx(1:n_nodes,1)
!        B(5,1:3*n_nodes-2:3) = dNdx(1:n_nodes,3)
!        B(5,3:3*n_nodes:3)   = dNdx(1:n_nodes,1)
!        B(6,2:3*n_nodes-1:3) = dNdx(1:n_nodes,3)
!        B(6,3:3*n_nodes:3)   = dNdx(1:n_nodes,2)
!
!        !convert B to Bbar:
!
!
!        if (element_identifier == 1002) then
!            do i = 1,3
!                B(i,1:3*n_nodes-2:3) = B(i,1:3*n_nodes-2:3)&
!                +(dNbardx(1:n_nodes,1)-dNdx(1:n_nodes,1))/3.d0
!                B(i,2:3*n_nodes-1:3) = B(i,2:3*n_nodes-1:3)&
!                +(dNbardx(1:n_nodes,2)-dNdx(1:n_nodes,2))/3.d0
!                B(i,3:3*n_nodes:3) = B(i,3:3*n_nodes:3)&
!                +(dNbardx(1:n_nodes,3)-dNdx(1:n_nodes,3))/3.d0
!            end do
!        end if
!
!        strain = matmul(B,dof_total)
!        dstrain = matmul(B,dof_increment)
!        if (n_properties==2) then
!           stress = matmul(D,strain+dstrain)
!        else
!           call hypoelastic(strain,dstrain,element_properties,n_properties,stress,D)
!        endif

        u = dof_increment + dof_total

        Iden(4:6) = 0.d0
        Iden(1:3) = 1.d0
        umat = reshape((u), (/3,n_nodes/))
        F = 0.d0
        F = eye3_D + matmul(umat,dNdx)

        call invert_small(F,Finv,J)
        dNdy = matmul(dNdx,Finv)
        Bhyper = matmul(F,transpose(F))

        B = 0.d0
        if (element_identifier == 1005) then

        B(1,1:3*n_nodes-2:3) = dNdy(1:n_nodes,1)
        B(2,2:3*n_nodes-1:3) = dNdy(1:n_nodes,2)
        B(3,3:3*n_nodes:3)   = dNdy(1:n_nodes,3)
        B(4,1:3*n_nodes-2:3) = dNdy(1:n_nodes,2)
        B(4,2:3*n_nodes-1:3) = dNdy(1:n_nodes,1)
        B(5,1:3*n_nodes-2:3) = dNdy(1:n_nodes,3)
        B(5,3:3*n_nodes:3)   = dNdy(1:n_nodes,1)
        B(6,2:3*n_nodes-1:3) = dNdy(1:n_nodes,3)
        B(6,3:3*n_nodes:3)   = dNdy(1:n_nodes,2)

        else
        B(1,1:3*n_nodes-2:3) = dNdx(1:n_nodes,1)
        B(2,2:3*n_nodes-1:3) = dNdx(1:n_nodes,2)
        B(3,3:3*n_nodes:3)   = dNdx(1:n_nodes,3)
        B(4,1:3*n_nodes-2:3) = dNdx(1:n_nodes,2)
        B(4,2:3*n_nodes-1:3) = dNdx(1:n_nodes,1)
        B(5,1:3*n_nodes-2:3) = dNdx(1:n_nodes,3)
        B(5,3:3*n_nodes:3)   = dNdx(1:n_nodes,1)
        B(6,2:3*n_nodes-1:3) = dNdx(1:n_nodes,3)
        B(6,3:3*n_nodes:3)   = dNdx(1:n_nodes,2)
        endif

        if (element_identifier == 1005) then
        Bhypervector = (/Bhyper(1,1),Bhyper(2,2),Bhyper(3,3),Bhyper(1,2),Bhyper(1,3),Bhyper(2,3)/)
        call invert_small(Bhyper,Bhyperinverse,useless)
        Bhypervectorinverse = (/Bhyperinverse(1,1),Bhyperinverse(2,2),Bhyper(3,3),Bhyper(1,2),Bhyper(1,3),Bhyper(2,3)/)

        outerib = spread(Iden,dim=2,ncopies=6)*spread(Bhypervectorinverse,dim=1,ncopies=6)
        outerii = spread(Iden,dim=2,ncopies=6)*spread(Iden,dim=1,ncopies=6)
        outerbb = spread(Bhypervector,dim=2,ncopies=6)*spread(Bhypervectorinverse,dim=1,ncopies=6)

        dsubmat = 0.d0
        dsubmat(1,1) = 1
        dsubmat(2,2) = 1
        dsubmat(3,3) = 1
        dsubmat(4,4) = 0.5
        dsubmat(5,5) = 0.5
        dsubmat(6,6) = 0.5

        D = 0.d0
        !CODE FOR D GOES HERE!

        G = 0.d0
        G(1,1) = 2*B(1,1)
        G(1,4) = 2*B(1,2)
        G(1,6) = 2*B(1,3)
        G(2,2) = 2*B(2,2)
        G(2,5) = 2*B(1,2)
        G(2,8) = 2*B(2,3)
        G(3,3) = 2*B(3,3)
        G(3,7) = 2*B(1,3)
        G(3,9) = 2*B(2,3)
        G(4,1) = 2*B(1,2)
        G(4,2) = 2*B(1,2)
        G(4,4) = 2*B(2,2)
        G(4,5) = 2*B(1,1)
        G(4,6) = 2*B(2,3)
        G(4,8) = 2*B(1,3)
        G(5,1) = 2*B(1,3)
        G(5,3) = 2*B(1,3)
        G(5,4) = 2*B(2,3)
        G(5,6) = 2*B(3,3)
        G(5,7) = 2*B(1,1)
        G(5,9) = 2*B(1,2)
        G(6,2) = 2*B(2,3)
        G(6,3) = 2*B(2,3)
        G(6,5) = 2*B(1,3)
        G(6,7) = 2*B(1,2)
        G(6,8) = 2*B(3,3)
        G(6,9) = 2*B(2,2)


        Bstar = 0.d0
        Bstar(1,1:3*n_nodes-2:3) = dNdy(1:n_nodes,1)
        Bstar(2,2:3*n_nodes-1:3) = dNdy(1:n_nodes,2)
        Bstar(3,3:3*n_nodes:3)   = dNdy(1:n_nodes,3)
        Bstar(4,1:3*n_nodes-2:3) = dNdy(1:n_nodes,2)
        Bstar(5,1:3*n_nodes-2:3) = dNdy(1:n_nodes,1)
        Bstar(6,2:3*n_nodes-1:3) = dNdy(1:n_nodes,2)
        Bstar(7,3:3*n_nodes:3)   = dNdx(1:n_nodes,1)

        S = reshape(matmul(transpose(B),stress),(/3,length_dof_array/3/))
         do i = 1,n_nodes
          Pvec = reshape(spread(transpose(dNdx(i:i,1:3)),dim=2,ncopies=n_nodes),(/3*n_nodes/))
          Pmat(3*i-2:3*i,1:3*n_nodes) = spread(Pvec,dim=1,ncopies=3)
          Svec = reshape(spread(S(1:3,i:i),dim=2,ncopies=n_nodes),(/3*n_nodes/))
          Smat(3*i-2:3*i,1:3*n_nodes) = spread(Svec,dim=1,ncopies=3)
         end do
        Sigma = Pmat*transpose(Smat)

        else
        element_residual(1:3*n_nodes) = element_residual(1:3*n_nodes) -&
         matmul(transpose(B),stress)*w(kint)*determinant
        element_stiffness(1:3*n_nodes,1:3*n_nodes) = element_stiffness(1:3*n_nodes,1:3*n_nodes) &
            + matmul(transpose(B(1:6,1:3*n_nodes)),matmul(D,B(1:6,1:3*n_nodes)))*w(kint)*determinant
        endif
    end do
    return
end subroutine el_linelast_3dbasic

!========================== SUBROUTINE hypoelastic =============================================
subroutine hypoelastic(strain, dstrain, element_properties, n_properties, stress, D)
    use Types
    use ParamIO

    implicit none

    integer, intent( in )         :: n_properties                                                ! # nodes on the element

    real( prec ), intent( in )    :: strain(6), dstrain(6)             ! Strain vector contains [e11, e22, e33, 2e12, 2e13, 2e23]
    real( prec ), intent( in )    :: element_properties(n_properties)                       ! Element or material properties, stored in order listed in input file
    real( prec ), intent( out )    :: D(6,6)
    real( prec ), intent( out )   :: stress(6)



!Local variables

    real( prec )  :: K
    real( prec )  :: epsilon_e
    real( prec )  :: e_ij(6)
    real( prec )  :: e_kk
    real( prec )  :: epsilon_ij(6)
    real( prec )  :: epsilon_0
    real( prec )  :: sigma_e
    real( prec )  :: sigma_0
    real( prec )  :: nn
    real( prec )  :: E_s
    real( prec )  :: E_t
    real( prec )  :: e_dyadic_e(6,6)
    real( prec )  :: one, two
    real( prec )  :: Matrix1(6,6)
    real( prec )  :: Matrix2(6,6)

    integer :: i

!    D = 0.d0
!    D(1:3,1:3) = 0.33d0*100.d0/(1.33d0)/(1.d0-0.66d0)*(1.d0-0.33d0)
!    do i =1,4
!        D(i,i) = (1.d0-0.33d0)*100.d0/(1.33d0)/(1.d0-0.66d0)*(1.d0-0.33d0)
!    end do
!    do i = 4,6
!      D(i,i) = 50.d0/(1.33d0)
!    end do
!    stress = matmul(D,strain+dstrain)
!    return

epsilon_0 = element_properties(1)
sigma_0 = element_properties(2)
K = element_properties(3)
nn = element_properties(4)

epsilon_ij = strain + dstrain
write(6,*) epsilon_ij
e_kk = strain(1) + strain(2) + strain(3) + dstrain(1) + dstrain(2) + dstrain(3)
e_ij(1:3) = epsilon_ij(1:3) -(1.d0/3.d0)*e_kk
e_ij(4:6) = 0.5d0*epsilon_ij(4:6)
epsilon_e = dot_product(e_ij(1:3),e_ij(1:3)) + 2.d0 * dot_product(e_ij(4:6),e_ij(4:6))
epsilon_e = dsqrt( 2.d0*epsilon_e/3.d0 )
write(6,*) epsilon_e,epsilon_0
if (epsilon_e == 0.d0) then
    stress = 0.d0
    stress(1:3) = +K*e_kk
    sigma_e = sigma_0*( dsqrt( (1.d0+nn*nn)/(nn-1.d0)**2.d0 &
                    - ( nn/(nn-1.d0) - epsilon_e/epsilon_0 )**2.d0 ) - 1.d0/(nn-1.d0) )
    E_t = (sigma_0/epsilon_0)*( nn/(nn-1.d0) - epsilon_e/epsilon_0 )/ &
    dsqrt( (1.d0+nn*nn)/(nn-1.d0)**2.d0 - ( nn/(nn-1.d0) - epsilon_e/epsilon_0 )**2.d0 )
    E_s = E_t
else if (epsilon_e <= epsilon_0) then
    sigma_e = sigma_0*( dsqrt( (1.d0+nn*nn)/(nn-1.d0)**2.d0 &
                    - ( nn/(nn-1.d0) - epsilon_e/epsilon_0 )**2.d0 ) - 1.d0/(nn-1.d0) )
    stress = (2.d0/3.d0)*sigma_e*(e_ij/epsilon_e)
    stress(1:3) = stress(1:3) + K*e_kk
    E_s = sigma_e/epsilon_e
    E_t = (sigma_0/epsilon_0)*( nn/(nn-1.d0) - epsilon_e/epsilon_0 )/ &
    dsqrt( (1.d0+nn*nn)/(nn-1.d0)**2.d0 - ( nn/(nn-1.d0) - epsilon_e/epsilon_0 )**2.d0 )
else
    sigma_e = (epsilon_e/epsilon_0)**(1/nn)
    stress = (2.d0/3.d0)*sigma_e*(e_ij/epsilon_e)
    stress(1:3) = stress(1:3) + K*e_kk

    E_s = sigma_e/epsilon_e
    E_t = sigma_0/(nn*(epsilon_0*epsilon_e)**(1/nn))
    E_t = sigma_e/(nn*epsilon_e)
endif

Matrix1 = 0.d0
    one = 1.d0
    two = 2.d0
    Matrix1(1,1) = two
    Matrix1(2,2) = two
    Matrix1(3,3) = two
    Matrix1(4,4) = one
    Matrix1(5,5) = one
    Matrix1(6,6) = one

Matrix2 = 0.d0
    Matrix2(1:3,1:3) = one

e_dyadic_e = spread(e_ij,dim=2,ncopies=6)*spread(e_ij,dim=1,ncopies=6)

if (epsilon_e == 0.d0) then
    D = (E_t/3.d0)*Matrix1 + (K - ((2.d0*E_s)/9.d0))*Matrix2

else
    D = 4.d0/(9.d0*epsilon_e**2)*(E_t - E_s)*e_dyadic_e + (E_s/3)*Matrix1 + (K - ((2*E_s)/9))*Matrix2
endif



end subroutine


!========================== SUBROUTINE calculate sigma and D =============================================
!subroutine stressD(stress, D, element_properties, n_properties, Bhypervec,Bhyperinversevec,J)
!    use Types
!    use ParamIO
!
!    implicit none
!
!    integer, intent( in )         :: n_properties                                                ! # nodes on the element
!
!    real( prec ), intent( in )    :: element_properties(n_properties)                       ! Element or material properties, stored in order listed in input file
!    real( prec ), intent( out )   :: D(6,6)
!    real( prec ), intent( out )   :: stress(6)
!
!    !local variables
!
!end subroutine stressD


!==========================SUBROUTINE el_linelast_3dbasic_dynamic ==============================
subroutine el_linelast_3dbasic_dynamic(lmn, element_identifier, n_nodes, node_property_list, &           ! Input variables
    n_properties, element_properties,element_coords, length_coord_array, &                               ! Input variables
    dof_increment, dof_total, length_dof_array,  &                                                       ! Input variables
    n_state_variables, initial_state_variables, &                                                        ! Input variables
    updated_state_variables,element_residual,element_deleted)                                            ! Output variables
    use Types
    use ParamIO
    use Mesh, only : node
    use Element_Utilities, only : N => shape_functions_3D
    use Element_Utilities, only:  dNdxi => shape_function_derivatives_3D
    use Element_Utilities, only:  dNdx => shape_function_spatial_derivatives_3D
    use Element_Utilities, only : xi => integrationpoints_3D, w => integrationweights_3D
    use Element_Utilities, only : dxdxi => jacobian_3D
    use Element_Utilities, only : initialize_integration_points
    use Element_Utilities, only : calculate_shapefunctions
    use Element_Utilities, only : invert_small
    implicit none

    integer, intent( in )         :: lmn                                                    ! Element number
    integer, intent( in )         :: element_identifier                                     ! Flag identifying element type (specified in .in file)
    integer, intent( in )         :: n_nodes                                                ! # nodes on the element
    integer, intent( in )         :: n_properties                                           ! # properties for the element
    integer, intent( in )         :: length_coord_array                                     ! Total # coords
    integer, intent( in )         :: length_dof_array                                       ! Total # DOF
    integer, intent( in )         :: n_state_variables                                      ! # state variables for the element

    type (node), intent( in )     :: node_property_list(n_nodes)                  ! Data structure describing storage for nodal variables - see below
    !  type node
    !      sequence
    !      integer :: flag                          ! Integer identifier
    !      integer :: coord_index                   ! Index of first coordinate in coordinate array
    !      integer :: n_coords                      ! Total no. coordinates for the node
    !      integer :: dof_index                     ! Index of first DOF in dof array
    !      integer :: n_dof                         ! Total no. of DOF for node
    !   end type node
    !   Access these using node_property_list(k)%n_coords eg to find the number of coords for the kth node on the element

    real( prec ), intent( in )    :: element_coords(length_coord_array)                     ! Coordinates, stored as x1,(x2),(x3) for each node in turn
    real( prec ), intent( in )    :: dof_increment(length_dof_array)                        ! DOF increment, stored as du1,du2,du3,du4... for each node in turn
    real( prec ), intent( in )    :: dof_total(length_dof_array)                            ! accumulated DOF, same storage as for increment

    real( prec ), intent( in )    :: element_properties(n_properties)                       ! Element or material properties, stored in order listed in input file
    real( prec ), intent( in )    :: initial_state_variables(n_state_variables)             ! Element state variables.  Defined in this routine
               
    real( prec ), intent( inout ) :: updated_state_variables(n_state_variables)             ! State variables at end of time step
    real( prec ), intent( out )   :: element_residual(length_dof_array)                     ! Element residual force (ROW)
          
    logical, intent( inout )      :: element_deleted                                        ! Set to .true. to delete element

    ! Local Variables
    integer      :: n_points,kint

    real (prec)  ::  strain(6), dstrain(6)             ! Strain vector contains [e11, e22, e33, 2e12, 2e13, 2e23]
    real (prec)  ::  stress(6)                         ! Stress vector contains [s11, s22, s33, s12, s13, s23]
    real (prec)  ::  D(6,6)                            ! stress = D*(strain+dstrain)  (NOTE FACTOR OF 2 in shear strain)
    real (prec)  ::  B(6,length_dof_array)             ! strain = B*(dof_total+dof_increment)
    real (prec)  ::  dxidx(3,3), determinant           ! Jacobian inverse and determinant
    real (prec)  ::  x(3,length_coord_array/3)         ! Re-shaped coordinate array x(i,a) is ith coord of ath node
    real (prec)  :: E, xnu, D44, D11, D12              ! Material properties
    !
    !     Subroutine to compute element force vector for a linear elastodynamic problem
    !     El props are:

    !     element_properties(1)         Young's modulus
    !     element_properties(2)         Poisson's ratio
    
    x = reshape(element_coords,(/3,length_coord_array/3/))

    if (n_nodes == 4) n_points = 1
    if (n_nodes == 10) n_points = 4
    if (n_nodes == 8) n_points = 8
    if (n_nodes == 20) n_points = 27

    call initialize_integration_points(n_points, n_nodes, xi, w)

    element_residual = 0.d0
	
    D = 0.d0
    E = element_properties(1)
    xnu = element_properties(2)
    d44 = 0.5D0*E/(1+xnu) 
    d11 = (1.D0-xnu)*E/( (1+xnu)*(1-2.D0*xnu) )
    d12 = xnu*E/( (1+xnu)*(1-2.D0*xnu) )
    D(1:3,1:3) = d12
    D(1,1) = d11
    D(2,2) = d11
    D(3,3) = d11
    D(4,4) = d44
    D(5,5) = d44
    D(6,6) = d44
  
    !     --  Loop over integration points
    do kint = 1, n_points
        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)
        B = 0.d0
        B(1,1:3*n_nodes-2:3) = dNdx(1:n_nodes,1)
        B(2,2:3*n_nodes-1:3) = dNdx(1:n_nodes,2)
        B(3,3:3*n_nodes:3)   = dNdx(1:n_nodes,3)
        B(4,1:3*n_nodes-2:3) = dNdx(1:n_nodes,2)
        B(4,2:3*n_nodes-1:3) = dNdx(1:n_nodes,1)
        B(5,1:3*n_nodes-2:3) = dNdx(1:n_nodes,3)
        B(5,3:3*n_nodes:3)   = dNdx(1:n_nodes,1)
        B(6,2:3*n_nodes-1:3) = dNdx(1:n_nodes,3)
        B(6,3:3*n_nodes:3)   = dNdx(1:n_nodes,2)



        strain = matmul(B,dof_total)
        dstrain = matmul(B,dof_increment)
      
        stress = matmul(D,strain+dstrain)
        element_residual(1:3*n_nodes) = element_residual(1:3*n_nodes) - matmul(transpose(B),stress)*w(kint)*determinant

    end do
  
    return
end subroutine el_linelast_3dbasic_dynamic


!==========================SUBROUTINE fieldvars_linelast_3dbasic ==============================
subroutine fieldvars_linelast_3dbasic(lmn, element_identifier, n_nodes, node_property_list, &           ! Input variables
    n_properties, element_properties,element_coords,length_coord_array, &                                ! Input variables
    dof_increment, dof_total, length_dof_array,  &                                                      ! Input variables
    n_state_variables, initial_state_variables,updated_state_variables, &                               ! Input variables
    n_field_variables,field_variable_names, &                                                           ! Field variable definition
    nodal_fieldvariables)      ! Output variables
    use Types
    use ParamIO
    use Mesh, only : node
    use Element_Utilities, only : N => shape_functions_3D
    use Element_Utilities, only: dNdxi => shape_function_derivatives_3D
    use Element_Utilities, only: dNdx => shape_function_spatial_derivatives_3D
    use Element_Utilities, only: dNbardx => vol_avg_shape_function_derivatives_3D
    use Element_Utilities, only : xi => integrationpoints_3D, w => integrationweights_3D
    use Element_Utilities, only : dxdxi => jacobian_3D
    use Element_Utilities, only : initialize_integration_points
    use Element_Utilities, only : calculate_shapefunctions
    use Element_Utilities, only : invert_small
    implicit none

    integer, intent( in )         :: lmn                                                    ! Element number
    integer, intent( in )         :: element_identifier                                     ! Flag identifying element type (specified in .in file)
    integer, intent( in )         :: n_nodes                                                ! # nodes on the element
    integer, intent( in )         :: n_properties                                           ! # properties for the element
    integer, intent( in )         :: length_coord_array                                     ! Total # coords
    integer, intent( in )         :: length_dof_array                                       ! Total # DOF
    integer, intent( in )         :: n_state_variables                                      ! # state variables for the element
    integer, intent( in )         :: n_field_variables                                      ! # field variables

    type (node), intent( in )     :: node_property_list(n_nodes)                  ! Data structure describing storage for nodal variables - see below
    !  type node
    !      sequence
    !      integer :: flag                          ! Integer identifier
    !      integer :: coord_index                   ! Index of first coordinate in coordinate array
    !      integer :: n_coords                      ! Total no. coordinates for the node
    !      integer :: dof_index                     ! Index of first DOF in dof array
    !      integer :: n_dof                         ! Total no. of DOF for node
    !   end type node
    !   Access these using node_property_list(k)%n_coords eg to find the number of coords for the kth node on the element

    character (len=100), intent(in) :: field_variable_names(n_field_variables)

    real( prec ), intent( in )    :: element_coords(length_coord_array)                     ! Coordinates, stored as x1,x2,(x3) for each node in turn
    real( prec ), intent( in )    :: dof_increment(length_dof_array)                        ! DOF increment, stored as du1,du2,du3,du4... for each node in turn
    real( prec ), intent( in )    :: dof_total(length_dof_array)                            ! accumulated DOF, same storage as for increment

    real( prec ), intent( in )    :: element_properties(n_properties)                       ! Element or material properties, stored in order listed in input file
    real( prec ), intent( in )    :: initial_state_variables(n_state_variables)             ! Element state variables.  Defined in this routine
    real( prec ), intent( in )    :: updated_state_variables(n_state_variables)             ! State variables at end of time step
             
    real( prec ), intent( out )   :: nodal_fieldvariables(n_field_variables,n_nodes)        ! Nodal field variables
  
    ! Local Variables
    logical      :: strcmp
  
    integer      :: n_points,kint,k,i

    real (prec)  ::  strain(6), dstrain(6)             ! Strain vector contains [e11, e22, e33, 2e12, 2e13, 2e23]
    real (prec)  ::  stress(6)                         ! Stress vector contains [s11, s22, s33, s12, s13, s23]
    real (prec)  ::  sdev(6)                           ! Deviatoric stress
    real (prec)  ::  D(6,6)                            ! stress = D*(strain+dstrain)  (NOTE FACTOR OF 2 in shear strain)
    real (prec)  ::  B(6,length_dof_array)             ! strain = B*(dof_total+dof_increment)
    real (prec)  ::  dxidx(3,3), determinant           ! Jacobian inverse and determinant
    real (prec)  ::  x(3,length_coord_array/3)         ! Re-shaped coordinate array x(i,a) is ith coord of ath node
    real (prec)  ::  E, xnu, D44, D11, D12              ! Material properties
    real (prec)  ::  p, smises                          ! Pressure and Mises stress
    real (prec)  ::  integralargumentdNbardx           !
    real (prec)  ::  integralargumentelvol
    real (prec)  ::  elvol
    !
    !     Subroutine to compute element contribution to project element integration point data to nodes

    !     element_properties(1)         Young's modulus
    !     element_properties(2)         Poisson's ratio

    x = reshape(element_coords,(/3,length_coord_array/3/))

    if (n_nodes == 4) n_points = 1
    if (n_nodes == 10) n_points = 4
    if (n_nodes == 8) n_points = 8
    if (n_nodes == 20) n_points = 27


    write(6,*) ' hello2'

    call initialize_integration_points(n_points, n_nodes, xi, w)
    elvol = 0.d0
    dNbardx = 0.d0
    !     --  Loop over integration points
    do kint = 1, n_points
        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        integralargumentelvol = determinant*w(kint)


        dNbardx = dNbardx + dNdx*w(kint)*determinant
        elvol = elvol + integralargumentelvol
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end do
    dNbardx = dNbardx/elvol
    write(6,*) dNbardx

    nodal_fieldvariables = 0.d0
	
    D = 0.d0
    E = element_properties(1)
    xnu = element_properties(2)
    d44 = 0.5D0*E/(1+xnu) 
    d11 = (1.D0-xnu)*E/( (1+xnu)*(1-2.D0*xnu) )
    d12 = xnu*E/( (1+xnu)*(1-2.D0*xnu) )
    D(1:3,1:3) = d12
    D(1,1) = d11
    D(2,2) = d11
    D(3,3) = d11
    D(4,4) = d44
    D(5,5) = d44
    D(6,6) = d44
  
    !     --  Loop over integration points
    do kint = 1, n_points
        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)
        B = 0.d0
        B(1,1:3*n_nodes-2:3) = dNdx(1:n_nodes,1)
        B(2,2:3*n_nodes-1:3) = dNdx(1:n_nodes,2)
        B(3,3:3*n_nodes:3)   = dNdx(1:n_nodes,3)
        B(4,1:3*n_nodes-2:3) = dNdx(1:n_nodes,2)
        B(4,2:3*n_nodes-1:3) = dNdx(1:n_nodes,1)
        B(5,1:3*n_nodes-2:3) = dNdx(1:n_nodes,3)
        B(5,3:3*n_nodes:3)   = dNdx(1:n_nodes,1)
        B(6,2:3*n_nodes-1:3) = dNdx(1:n_nodes,3)
        B(6,3:3*n_nodes:3)   = dNdx(1:n_nodes,2)

        !convert B to Bbar:

        if (element_identifier == 1002) then
          do i = 1,3
          B(i,1:3*n_nodes-2:3) = B(i,1:3*n_nodes-2:3) &
          +(dNbardx(1:n_nodes,1)-dNdx(1:n_nodes,1))/3.d0
          B(i,2:3*n_nodes-1:3) = B(i,2:3*n_nodes-1:3)&
          +(dNbardx(1:n_nodes,2)-dNdx(1:n_nodes,2))/3.d0
          B(i,3:3*n_nodes:3) = B(i,3:3*n_nodes:3)&
          +(dNbardx(1:n_nodes,3)-dNdx(1:n_nodes,3))/3.d0
        end do
        endif

        strain = matmul(B,dof_total)
        dstrain = matmul(B,dof_increment)
        ! Use the same if statement to find stress for hypoelastic material

        if (n_properties==2) then
           stress = matmul(D,strain+dstrain)
 !      write(6,*) stress
        else
           write(6,*) "Hypoelastic"
           call hypoelastic(strain,dstrain,element_properties,n_properties,stress,D)
        endif

        p = sum(stress(1:3))/3.d0
        sdev = stress
        sdev(1:3) = sdev(1:3)-p
        smises = dsqrt( dot_product(sdev(1:3),sdev(1:3)) + 2.d0*dot_product(sdev(4:6),sdev(4:6)) )*dsqrt(1.5d0)

        ! In the code below the strcmp( string1, string2, nchar) function returns true if the first nchar characters in strings match
        do k = 1,n_field_variables
            if (strcmp(field_variable_names(k),'S11',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(1)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'S22',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(2)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'S33',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(3)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'S12',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(4)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'S13',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(5)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'S23',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(6)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'SMISES',6) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + smises*N(1:n_nodes)*determinant*w(kint)
            endif
        end do
 
    end do
  
    return
end subroutine fieldvars_linelast_3dbasic

