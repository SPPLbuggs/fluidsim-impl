module petsc_lib
#include <petsc/finclude/petscksp.h>
    use petscksp
    use props
    implicit none
    
    ! Petsc Variables
    Mat A1, A2, A3
    Vec b1, b2, b3, x1, x2, x3
    KSP ksp
    PetscInt Istart, Iend, ii, jj, nn(1)
    
    abstract interface
        subroutine subIn(g, i, j, n, m, dof, f, b_temp)
            use props
            type(grid), intent(in) :: g
            integer, intent(in) :: i, j, n, m, dof
            real(8), intent(in) :: f(n,m,dof)
            real(8), intent(out) :: b_temp(dof)
        end subroutine
    end interface
    
contains
    
! *** Create PETSc Objects ***
    subroutine petsc_create(g, A, b, x)
    type(grid), intent(in) :: g
    Mat A
    Vec b, x

    ! Petsc Objects A and b
    call MatCreate(comm, A, ierr)
    call MatSetSizes(A, g%nloc, g%nloc, g%nglob, g%nglob, ierr)
    call MatSetUp(A, ierr)
    call MatSetFromOptions(A, ierr)
    call MatSeqAIJSetPreallocation(A, g%dof*5, petsc_null_integer, ierr)
    call MatSetOption(A, mat_ignore_zero_entries, petsc_true, ierr)

    ! Find parallel partitioning range
    call MatGetOwnershipRange(A, Istart, Iend, ierr)

    ! Create parallel vectors
    call VecCreateMPI(comm, g%nloc, g%nglob, b, ierr)
    call VecSetFromOptions(b, ierr)
    call VecSetOption(b, vec_ignore_negative_indices, petsc_true, ierr)
    
    call VecCreateMPI(comm, g%nloc, g%nglob, x, ierr)
    call VecSetFromOptions(x, ierr)
    call VecSetOption(x, vec_ignore_negative_indices, petsc_true, ierr)

    ! Initialize Solution Vector to Zero
    call VecSet(x, 0d0, ierr)
    call VecAssemblyBegin(x, ierr)
    call VecAssemblyEnd(x, ierr)
    
    ! Create Linear Solver
    call KSPCreate(comm, ksp, ierr)
    call KSPSetOperators(ksp, A, A, ierr)
    !call KSPSetInitialGuessNonzero(ksp, PETSC_TRUE, ierr) !can be helpful
    call KSPSetTYpe(ksp, KSPIBCGS, ierr) !works well for poisson
    call KSPSetFromOptions(ksp, ierr)
    
    end subroutine

! *** Integration Step of f_pl using fEval ***
    subroutine petsc_step(g, A, b, x, f_pl, fEval, assem)
        type(grid), intent(in) :: g
        Mat A
        Vec b, x
        procedure(subIn) :: fEval
        real(8), intent(inout) :: f_pl(:,:,:)
        logical, intent(inout) :: assem
        integer :: iter, conv
        real(8) :: relErr
        
        if (assem) call petsc_create(g, A, b, x)
    
        do iter = 1, 100
            ! Assemble jacobian and RHS
            call assem_Ab(g, A, b, f_pl, fEval)
            if (assem) call MatSetOption(A, Mat_New_Nonzero_Locations, &
                                         PETSc_False, ierr)
            
            !call view(A, b)
            
            ! Check norm of residual
            call VecNorm(b, norm_2, relErr, ierr)
            if (relErr < 1d-8) exit
            
            ! Solve system:
            call KSPSetOperators(ksp, A, A, ierr)
            call KSPSolve(ksp, b, x, ierr)
            
            call KSPGetConvergedReason(ksp, conv, ierr)
            if ((my_id == 0) .and. (conv < 0)) call ksp_div(conv)
            
            ! Update variables with solution:
            call upd_soln(g, x, f_pl)
        end do
        
        assem = .False.
    end subroutine

! *** Assemble A and b ***
    subroutine assem_Ab(g, A, b, f, feval)
    type(grid), intent(in) :: g
    Mat A
    Vec b
    real(8), intent(inout) :: f(:,:,:)
    procedure(subIn) :: feval
    integer, allocatable :: cols(:), rows(:)
    real(8), allocatable :: b_temp(:), A_temp(:,:)
    integer :: i, j
    
    allocate(b_temp(g%dof), A_temp(g%dof*5, g%dof), cols(g%dof*5), rows(g%dof))
    
    ! Assemble A and b
    do j = 2, g%by+1
        do i = 2, g%bx+1
            cols = -1
            rows = g%node(i,j,:)
            
            call feval(g, i, j, g%bx+2, g%by+2, g%dof, f, b_temp)
            call jacob(g, i, j, cols, f, feval, b_temp, A_temp)
            
            ii = g%dof
            jj = g%dof*5
            call VecSetValues(b, ii, rows, -b_temp, Insert_Values, ierr)
            call MatSetValues(A, ii, rows, jj, cols, A_temp, &
                              Insert_Values, ierr)
      end do
    end do
    
    call VecAssemblyBegin(b, ierr)
    call VecAssemblyEnd(b, ierr)

    call MatAssemblyBegin(A, Mat_Final_Assembly, ierr)
    call MatAssemblyEnd(A, Mat_Final_Assembly, ierr)
    
    end subroutine
   
! *** Assemble b ***
    subroutine assem_b(g, b, f, feval)
    type(grid), intent(in) :: g
    Vec b
    real(8), intent(inout) :: f(:,:,:)
    procedure(subIn) :: feval
    integer, allocatable :: rows(:)
    real(8), allocatable :: b_temp(:)
    integer :: i, j
    
    allocate(b_temp(g%dof), rows(g%dof))
    
    ! Assemble A and b
    do j = 2, g%by+1
        do i = 2, g%bx+1
            rows = g%node(i,j,:)
            
            call feval(g, i, j, g%bx+2, g%by+2, g%dof, f, b_temp)
            
            ii = g%dof
            call VecSetValues(b, ii, rows, -b_temp, Insert_Values, ierr)
      end do
    end do
    
    call VecAssemblyBegin(b, ierr)
    call VecAssemblyEnd(b, ierr)
    end subroutine
    
! *** Numerical Jacobian ***
    subroutine jacob(g, i_loc, j_loc, cols, f, feval, b_temp, A_temp)
    type(grid), intent(in) :: g
    integer, intent(in):: i_loc, j_loc
    integer, intent(inout):: cols(:)
    procedure(subIn) :: feval
    real(8), intent(in):: b_temp(:)
    real(8), intent(inout):: f(:,:,:), A_temp(:,:)
    real(8) :: perturb, temp, b_pert(size(b_temp))
    integer :: i,j,k,d, width, k_start, k_stop
    integer, dimension(5,2):: stencil

    ! initialize
    temp = 0
    width = 0
    perturb = 1e-4
    b_pert = 0
        
    k_start = 1
    k_stop  = 5
    
    if (g%ny == 1) k_stop  = 3
    if (g%nx == 1) k_start = 3
    
    do K = k_start, k_stop
        if ((K .eq. 1) .and. (g%type_x(i_loc-1,j_loc-1) .ge. 0)) then
            width = width + 1
            stencil(width,1) = -1
            stencil(width,2) =  0
        else if ((K .eq. 2) .and. (g%type_x(i_loc-1,j_loc-1) .le. 0)) then
            width = width + 1
            stencil(width,1) = 1
            stencil(width,2) = 0
        else if (K .eq. 3) then
            width = width + 1
            stencil(width,1) = 0
            stencil(width,2) = 0
        else if ((K .eq. 4) .and. (g%type_y(i_loc-1,j_loc-1) .le. 0)) then
            width = width + 1
            stencil(width,1) = 0
            stencil(width,2) = 1
        else if ((K .eq. 5) .and. (g%type_y(i_loc-1,j_loc-1) .ge. 0)) then
            width = width + 1
            stencil(width,1) =  0
            stencil(width,2) = -1
        end if
    end do
    
    do k = 1, width
        i = i_loc + stencil(k,1)
        j = j_loc + stencil(k,2)
        
        do d = 1, g%dof
            temp = f(i,j,d)
            f(i,j,d) = f(i,j,d) + perturb
                
            call feval(g, i_loc, j_loc, g%bx+2, g%by+2, g%dof, f, b_pert)
            
            cols(g%dof * (k-1) + d) = g%node(i,j,d)
            A_temp(g%dof * (k-1) + d, :) = (b_pert - b_temp) / perturb
            
            f(i,j,d) = temp
        end do
    end do
    end subroutine
    
! *** Update Solution ***
    subroutine upd_soln(g, x, f)
    type(grid), intent(in) :: g
    Vec x
    real(8), intent(inout) :: f(:,:,:)
    integer :: i,j,d
    real(8) :: soln(1)
    
    do d = 1, g%dof
        do j = 2, g%by+1
            do i = 2, g%bx+1
                nn = g%node(i,j, d)
                call VecGetValues(x, 1, nn, soln, ierr)        
                f(i,j,d) = f(i,j,d) + soln(1)
            end do
        end do
        
        call comm_real(g%bx, g%by, f(:,:,d))
    end do
    
    end subroutine
    
! *** Destroy PETSc Objects ***
    subroutine petsc_destroy(A, b, x)
    Mat A
    Vec b, x
    call KSPDestroy(ksp,ierr)
    call VecDestroy(b,ierr)
    call VecDestroy(x,ierr)
    call MatDestroy(A,ierr)
    end subroutine

! *** View Matrix and Vector ***
    subroutine view(A, b)
    Mat A
    Vec b
    integer :: wait
    call MatView(A,PETSC_VIEWER_STDOUT_WORLD,ierr)
    call VecView(b,PETSC_VIEWER_STDOUT_WORLD,ierr)
    if (my_id == 0) read(*,*) wait
    call MPI_Barrier(comm, ierr)
    end subroutine
    
! *** Get KSP Diverged Reason
    subroutine ksp_div(val)
        integer, intent(in) :: val
        
        if (val == -2) then
            write(*,2)
        else if (val == -3) then
            write(*,3)
        else if (val == -4) then
            write(*,4)
        else if (val == -5) then
            write(*,5)
        else if (val == -6) then
            write(*,6)
        else if (val == -7) then
            write(*,7)
        else if (val == -8) then
            write(*,8)
        else if (val == -9) then
            write(*,9)
        else if (val == -10) then
            write(*,10)
        else if (val == -11) then
            write(*,11)
        end if
            
    2 format('KSP diverged due to: Null')
    3 format('KSP diverged due to: Max Iterations')
    4 format('KSP diverged due to: Divergence Tol.')
    5 format('KSP diverged due to: Breakdown')
    6 format('KSP diverged due to: Breakdown BICG')
    7 format('KSP diverged due to: Non-Symmetric')
    8 format('KSP diverged due to: Indefinite PC')
    9 format('KSP diverged due to: NaN or Infinite')
    10 format('KSP diverged due to: Indefinite Matrix')
    11 format('KSP diverged due to: PC Setup Failed')
    end subroutine
end module
