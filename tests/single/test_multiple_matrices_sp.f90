program test_multiple_matrices_sp  ! Scopo: esegue un test numerico di eigen_small e, dove previsto, confronta con LAPACK.
  ! Scopo: esegue un test numerico di eigen_small e, dove previsto, confronta con LAPACK.
  use iso_fortran_env, only: sp => real32, output_unit
  use eigen_small
  implicit none

  interface
    subroutine ssyev(jobz, uplo, n, a, lda, w, work, lwork, info)
      import sp
      character,  intent(in)    :: jobz, uplo
      integer,    intent(in)    :: n, lda, lwork
      real(sp),   intent(inout) :: a(lda, *)
      real(sp),   intent(out)   :: w(*)
      real(sp),   intent(out)   :: work(*)
      integer,    intent(out)   :: info
    end subroutine ssyev
  end interface

  integer, parameter :: MATRIX_SIZE = 3
  integer, parameter :: N_MATRICES = 1000  ! Numero di matrici (può essere modificato)
  real(sp), parameter :: THRESHOLD = 1.0e-10

  real(sp), allocatable :: A(:,:,:)           ! A(3,3,n_matrices)
  real(sp), allocatable :: B(:,:)             ! B(3,n_matrices) - vettori termini noti
  real(sp), allocatable :: X_eig(:,:)         ! X_eig(3,n_matrices) - soluzioni eig
  real(sp), allocatable :: X_lapack(:,:)      ! X_lapack(3,n_matrices) - soluzioni lapack
  real(sp), allocatable :: residuals_eig(:)   ! residui eig
  real(sp), allocatable :: residuals_lapack(:) ! residui lapack

  real(sp) :: U_eig(MATRIX_SIZE,MATRIX_SIZE), lambda_eig(MATRIX_SIZE)
  real(sp) :: U_lapack(MATRIX_SIZE,MATRIX_SIZE), lambda_lapack(MATRIX_SIZE)
  real(sp) :: A_save(MATRIX_SIZE,MATRIX_SIZE)
  real(sp) :: A_pinv_eig(MATRIX_SIZE,MATRIX_SIZE), A_pinv_lapack(MATRIX_SIZE,MATRIX_SIZE)
  real(sp) :: work_array(64*MATRIX_SIZE)
  real(sp) :: random_val
  integer :: info, i, j, k, seed

  !---------------------------------------------------------------------------
  ! Allocazione degli array
  !---------------------------------------------------------------------------
  allocate(A(MATRIX_SIZE, MATRIX_SIZE, N_MATRICES))
  allocate(B(MATRIX_SIZE, N_MATRICES))
  allocate(X_eig(MATRIX_SIZE, N_MATRICES))
  allocate(X_lapack(MATRIX_SIZE, N_MATRICES))
  allocate(residuals_eig(N_MATRICES))
  allocate(residuals_lapack(N_MATRICES))

  !---------------------------------------------------------------------------
  ! Intestazione
  !---------------------------------------------------------------------------
  write(output_unit,'(/,a)') repeat('=',80)
  write(output_unit,'(a,i0,a)') 'TEST: Risoluzione di ', N_MATRICES, ' sistemi lineari simmetrici (3x3)'
  write(output_unit,'(a,i0)') 'Numero di matrici (N): ', N_MATRICES
  write(output_unit,'(a)') repeat('=',80)

  !---------------------------------------------------------------------------
  ! Generazione di matrici simmetriche casuali e vettori b
  !---------------------------------------------------------------------------
  write(output_unit,'(/,a)') '--- Generazione matrici simmetriche casuali ---'
  
  ! Inizializzazione generatore di numeri casuali
  seed = 12345
  call random_seed(put=[(i, i=1, 33)])

  do k = 1, N_MATRICES
    ! Generare matrice casuali 3x3 non simmetrica
    do i = 1, MATRIX_SIZE
      do j = 1, MATRIX_SIZE
        call random_number(random_val)
        A(i,j,k) = 2.0_sp * random_val - 1.0_sp  ! Valori tra -1 e 1
      end do
    end do

    ! Simmetrizzare: A = 0.5*(A + A^T)
    do i = 1, MATRIX_SIZE
      do j = i+1, MATRIX_SIZE
        A(i,j,k) = 0.5_sp * (A(i,j,k) + A(j,i,k))
        A(j,i,k) = A(i,j,k)
      end do
    end do
! facciamola semidefinita positiva moltiplicando A per la sua trasposta
    A(:,:,k) = matmul(A(:,:,k), transpose(A(:,:,k)))
    
    ! Generare vettore b casuale
    do i = 1, MATRIX_SIZE
      call random_number(random_val)
      B(i,k) = 2.0_sp * random_val - 1.0_sp  ! Valori tra -1 e 1
    end do
  end do

  write(output_unit,'(a,i0,a)') 'Generati ', N_MATRICES, ' matrici simmetriche 3x3 e relativi vettori b'

  !---------------------------------------------------------------------------
  ! Ciclo principale: risolvere ogni sistema
  !---------------------------------------------------------------------------
  write(output_unit,'(/,a)') repeat('-',80)

  do k = 1, N_MATRICES
    write(output_unit,'(/,a,i0,a,i0)') 'SISTEMA ', k, ' di ', N_MATRICES
    write(output_unit,'(a)') repeat('-',80)

    ! Salvare la matrice originale
    A_save = A(:,:,k)

    !-----------------------------------------------------------------------
    ! Stampa della matrice e del vettore b
    !-----------------------------------------------------------------------
    write(output_unit,'(a)') 'Matrice A:'
    do i = 1, MATRIX_SIZE
      write(output_unit,'(3f12.6)') A_save(i,:)
    end do

    write(output_unit,'(/,a)') 'Vettore b:'
    do i = 1, MATRIX_SIZE
      write(output_unit,'(f12.6)') B(i,k)
    end do

    !-----------------------------------------------------------------------
    ! Metodo 1: eig_inplace
    !-----------------------------------------------------------------------
    write(output_unit,'(/,a)') '  Metodo 1: eig_inplace'
    U_eig = A_save
    call eig_inplace(U_eig, lambda_eig)

    write(output_unit,'(a)') '    Autovalori:'
    do i = 1, MATRIX_SIZE
      write(output_unit,'(a,f14.8)') '      ', lambda_eig(i)
    end do

    ! Costruire pseudo-inversa
    A_pinv_eig = construct_pseudoinverse(U_eig, lambda_eig, THRESHOLD)

    ! Risolvere sistema
    X_eig(:,k) = matmul(A_pinv_eig, B(:,k))
    residuals_eig(k) = norm2(matmul(A_save, X_eig(:,k)) - B(:,k))

    write(output_unit,'(a)') '    Soluzione x:'
    do i = 1, MATRIX_SIZE
      write(output_unit,'(a,f14.8)') '      ', X_eig(i,k)
    end do
    write(output_unit,'(a,e14.6)') '    Residuo ||Ax - b||: ', residuals_eig(k)

    !-----------------------------------------------------------------------
    ! Metodo 2: LAPACK ssyev
    !-----------------------------------------------------------------------
    write(output_unit,'(/,a)') '  Metodo 2: LAPACK ssyev'
    U_lapack = A_save
    call ssyev('V', 'U', MATRIX_SIZE, U_lapack, MATRIX_SIZE, lambda_lapack, work_array, 64*MATRIX_SIZE, info)
    
    if (info /= 0) then
      write(output_unit,'(a,i0)') '    ERRORE: ssyev ha restituito info=', info
      cycle
    end if

    write(output_unit,'(a)') '    Autovalori:'
    do i = 1, MATRIX_SIZE
      write(output_unit,'(a,f14.8)') '      ', lambda_lapack(i)
    end do

    ! Costruire pseudo-inversa
    A_pinv_lapack = construct_pseudoinverse(U_lapack, lambda_lapack, THRESHOLD)

    ! Risolvere sistema
    X_lapack(:,k) = matmul(A_pinv_lapack, B(:,k))
    residuals_lapack(k) = norm2(matmul(A_save, X_lapack(:,k)) - B(:,k))

    write(output_unit,'(a)') '    Soluzione x:'
    do i = 1, MATRIX_SIZE
      write(output_unit,'(a,f14.8)') '      ', X_lapack(i,k)
    end do
    write(output_unit,'(a,e14.6)') '    Residuo ||Ax - b||: ', residuals_lapack(k)

    !-----------------------------------------------------------------------
    ! Confronto
    !-----------------------------------------------------------------------
    write(output_unit,'(/,a,e14.6)') '  Differenza max soluzioni: ', maxval(abs(X_eig(:,k) - X_lapack(:,k)))

  end do

  !---------------------------------------------------------------------------
  ! Riepilogo finale
  !---------------------------------------------------------------------------
  write(output_unit,'(/,a)') repeat('=',80)
  write(output_unit,'(a)') 'RIEPILOGO GENERALE'
  write(output_unit,'(a)') repeat('=',80)

  write(output_unit,'(/,a)') 'Residui eig_inplace:'
  do k = 1, N_MATRICES
    write(output_unit,'(a,i0,a,e14.6)') '  Sistema ', k, ': ', residuals_eig(k)
  end do

  write(output_unit,'(/,a)') 'Residui LAPACK:'
  do k = 1, N_MATRICES
    write(output_unit,'(a,i0,a,e14.6)') '  Sistema ', k, ': ', residuals_lapack(k)
  end do

  write(output_unit,'(/,a)') 'Statistiche:'
  write(output_unit,'(a,e14.6)') '  Residuo medio (eig_inplace):    ', sum(residuals_eig)/N_MATRICES
  write(output_unit,'(a,e14.6)') '  Residuo max (eig_inplace):      ', maxval(residuals_eig)
  write(output_unit,'(a,e14.6)') '  Residuo medio (LAPACK):         ', sum(residuals_lapack)/N_MATRICES
  write(output_unit,'(a,e14.6)') '  Residuo max (LAPACK):           ', maxval(residuals_lapack)

  write(output_unit,'(/,a)') repeat('=',80)

  ! Deallocazione
  deallocate(A, B, X_eig, X_lapack, residuals_eig, residuals_lapack)

contains

  !---------------------------------------------------------------------------
  ! Funzione ausiliaria: costruisce la pseudo-inversa di Neumann
  !---------------------------------------------------------------------------
  function construct_pseudoinverse(U, lambda, threshold) result(A_pinv)
    real(sp), intent(in) :: U(:,:), lambda(:), threshold
    real(sp) :: A_pinv(size(U,1), size(U,1))
    real(sp) :: Lambda_inv(size(lambda)), norm_lambda
    real(sp) :: UT(size(U,1), size(U,1))
    integer :: i, n

    n = size(lambda)
    Lambda_inv = 0.0_sp
    norm_lambda = maxval(abs(lambda))

    ! Costruire Λ^+ (inversa diagonale con soglia)
    do i = 1, n
      if (norm_lambda > 0.0_sp .and. abs(lambda(i))/norm_lambda > threshold) then
        Lambda_inv(i) = 1.0_sp / lambda(i)
      else
        Lambda_inv(i) = 0.0_sp
      end if
    end do

    ! A^+ = U * Λ^+ * U^T
    do i = 1, n
      UT(:,i) = U(:,i) * sqrt(abs(Lambda_inv(i)))
    end do

    A_pinv = matmul(UT, transpose(UT))

  end function construct_pseudoinverse

end program test_multiple_matrices_sp