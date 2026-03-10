program test_pseudoinverse_sp  ! Scopo: esegue un test numerico di eigen_small e, dove previsto, confronta con LAPACK.
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

  integer, parameter :: N = 3
  real(sp), parameter :: THRESHOLD = 1.0e-10

  real(sp) :: A(N,N), A_save(N,N)
  real(sp) :: U_eig(N,N), U_lapack(N,N)
  real(sp) :: lambda_eig(N), lambda_lapack(N)
  real(sp) :: B(N), x_eig(N), x_lapack(N)
  real(sp) :: A_pinv_eig(N,N), A_pinv_lapack(N,N)
  real(sp) :: residual_eig, residual_lapack
  real(sp) :: work_array(64*N)
  integer :: info, i

  !---------------------------------------------------------------------------
  ! STEP 1: Creare matrice simmetrica 3x3
  !---------------------------------------------------------------------------
  write(output_unit,'(/,a)') repeat('=',70)
  write(output_unit,'(a)') 'TEST: Pseudo-inversa di Neumann per risolvere Ax=b'
  write(output_unit,'(a)') repeat('=',70)

  ! Creare matrice simmetrica A
  A(1,1) =  0.0_sp;   A(1,2) =  1.0_sp;   A(1,3) =  0.5_sp
  A(2,1) =  1.0_sp;   A(2,2) =  3.0_sp;   A(2,3) =  0.2_sp
  A(3,1) =  0.5_sp;   A(3,2) =  0.2_sp;   A(3,3) =  2.0_sp

  ! Simmetrizzazione
  A = 0.5_sp * (A + transpose(A))
  A_save = A

  write(output_unit,'(/,a)') 'Matrice A (3x3 simmetrica):'
  do i = 1, N
    write(output_unit,'(3f12.6)') A(i,:)
  end do

  !---------------------------------------------------------------------------
  ! STEP 2: Vettore di termini noti B
  !---------------------------------------------------------------------------
  B = [1.0_sp, 2.0_sp, 1.5_sp]

  write(output_unit,'(/,a)') 'Vettore B:'
  do i = 1, N
    write(output_unit,'(f12.6)') B(i)
  end do

  !---------------------------------------------------------------------------
  ! STEP 3: Diagonalizzazione con eigen_small
  !---------------------------------------------------------------------------
  write(output_unit,'(/,a)') '--- Metodo 1: eigen_small (Smith-Kopp 3x3) ---'

  U_eig = A_save
  call eig_safe(A_save, lambda_eig, U_eig)

  write(output_unit,'(/,a)') 'Autovalori (eigen_small):'
  do i = 1, N
    write(output_unit,'(f14.8)') lambda_eig(i)
  end do

  write(output_unit,'(/,a)') 'Autovettori (colonne di U_eig):'
  do i = 1, N
    write(output_unit,'(3f14.8)') U_eig(i,:)
  end do

  !---------------------------------------------------------------------------
  ! STEP 4: Diagonalizzazione con LAPACK ssyev
  !---------------------------------------------------------------------------
  write(output_unit,'(/,a)') '--- Metodo 2: LAPACK ssyev ---'

  U_lapack = A_save
  call ssyev('V', 'U', N, U_lapack, N, lambda_lapack, work_array, 64*N, info)
  if (info /= 0) then
    write(output_unit,'(a,i0)') 'ERRORE: ssyev ha restituito info=', info
    stop
  end if

  write(output_unit,'(/,a)') 'Autovalori (ssyev):'
  do i = 1, N
    write(output_unit,'(f14.8)') lambda_lapack(i)
  end do

  write(output_unit,'(/,a)') 'Autovettori (colonne di U_lapack):'
  do i = 1, N
    write(output_unit,'(3f14.8)') U_lapack(i,:)
  end do

  !---------------------------------------------------------------------------
  ! STEP 5: Costruire pseudo-inversa di Neumann con soglia (eigen_small)
  !---------------------------------------------------------------------------
  write(output_unit,'(/,a)') '--- Pseudo-inversa di Neumann (eigen_small) ---'
  write(output_unit,'(a,e12.4)') 'Soglia per autovalori:', THRESHOLD

  A_pinv_eig = construct_pseudoinverse(U_eig, lambda_eig, THRESHOLD)

  write(output_unit,'(/,a)') 'Pseudo-inversa A^+ (eigen_small):'
  do i = 1, N
    write(output_unit,'(3f14.8)') A_pinv_eig(i,:)
  end do

  !---------------------------------------------------------------------------
  ! STEP 6: Costruire pseudo-inversa di Neumann con soglia (LAPACK)
  !---------------------------------------------------------------------------
  write(output_unit,'(/,a)') '--- Pseudo-inversa di Neumann (LAPACK) ---'

  A_pinv_lapack = construct_pseudoinverse(U_lapack, lambda_lapack, THRESHOLD)

  write(output_unit,'(/,a)') 'Pseudo-inversa A^+ (LAPACK):'
  do i = 1, N
    write(output_unit,'(3f14.8)') A_pinv_lapack(i,:)
  end do

  !---------------------------------------------------------------------------
  ! STEP 7: Risolvere il sistema Ax=b usando le pseudo-inverse
  !---------------------------------------------------------------------------
  write(output_unit,'(/,a)') '--- Soluzione del sistema Ax=b ---'

  ! x = A^+ * b
  x_eig = matmul(A_pinv_eig, B)
  x_lapack = matmul(A_pinv_lapack, B)

  write(output_unit,'(/,a)') 'Soluzione x (eigen_small):'
  do i = 1, N
    write(output_unit,'(f14.8)') x_eig(i)
  end do

  write(output_unit,'(/,a)') 'Soluzione x (LAPACK):'
  do i = 1, N
    write(output_unit,'(f14.8)') x_lapack(i)
  end do

  !---------------------------------------------------------------------------
  ! STEP 8: Verificare i residui ||Ax - b||
  !---------------------------------------------------------------------------
  write(output_unit,'(/,a)') '--- Verifica dei residui ---'

  residual_eig = norm2(matmul(A_save, x_eig) - B)
  residual_lapack = norm2(matmul(A_save, x_lapack) - B)

  write(output_unit,'(/,a,e14.6)') 'Residuo ||Ax - b|| (eigen_small): ', residual_eig
  write(output_unit,'(a,e14.6)') 'Residuo ||Ax - b|| (LAPACK):      ', residual_lapack

  !---------------------------------------------------------------------------
  ! STEP 9: Differenza tra le due soluzioni
  !---------------------------------------------------------------------------
  write(output_unit,'(/,a)') '--- Confronto delle soluzioni ---'

  write(output_unit,'(/,a,e14.6)') 'Differenza max |x_eig - x_lapack|: ', maxval(abs(x_eig - x_lapack))
  write(output_unit,'(a,e14.6)') 'Norma ||x_eig - x_lapack||:       ', norm2(x_eig - x_lapack)

  write(output_unit,'(/,a)') repeat('=',70)

contains

  !---------------------------------------------------------------------------
  ! Funzione ausiliaria: costruisce la pseudo-inversa di Neumann
  !
  ! Dato: A = U * Λ * U^T
  ! Ritorna: A^+ = U * Λ^+ * U^T
  ! dove Λ^+ ha 1/λ_i sulla diagonale se |λ_i| > soglia, altrimenti 0
  !---------------------------------------------------------------------------
  function construct_pseudoinverse(U, lambda, threshold) result(A_pinv)
    real(sp), intent(in) :: U(:,:), lambda(:), threshold
    real(sp) :: A_pinv(size(U,1), size(U,1))
    real(sp) :: Lambda_inv(size(lambda)),norm_lambda
    real(sp) :: UT(size(U,1), size(U,1))
    integer :: i, n

    n = size(lambda)
    Lambda_inv = 0.0_sp
    norm_lambda = maxval(lambda)
    
    ! Costruire Λ^+ (inversa diagonale con soglia)
    do i = 1, n
      if ( lambda(i)/norm_lambda > threshold ) then
        Lambda_inv(i) = 1.0_sp / lambda(i)
      else
        Lambda_inv(i) = 0.0_sp
      end if
    end do

    ! A^+ = U * Λ^+ * U^T
       DO i=1,n
          UT(:,i) = U(:,i) * SQRT( lambda_inv(i) )
       END DO

    A_pinv = matmul(UT, transpose(UT))

  end function construct_pseudoinverse

end program test_pseudoinverse_sp