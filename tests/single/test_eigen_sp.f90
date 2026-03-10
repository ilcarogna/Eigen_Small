!==============================================================================
! test_eigen_sp.f90 – programma di test per il modulo eigen_small (sp)
!
! STRATEGIA DI VERIFICA:
!   Non vengono usati valori di riferimento scritti a mano (eccetto per
!   casi banali come la matrice identita') perche' e' facile sbagliare
!   i valori esatti.  Per ogni test si chiama LAPACK come oracolo e si
!   misura la distanza tra eigen_small e LAPACK.
!   L'unica eccezione sono i test di simmetrizzazione (uplo) dove si
!   verifica la CONSISTENZA tra uplo='F','U','L' sulla stessa matrice.
!
! Test eseguiti:
!   1)  Correttezza di symmetrize: uplo='U','L' concordano con uplo='F'
!   2)  Casi noti 2x2: identita', diagonale, densa, autovalore doppio
!   3)  Casi noti 3x3: identita', diagonale, autovalori ripetuti,
!       densa generica, quasi-singolare
!       --> in tutti i casi si confronta con LAPACK, NON con valori a mano
!   4)  Casi degeneri: matrice nulla, autovalori tripli, matrice scalata
!   5)  Consistenza delle interfacce: inplace vs safe vs uplo='U','L'
!   6)  Test statistico su N matrici SPD casuali: max errore vs LAPACK
!   7)  Benchmark: eigen_small (tre uplo) vs LAPACK
!
! Compilazione:
!   gfortran -O3 -march=native eigen_small.f90 test_eigen_sp.f90 \
!            -llapack -lblas -o test_eigen_sp && ./test_eigen_sp
!   ifx -O3 -xHost eigen_small.f90 test_eigen_sp.f90 -mkl -o test_eigen_sp
!==============================================================================
program test_eigen_sp
  ! Scopo: esegue un test numerico di eigen_small e, dove previsto, confronta con LAPACK.
  use iso_fortran_env, only: sp => real32, dp => real64, output_unit
  use eigen_small
  implicit none

  !----------------------------------------------------------------------------
  ! Interfacce LAPACK (dsyev per dp, ssyev per sp)
  !----------------------------------------------------------------------------
  interface
    subroutine dsyev(jobz, uplo, n, a, lda, w, work, lwork, info)
      import dp
      character,  intent(in)    :: jobz, uplo
      integer,    intent(in)    :: n, lda, lwork
      real(dp),   intent(inout) :: a(lda, *)
      real(dp),   intent(out)   :: w(*)
      real(dp),   intent(out)   :: work(*)
      integer,    intent(out)   :: info
    end subroutine dsyev

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

  integer, parameter :: N_STAT  = 2000000
  integer, parameter :: N_BENCH = 1000000

  ! Contatori globali per il resoconto finale
  integer :: n_ok = 0, n_fail = 0

  write(output_unit,'(/,a)') repeat('=',72)
  write(output_unit,'(a)')   '  TEST eigen_small_sp – confronto sistematico con LAPACK'
  write(output_unit,'(a,/)') repeat('=',72)

  call test_symmetrize()
  call test_casi_noti_2x2()
  call test_casi_noti_3x3()
  call test_casi_degeneri()
  call test_autovalori_negativi()
  call test_interfacce()
  call test_statistico_sp(N_STAT)
  call benchmark(N_BENCH)

  !  Resoconto finale
  write(output_unit,'(/,a)') repeat('=',72)
  write(output_unit,'(a,i0,a,i0,a)') &
    '  TOTALE: ', n_ok, ' OK,  ', n_fail, ' FAIL'
  write(output_unit,'(a,/)') repeat('=',72)

contains

  !=============================================================================
  ! UTILITA' DEL PROGRAMMA DI TEST
  !=============================================================================

  !> Separatore con titolo sezione
  subroutine sep(title)
    character(len=*), intent(in) :: title
    write(output_unit,'(/,a)')  repeat('-',72)
    write(output_unit,'(2a)')   '  ', trim(title)
    write(output_unit,'(a)')    repeat('-',72)
  end subroutine sep

  !> Controlla che val <= tol; stampa OK o FAIL e aggiorna i contatori
  subroutine chk(label, val, tol)
    character(len=*), intent(in) :: label
    real(sp),         intent(in) :: val, tol
    if (val <= tol) then
      write(output_unit,'(4x,a,a,e10.2)') 'OK   ', trim(label)//':', val
      n_ok   = n_ok   + 1
    else
      write(output_unit,'(4x,a,a,e10.2,a,e10.2)') &
        'FAIL ', trim(label)//':', val, '  tol=', tol
      n_fail = n_fail + 1
    end if
  end subroutine chk

  !> Errore di ortogonalita' max|U^T U - I|
  function ort_sp(U, n) result(e)
    integer,  intent(in) :: n
    real(sp), intent(in) :: U(n,n)
    real(sp) :: G(n,n), e
    integer :: i
    G = matmul(transpose(U), U)
    do i = 1, n; G(i,i) = G(i,i) - 1.0_sp; end do
    e = maxval(abs(G))
  end function ort_sp

  !> Chiama LAPACK SSYEV su una copia di A (n x n, sp)
  subroutine lapack_sp(A, n, lam, U)
    integer,  intent(in)  :: n
    real(sp), intent(in)  :: A(n,n)
    real(sp), intent(out) :: lam(n), U(n,n)
    real(sp) :: work(64*n)
    integer  :: info
    U = A
    call ssyev('V', 'U', n, U, n, lam, work, 64*n, info)
    if (info /= 0) then
      write(output_unit,'(4x,a,i0)') 'WARN: ssyev info = ', info
    end if
  end subroutine lapack_sp

  !> Stampa e controlla autovalori nostri vs LAPACK; controlla ortogonalita'
  !  Restituisce .true. se tutto passa
  subroutine cmp2(label, lam_our, U_our, lam_lap, tol_lam, tol_ort)
    character(len=*), intent(in) :: label
    real(sp),         intent(in) :: lam_our(2), U_our(2,2)
    real(sp),         intent(in) :: lam_lap(2)
    real(sp),         intent(in) :: tol_lam, tol_ort

    write(output_unit,'(/,4x,a)') label
    write(output_unit,'(6x,a,2f14.8)') 'eigen_small:', lam_our
    write(output_unit,'(6x,a,2f14.8)') 'LAPACK     :', lam_lap

    call chk(label//' |Dlam|', real(maxval(abs(lam_our - lam_lap)), sp), tol_lam)
    call chk(label//' |U^TU-I|', ort_sp(U_our, 2), tol_ort)
  end subroutine cmp2

  subroutine cmp3(label, lam_our, U_our, lam_lap, tol_lam, tol_ort)
    character(len=*), intent(in) :: label
    real(sp),         intent(in) :: lam_our(3), U_our(3,3)
    real(sp),         intent(in) :: lam_lap(3)
    real(sp),         intent(in) :: tol_lam, tol_ort

    write(output_unit,'(/,4x,a)') label
    write(output_unit,'(6x,a,3f14.8)') 'eigen_small:', lam_our
    write(output_unit,'(6x,a,3f14.8)') 'LAPACK     :', lam_lap

    call chk(label//' |Dlam|', real(maxval(abs(lam_our - lam_lap)), sp), tol_lam)
    call chk(label//' |U^TU-I|', ort_sp(U_our, 3), tol_ort)
  end subroutine cmp3

  !> Generatore LCG in [0,1)
  function rnd(seed) result(r)
    integer, intent(inout) :: seed
    real(sp) :: r
    seed = mod(seed * 1664525 + 1013904223, 2147483647)
    r = (seed + 1073741824) / 2147483648.0_sp
  end function rnd

  !> Matrice SPD n×n come A = L L^T
  subroutine spd_sp(A, n, seed)
    integer,  intent(in)    :: n
    real(sp), intent(out)   :: A(n,n)
    integer,  intent(inout) :: seed
    real(sp) :: L(n,n)
    integer :: i, j
    L = 0.0_sp
    do j = 1, n
      do i = j, n
        L(i,j) = rnd(seed)
      end do
    end do
    A = matmul(L, transpose(L))
  end subroutine spd_sp

  !> Converte intero in stringa
  function i2s(n) result(s)
    integer, intent(in) :: n
    character(len=20) :: s
    write(s, '(i0)') n
    s = adjustl(s)
  end function i2s

  !> Stampa riga benchmark
  subroutine print_time(label, t0, t1, n, acc)
    character(len=*), intent(in) :: label
    real(sp),         intent(in) :: t0, t1, acc
    integer,          intent(in) :: n
    write(output_unit,'(4x,a,f9.3,a,f8.2,a,e11.4,a)') &
      label, (t1-t0)*1.0e3_sp, ' ms  |', &
      (t1-t0)*1.0e9_sp/n, ' ns/mat  (acc=', acc, ')'
  end subroutine print_time

  !=============================================================================
  ! TEST 1 – CORRETTEZZA DI symmetrize (uplo)
  !
  ! La stessa matrice simmetrica viene passata con:
  !   uplo='F': matrice completa (riferimento)
  !   uplo='U': solo triangolo superiore, inferiore = 999 (spazzatura)
  !   uplo='L': solo triangolo inferiore, superiore = 999 (spazzatura)
  ! I tre set di autovalori devono essere identici a meno di eps_macchina.
  !=============================================================================

  subroutine test_symmetrize()
    real(sp), parameter :: JUNK = 999.0_sp
    real(sp), parameter :: M(3,3) = reshape( &
      [6.0_sp, 2.0_sp, 1.0_sp,  &
       2.0_sp, 5.0_sp, 3.0_sp,  &
       1.0_sp, 3.0_sp, 4.0_sp], [3,3])

    real(sp) :: M_U(3,3), M_L(3,3)
    real(sp) :: lF(3), lU(3), lL(3), U3(3,3)
    real(sp) :: lF2(2), lU2(2), lL2(2), U2(2,2)
    real(sp), parameter :: M2(2,2) = reshape( &
      [3.0_sp, 1.5_sp, 1.5_sp, 7.0_sp], [2,2])
    real(sp) :: M2_U(2,2), M2_L(2,2)

    call sep("TEST 1 – Correttezza symmetrize: uplo='U','L' == uplo='F'")

    !── 3x3 ────────────────────────────────────────────────────────────────────
    M_U = JUNK
    M_U(1,1)=M(1,1); M_U(2,2)=M(2,2); M_U(3,3)=M(3,3)
    M_U(1,2)=M(1,2); M_U(1,3)=M(1,3); M_U(2,3)=M(2,3)   ! tri sup

    M_L = JUNK
    M_L(1,1)=M(1,1); M_L(2,2)=M(2,2); M_L(3,3)=M(3,3)
    M_L(2,1)=M(2,1); M_L(3,1)=M(3,1); M_L(3,2)=M(3,2)   ! tri inf

    call eig_safe(M,   lF, U3)
    call eig_safe(M_U, lU, U3, uplo='U')
    call eig_safe(M_L, lL, U3, uplo='L')

    write(output_unit,'(/,4x,a,3f12.6)') '3x3 lam F:', lF
    write(output_unit,'(4x,a,3f12.6)')   '3x3 lam U:', lU
    write(output_unit,'(4x,a,3f12.6)')   '3x3 lam L:', lL

    call chk('3x3 |lU-lF| max', maxval(abs(lU-lF)), 1.0e-5_sp)
    call chk('3x3 |lL-lF| max', maxval(abs(lL-lF)), 1.0e-5_sp)

    !── 2x2 ────────────────────────────────────────────────────────────────────
    M2_U = JUNK;  M2_U(1,1)=M2(1,1); M2_U(2,2)=M2(2,2); M2_U(1,2)=M2(1,2)
    M2_L = JUNK;  M2_L(1,1)=M2(1,1); M2_L(2,2)=M2(2,2); M2_L(2,1)=M2(2,1)

    call eig_safe(M2,   lF2, U2)
    call eig_safe(M2_U, lU2, U2, uplo='U')
    call eig_safe(M2_L, lL2, U2, uplo='L')

    call chk('2x2 |lU-lF| max', maxval(abs(lU2-lF2)), 1.0e-6_sp)
    call chk('2x2 |lL-lF| max', maxval(abs(lL2-lF2)), 1.0e-6_sp)
  end subroutine test_symmetrize

  !=============================================================================
  ! TEST 2 – CASI NOTI 2x2 (sp)
  !  Confronto vs LAPACK su ciascuna matrice.
  !  Tolleranze: 100*eps per gli autovalori, 10*eps per l'ortogonalita'.
  !=============================================================================

  subroutine test_casi_noti_2x2()
    real(sp) :: M(2,2), lam(2), U(2,2), lL(2), UL(2,2)
    real(sp), parameter :: TOL_L = 1.0e-5_sp, TOL_O = 1.0e-6_sp

    call sep("TEST 2 – Casi noti 2x2 sp vs LAPACK")

    !── Identita' ───────────────────────────────────────────────────────────────
    M = reshape([1.0_sp,0.0_sp,0.0_sp,1.0_sp],[2,2])
    call eig_safe(M, lam, U)
    call lapack_sp(M, 2, lL, UL)
    call cmp2("Identita' 2x2", lam, U, lL, TOL_L, TOL_O)

    !── Diagonale non uniforme ──────────────────────────────────────────────────
    M = reshape([3.0_sp,0.0_sp,0.0_sp,7.0_sp],[2,2])
    call eig_safe(M, lam, U)
    call lapack_sp(M, 2, lL, UL)
    call cmp2('Diag(3,7)', lam, U, lL, TOL_L, TOL_O)

    !── Densa simmetrica ────────────────────────────────────────────────────────
    M = reshape([4.0_sp,1.5_sp,1.5_sp,6.0_sp],[2,2])
    call eig_safe(M, lam, U)
    call lapack_sp(M, 2, lL, UL)
    call cmp2('Densa [[4,1.5],[1.5,6]]', lam, U, lL, TOL_L, TOL_O)

    !── Autovalore doppio ────────────────────────────────────────────────────────
    M = reshape([5.0_sp,0.0_sp,0.0_sp,5.0_sp],[2,2])
    call eig_safe(M, lam, U)
    call lapack_sp(M, 2, lL, UL)
    call cmp2('5*I (autovalore doppio)', lam, U, lL, TOL_L, TOL_O)

    !── Matrice nulla ────────────────────────────────────────────────────────────
    M = 0.0_sp
    call eig_safe(M, lam, U)
    call lapack_sp(M, 2, lL, UL)
    call cmp2('Matrice nulla', lam, U, lL, TOL_L, TOL_O)
  end subroutine test_casi_noti_2x2

  !=============================================================================
  ! TEST 3 – CASI NOTI 3x3 (sp)
  !  Tutti i confronti sono vs LAPACK, NO valori a mano.
  !  Questo evita errori come quello della matrice "densa" con lambda=[2,4,12]
  !  che erano sbagliati (i valori corretti sono [7-sqrt(33), 4, 7+sqrt(33)]).
  !=============================================================================

  subroutine test_casi_noti_3x3()
    real(sp) :: M(3,3), lam(3), U(3,3), lL(3), UL(3,3)
    real(sp), parameter :: TOL_L = 1.0e-4_sp, TOL_O = 1.0e-6_sp

    call sep("TEST 3 – Casi noti 3x3 sp vs LAPACK (NO valori a mano)")

    !── Identita' ───────────────────────────────────────────────────────────────
    M = 0.0_sp; M(1,1)=1; M(2,2)=1; M(3,3)=1
    call eig_safe(M, lam, U)
    call lapack_sp(M, 3, lL, UL)
    call cmp3("Identita' 3x3", lam, U, lL, TOL_L, TOL_O)

    !── Diagonale disordinata ───────────────────────────────────────────────────
    M = 0.0_sp; M(1,1)=9; M(2,2)=1; M(3,3)=4   ! ordine non crescente
    call eig_safe(M, lam, U)
    call lapack_sp(M, 3, lL, UL)
    call cmp3('Diag(9,1,4)', lam, U, lL, TOL_L, TOL_O)

    !── Densa generica (questo era il caso che falliva con valori a mano) ────────
    !  [[6,4,2],[4,6,4],[2,4,6]]
    !  autovalori VERI: 7-sqrt(33) ~ 1.255, 4, 7+sqrt(33) ~ 12.745
    !  (NON 2, 4, 12 come scritto erroneamente in precedenza)
    M = reshape([6.0_sp,4.0_sp,2.0_sp, &
                 4.0_sp,6.0_sp,4.0_sp, &
                 2.0_sp,4.0_sp,6.0_sp],[3,3])
    call eig_safe(M, lam, U)
    call lapack_sp(M, 3, lL, UL)
    call cmp3('Densa [[6,4,2],[4,6,4],[2,4,6]]', lam, U, lL, TOL_L, TOL_O)
    write(output_unit,'(6x,a)') &
      '  Nota: valori corretti ~[1.2554, 4.0000, 12.7446], NON [2,4,12]'

    !── Seconda matrice densa ────────────────────────────────────────────────────
    M = reshape([4.0_sp,2.0_sp,1.0_sp, &
                 2.0_sp,6.0_sp,3.0_sp, &
                 1.0_sp,3.0_sp,5.0_sp],[3,3])
    call eig_safe(M, lam, U)
    call lapack_sp(M, 3, lL, UL)
    call cmp3('Densa [[4,2,1],[2,6,3],[1,3,5]]', lam, U, lL, TOL_L, TOL_O)

    !── Autovalori ripetuti: lambda1=lambda2 < lambda3 ─────────────────────────
    M = 0.0_sp; M(1,1)=2; M(2,2)=2; M(3,3)=7
    call eig_safe(M, lam, U)
    call lapack_sp(M, 3, lL, UL)
    call cmp3('Rip lambda1=lambda2=2, lambda3=7', lam, U, lL, TOL_L, TOL_O)

    !── Autovalori ripetuti: lambda1 < lambda2=lambda3 ─────────────────────────
    M = 0.0_sp; M(1,1)=1; M(2,2)=5; M(3,3)=5
    call eig_safe(M, lam, U)
    call lapack_sp(M, 3, lL, UL)
    call cmp3('Rip lambda1=1, lambda2=lambda3=5', lam, U, lL, TOL_L, TOL_O)

    !── Quasi-singolare (autovalore piccolo) ────────────────────────────────────
    M = reshape([1.0e-6_sp, 0.0_sp,    0.0_sp, &
                 0.0_sp,     1.0_sp,    0.5_sp, &
                 0.0_sp,     0.5_sp,    2.0_sp],[3,3])
    call eig_safe(M, lam, U)
    call lapack_sp(M, 3, lL, UL)
    write(output_unit,'(/,4x,a)')           'Quasi-singolare (lam_min ~ 1e-6)'
    write(output_unit,'(6x,a,3e16.8)') 'eigen_small:', lam
    write(output_unit,'(6x,a,3e16.8)') 'LAPACK     :', lL
    call chk('Quasi-sing |Dlam|',   maxval(abs(lam-lL)), 1.0e-6_sp)
    call chk('Quasi-sing |U^TU-I|', ort_sp(U,3),         1.0e-6_sp)

    !── Matrice con elementi molto grandi (test di scaling) ─────────────────────
    M = reshape([6.0_sp,4.0_sp,2.0_sp, &
                 4.0_sp,6.0_sp,4.0_sp, &
                 2.0_sp,4.0_sp,6.0_sp],[3,3]) * 1.0e6_sp
    call eig_safe(M, lam, U)
    call lapack_sp(M, 3, lL, UL)
    ! Tolleranza relativa: eps_sp * |lambda_max|, piu' un margine di 100x
    call cmp3('Scala grande (1e6)', lam, U, lL, &
              100.0_sp * epsilon(1.0_sp) * maxval(abs(lL)), TOL_O)

    !── Matrice con elementi molto piccoli ──────────────────────────────────────
    M = reshape([6.0_sp,4.0_sp,2.0_sp, &
                 4.0_sp,6.0_sp,4.0_sp, &
                 2.0_sp,4.0_sp,6.0_sp],[3,3]) * 1.0e-6_sp
    call eig_safe(M, lam, U)
    call lapack_sp(M, 3, lL, UL)
    ! Stessa logica: tolleranza relativa alla scala effettiva degli autovalori
    call cmp3('Scala piccola (1e-6)', lam, U, lL, &
              100.0_sp * epsilon(1.0_sp) * maxval(abs(lL)), TOL_O)
  end subroutine test_casi_noti_3x3

  !=============================================================================
  ! TEST 4 – CASI DEGENERI
  !=============================================================================

  subroutine test_casi_degeneri()
    real(sp) :: M(3,3), lam(3), U(3,3), lL(3), UL(3,3)

    call sep("TEST 4 – Casi degeneri 3x3 sp vs LAPACK")

    !── Matrice nulla ────────────────────────────────────────────────────────────
    M = 0.0_sp
    call eig_safe(M, lam, U)
    call lapack_sp(M, 3, lL, UL)
    call cmp3('Matrice nulla', lam, U, lL, 1.0e-6_sp, 1.0e-6_sp)

    !── Matrice scalare: alpha * I ───────────────────────────────────────────────
    M = 0.0_sp; M(1,1)=3; M(2,2)=3; M(3,3)=3
    call eig_safe(M, lam, U)
    call lapack_sp(M, 3, lL, UL)
    call cmp3('3*I (autovalore triplo)', lam, U, lL, 1.0e-6_sp, 1.0e-6_sp)

    !── Matrice con una riga/colonna nulla ───────────────────────────────────────
    M = 0.0_sp; M(2,2)=4; M(3,3)=9; M(2,3)=2; M(3,2)=2
    call eig_safe(M, lam, U)
    call lapack_sp(M, 3, lL, UL)
    call cmp3('Block-diag [0, [[4,2],[2,9]]]', lam, U, lL, &
              1.0e-5_sp, 1.0e-5_sp)

    !── Matrice con autovalori molto vicini ma non uguali ───────────────────────
    M = 0.0_sp
    M(1,1) = 1.0_sp;         M(2,2) = 1.0_sp + 1.0e-4_sp
    M(3,3) = 10.0_sp;        M(1,2) = 1.0e-5_sp; M(2,1) = 1.0e-5_sp
    call eig_safe(M, lam, U)
    call lapack_sp(M, 3, lL, UL)
    write(output_unit,'(/,4x,a)') 'Autovalori quasi-coincidenti (gap=1e-4)'
    write(output_unit,'(6x,a,3e18.10)') 'eigen_small:', lam
    write(output_unit,'(6x,a,3e18.10)') 'LAPACK     :', lL
    call chk('Quasi-coi |Dlam|',   maxval(abs(lam-lL)), 1.0e-5_sp)
    call chk('Quasi-coi |U^TU-I|', ort_sp(U,3),         1.0e-5_sp)
  end subroutine test_casi_degeneri

  !=============================================================================
  ! TEST 5 – MATRICI CON AUTOVALORI NEGATIVI
  !  Test di matrici simmetriche 2x2 e 3x3 con uno o piu' autovalori negativi
  !  Confronto vs LAPACK per verificare l'accuratezza anche nel caso negativo.
  !=============================================================================

  subroutine test_autovalori_negativi()
    real(sp) :: M(3,3), lam(3), U(3,3), lL(3), UL(3,3)
    real(sp) :: M2(2,2), lam2(2), U2(2,2), lL2(2), UL2(2,2)
    real(sp), parameter :: TOL_L = 1.0e-5_sp, TOL_O = 1.0e-6_sp

    call sep("TEST 5 – Matrici con autovalori negativi 2x2 e 3x3 vs LAPACK")

    !── TEST 2x2: Diagonale con due autovalori negativi ──────────────────────────
    M2 = reshape([-2.0_sp, 0.0_sp, 0.0_sp, -3.0_sp], [2,2])
    call eig_safe(M2, lam2, U2)
    call lapack_sp(M2, 2, lL2, UL2)
    call cmp2('Diag(-2,-3)', lam2, U2, lL2, TOL_L, TOL_O)

    !── TEST 2x2: Diagonale con uno positivo e uno negativo ──────────────────────
    M2 = reshape([2.0_sp, 0.0_sp, 0.0_sp, -1.0_sp], [2,2])
    call eig_safe(M2, lam2, U2)
    call lapack_sp(M2, 2, lL2, UL2)
    call cmp2('Diag(2,-1)', lam2, U2, lL2, TOL_L, TOL_O)

    !── TEST 2x2: Densa simmetrica con autovalori negativi ────────────────────────
    !  [[−1, 0.5], [0.5, −2]] ha autovalori -0.5 e -2.5
    M2 = reshape([-1.0_sp, 0.5_sp, 0.5_sp, -2.0_sp], [2,2])
    call eig_safe(M2, lam2, U2)
    call lapack_sp(M2, 2, lL2, UL2)
    call cmp2('Densa [[-1,0.5],[0.5,-2]]', lam2, U2, lL2, TOL_L, TOL_O)

    !── TEST 3x3: Diagonale con tre autovalori negativi ────────────────────────────
    M = 0.0_sp; M(1,1)=-1; M(2,2)=-2; M(3,3)=-3
    call eig_safe(M, lam, U)
    call lapack_sp(M, 3, lL, UL)
    call cmp3('Diag(-1,-2,-3)', lam, U, lL, TOL_L, TOL_O)

    !── TEST 3x3: Diagonale con autovalori misti (pos, neg) ──────────────────────
    M = 0.0_sp; M(1,1)=1; M(2,2)=-2; M(3,3)=3
    call eig_safe(M, lam, U)
    call lapack_sp(M, 3, lL, UL)
    call cmp3('Diag(1,-2,3)', lam, U, lL, TOL_L, TOL_O)

    !── TEST 3x3: Diagonale con due negativi e uno positivo ──────────────────────
    M = 0.0_sp; M(1,1)=-1; M(2,2)=-4; M(3,3)=2
    call eig_safe(M, lam, U)
    call lapack_sp(M, 3, lL, UL)
    call cmp3('Diag(-1,-4,2)', lam, U, lL, TOL_L, TOL_O)

    !── TEST 3x3: Densa simmetrica con autovalori negativi ────────────────────────
    !  [[-2, 1, 0.5], [1, -3, 0.5], [0.5, 0.5, -1]]
    M = reshape([-2.0_sp, 1.0_sp, 0.5_sp, &
                 1.0_sp, -3.0_sp, 0.5_sp, &
                 0.5_sp, 0.5_sp, -1.0_sp], [3,3])
    call eig_safe(M, lam, U)
    call lapack_sp(M, 3, lL, UL)
    write(output_unit,'(/,4x,a)') 'Densa [[-2,1,0.5],[1,-3,0.5],[0.5,0.5,-1]]'
    write(output_unit,'(6x,a,3f14.8)') 'eigen_small:', lam
    write(output_unit,'(6x,a,3f14.8)') 'LAPACK     :', lL
    call chk('Densa neg |Dlam|',   maxval(abs(lam-lL)), TOL_L)
    call chk('Densa neg |U^TU-I|', ort_sp(U,3),         TOL_O)

    !── TEST 3x3: Densa con autovalori misti ────────────────────────────────────
    !  [[2, 1, 0], [1, -2, 0.5], [0, 0.5, 3]]
    M = reshape([2.0_sp, 1.0_sp, 0.0_sp, &
                 1.0_sp, -2.0_sp, 0.5_sp, &
                 0.0_sp, 0.5_sp, 3.0_sp], [3,3])
    call eig_safe(M, lam, U)
    call lapack_sp(M, 3, lL, UL)
    write(output_unit,'(/,4x,a)') 'Densa mista [[2,1,0],[1,-2,0.5],[0,0.5,3]]'
    write(output_unit,'(6x,a,3f14.8)') 'eigen_small:', lam
    write(output_unit,'(6x,a,3f14.8)') 'LAPACK     :', lL
    call chk('Densa mista |Dlam|',   maxval(abs(lam-lL)), TOL_L)
    call chk('Densa mista |U^TU-I|', ort_sp(U,3),         TOL_O)
  end subroutine test_autovalori_negativi

  !=============================================================================
  ! TEST 6 – CONSISTENZA INTERFACCE
  !  inplace vs safe, uplo='F','U','L' devono dare lo stesso risultato.
  !=============================================================================

  subroutine test_interfacce()
    real(sp), parameter :: M_sym(3,3) = reshape( &
      [4.0_sp,2.0_sp,1.0_sp, &
       2.0_sp,6.0_sp,3.0_sp, &
       1.0_sp,3.0_sp,5.0_sp],[3,3])

    real(sp) :: M_ip(3,3), M_U(3,3), M_L(3,3)
    real(sp) :: lF(3), lIP(3), lU(3), lL(3)
    real(sp) :: UF(3,3), UU(3,3), UL_m(3,3)

    call sep("TEST 6 – Consistenza interfacce: inplace vs safe vs uplo")

    ! Versione con solo triangolo superiore valido (inferiore = 0)
    M_U = 0.0_sp
    M_U(1,1)=M_sym(1,1); M_U(2,2)=M_sym(2,2); M_U(3,3)=M_sym(3,3)
    M_U(1,2)=M_sym(1,2); M_U(1,3)=M_sym(1,3); M_U(2,3)=M_sym(2,3)

    ! Versione con solo triangolo inferiore valido (superiore = 0)
    M_L = 0.0_sp
    M_L(1,1)=M_sym(1,1); M_L(2,2)=M_sym(2,2); M_L(3,3)=M_sym(3,3)
    M_L(2,1)=M_sym(2,1); M_L(3,1)=M_sym(3,1); M_L(3,2)=M_sym(3,2)

    ! Riferimento: safe con matrice completa
    call eig_safe(M_sym, lF, UF)

    ! inplace (matrice completa, uplo default)
    M_ip = M_sym
    call eig_inplace(M_ip, lIP)

    ! safe uplo='U' e 'L'
    call eig_safe(M_U, lU, UU, uplo='U')
    call eig_safe(M_L, lL, UL_m, uplo='L')

    write(output_unit,'(/,4x,a,3f11.5)') 'safe F  lambda:', lF
    write(output_unit,'(4x,a,3f11.5)')   'inplace lambda:', lIP
    write(output_unit,'(4x,a,3f11.5)')   'safe U  lambda:', lU
    write(output_unit,'(4x,a,3f11.5)')   'safe L  lambda:', lL

    call chk('|lIP - lF|',   maxval(abs(lIP - lF)),  1.0e-6_sp)
    call chk('|lU  - lF|',   maxval(abs(lU  - lF)),  1.0e-6_sp)
    call chk('|lL  - lF|',   maxval(abs(lL  - lF)),  1.0e-6_sp)
    call chk('ort safe F',    ort_sp(UF,   3), 1.0e-6_sp)
    call chk('ort inplace',   ort_sp(M_ip, 3), 1.0e-6_sp)
    call chk('ort safe U',    ort_sp(UU,   3), 1.0e-6_sp)
    call chk('ort safe L',    ort_sp(UL_m, 3), 1.0e-6_sp)
  end subroutine test_interfacce

  !=============================================================================
  ! TEST 7 – STATISTICO sp: N matrici SPD casuali, 2x2 e 3x3
  !  Per ogni matrice: eigen_small (uplo=F,U,L) vs LAPACK.
  !  Riporta il massimo errore su autovalori e ortogonalita'.
  !=============================================================================

  subroutine test_statistico_sp(n_mat)
    integer, intent(in) :: n_mat

    real(sp) :: A2(2,2), A3(3,3), A_U(3,3), A_L(3,3)
    real(sp) :: l2(2), l3F(3), l3U(3), l3L(3)
    real(sp) :: lL2(2), lL3(3)
    real(sp) :: U2(2,2), U3(3,3), U3U(3,3), U3L(3,3), UL2(2,2), UL3(3,3)
    real(sp) :: mxL2, mxL3F, mxL3U, mxL3L
    real(sp) :: mxO2, mxO3F, mxO3U, mxO3L, mxOLap2, mxOLap3
    integer  :: seed, k

    call sep("TEST 7 – Statistico sp vs LAPACK su " // &
             trim(i2s(n_mat)) // " matrici SPD")

    mxL2=0; mxL3F=0; mxL3U=0; mxL3L=0
    mxO2=0; mxO3F=0; mxO3U=0; mxO3L=0; mxOLap2=0; mxOLap3=0
    seed = 271828

    do k = 1, n_mat
      !── 2x2 ──────────────────────────────────────────────────────────────────
      call spd_sp(A2, 2, seed)
      call eig_safe(A2, l2, U2)
      call lapack_sp(A2, 2, lL2, UL2)
      mxL2    = max(mxL2,    maxval(abs(l2  - lL2)))
      mxO2    = max(mxO2,    ort_sp(U2,  2))
      mxOLap2 = max(mxOLap2, ort_sp(UL2, 2))

      !── 3x3 uplo='F' ─────────────────────────────────────────────────────────
      call spd_sp(A3, 3, seed)
      call eig_safe(A3, l3F, U3)
      call lapack_sp(A3, 3, lL3, UL3)
      mxL3F   = max(mxL3F,   maxval(abs(l3F - lL3)))
      mxO3F   = max(mxO3F,   ort_sp(U3,  3))
      mxOLap3 = max(mxOLap3, ort_sp(UL3, 3))

      !── 3x3 uplo='U' ─────────────────────────────────────────────────────────
      A_U = 0.0_sp
      A_U(1,1)=A3(1,1); A_U(2,2)=A3(2,2); A_U(3,3)=A3(3,3)
      A_U(1,2)=A3(1,2); A_U(1,3)=A3(1,3); A_U(2,3)=A3(2,3)
      call eig_safe(A_U, l3U, U3U, uplo='U')
      mxL3U = max(mxL3U, maxval(abs(l3U - lL3)))
      mxO3U = max(mxO3U, ort_sp(U3U, 3))

      !── 3x3 uplo='L' ─────────────────────────────────────────────────────────
      A_L = 0.0_sp
      A_L(1,1)=A3(1,1); A_L(2,2)=A3(2,2); A_L(3,3)=A3(3,3)
      A_L(2,1)=A3(2,1); A_L(3,1)=A3(3,1); A_L(3,2)=A3(3,2)
      call eig_safe(A_L, l3L, U3L, uplo='L')
      mxL3L = max(mxL3L, maxval(abs(l3L - lL3)))
      mxO3L = max(mxO3L, ort_sp(U3L, 3))
    end do

    write(output_unit,'(/,4x,a)') &
      'dim  uplo   max|Dlambda vs LAPACK|   max|U^TU-I| nostra  max|U^TU-I| LAPACK'
    write(output_unit,'(4x,a)') repeat('-',70)
    write(output_unit,'(4x,a,2e22.6,e18.6)') &
      '2x2  F    ', mxL2,  mxO2,  mxOLap2
    write(output_unit,'(4x,a,2e22.6,e18.6)') &
      '3x3  F    ', mxL3F, mxO3F, mxOLap3
    write(output_unit,'(4x,a,2e22.6)') &
      '3x3  U    ', mxL3U, mxO3U
    write(output_unit,'(4x,a,2e22.6)') &
      '3x3  L    ', mxL3L, mxO3L
    write(output_unit,'(4x,a)') repeat('-',70)

    ! Soglie: confrontiamo l'errore massimo con ~100*epsilon_sp ~ 1.2e-5
    call chk('stat sp 2x2 F |Dlam|',  mxL2,  1.0e-4_sp)
    call chk('stat sp 3x3 F |Dlam|',  mxL3F, 1.0e-4_sp)
    call chk('stat sp 3x3 U |Dlam|',  mxL3U, 1.0e-4_sp)
    call chk('stat sp 3x3 L |Dlam|',  mxL3L, 1.0e-4_sp)
    call chk('stat sp 2x2   ort',      mxO2,  1.0e-5_sp)
    call chk('stat sp 3x3 F ort',      mxO3F, 1.0e-5_sp)
    call chk('stat sp 3x3 U ort',      mxO3U, 1.0e-5_sp)
    call chk('stat sp 3x3 L ort',      mxO3L, 1.0e-5_sp)
  end subroutine test_statistico_sp

  !=============================================================================
  ! TEST 8 – BENCHMARK: eigen_small (tre uplo) vs LAPACK su N_BENCH matrici
  !=============================================================================

  subroutine benchmark(n_mat)
    integer, intent(in) :: n_mat

    real(sp) :: A(3,3), A_U(3,3), U(3,3), lambda(3)
    real(sp) :: t0, t1, acc
    integer  :: seed, k

    call sep("TEST 8 – Benchmark 3x3 sp: eigen_small vs LAPACK (N=" // &
             trim(i2s(n_mat)) // ")")

    ! warm-up
    seed = 1
    do k = 1, 500
      call spd_sp(A, 3, seed)
      call eig_safe(A, lambda, U)
    end do

    !── eigen_small uplo='F' ───────────────────────────────────────────────────
    seed = 42;  acc = 0.0_sp
    call cpu_time(t0)
    do k = 1, n_mat
      call spd_sp(A, 3, seed)
      call eig_safe(A, lambda, U)
      acc = acc + lambda(2)
    end do
    call cpu_time(t1)
    call print_time("eigen_small F  ", t0, t1, n_mat, acc)

    !── eigen_small uplo='U' ───────────────────────────────────────────────────
    seed = 42;  acc = 0.0_sp
    call cpu_time(t0)
    do k = 1, n_mat
      call spd_sp(A, 3, seed)
      A_U = 0.0_sp
      A_U(1,1)=A(1,1); A_U(2,2)=A(2,2); A_U(3,3)=A(3,3)
      A_U(1,2)=A(1,2); A_U(1,3)=A(1,3); A_U(2,3)=A(2,3)
      call eig_safe(A_U, lambda, U, uplo='U')
      acc = acc + lambda(2)
    end do
    call cpu_time(t1)
    call print_time("eigen_small U  ", t0, t1, n_mat, acc)

    !── eigen_small uplo='L' ───────────────────────────────────────────────────
    seed = 42;  acc = 0.0_sp
    call cpu_time(t0)
    do k = 1, n_mat
      call spd_sp(A, 3, seed)
      A_U = 0.0_sp
      A_U(1,1)=A(1,1); A_U(2,2)=A(2,2); A_U(3,3)=A(3,3)
      A_U(2,1)=A(2,1); A_U(3,1)=A(3,1); A_U(3,2)=A(3,2)
      call eig_safe(A_U, lambda, U, uplo='L')
      acc = acc + lambda(2)
    end do
    call cpu_time(t1)
    call print_time("eigen_small L  ", t0, t1, n_mat, acc)

    !── LAPACK SSYEV ──────────────────────────────────────────────────────────
    seed = 42;  acc = 0.0_sp
    call cpu_time(t0)
    do k = 1, n_mat
      call spd_sp(A, 3, seed)
      call lapack_sp(A, 3, lambda, U)
      acc = acc + lambda(2)
    end do
    call cpu_time(t1)
    call print_time("LAPACK SSYEV   ", t0, t1, n_mat, acc)

    write(output_unit,'(/,4x,a)') &
      'Nota: il tempo include la generazione della matrice SPD.'
    write(output_unit,'(4x,a)') &
      'I tre uplo di eigen_small devono essere comparabili tra loro.'
  end subroutine benchmark

end program test_eigen_sp
