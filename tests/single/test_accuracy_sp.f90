!==============================================================================
! test_accuracy_sp – confronto di accuratezza eigen_small vs LAPACK (sp)
!
! Per ogni famiglia di matrici (SPD casuali, indefinite, due autovalori
! quasi coincidenti, autovalori su scale diverse) e per n = 2, 3 si misurano,
! su NM matrici, gli errori massimi di eigen_small e di LAPACK rispetto a un
! riferimento indipendente: una Jacobi ciclica in precisione quadrupla
! (real128) sulla stessa matrice sp.
!
! Grandezze (massimo sulle matrici di una famiglia):
!   autovalori   max|lambda - lambda_ref| / max|lambda_ref|     [eps]
!   autovettori  angolo(v, v_ref) * gap / (eps * max|lambda|)   [1]
!                (errore inevitabile ~ eps/gap per un metodo backward-stable;
!                 si salta il vettore se gap < 100*eps*|A|, mal definito)
!   residuo      max|A U - U Lambda| / max|lambda_ref|          [eps]
!
! Criterio di superamento, per ogni famiglia e grandezza:
!   errore(eigen_small) <= 1.5 * errore(LAPACK) + 1   e   <= 20
! (unita' sopra), cioe' eigen_small deve essere accurato quanto LAPACK entro
! un piccolo margine e comunque di pochi eps.
!
! Il generatore di numeri casuali e' deterministico: i risultati sono
! riproducibili e confrontabili con la tabella del README.
!
! Compilazione:
!   gfortran -O3 eigen_small.f90 test_accuracy_sp.f90 -llapack -lblas
!   ifx -O3 eigen_small.f90 test_accuracy_sp.f90 -qmkl=sequential
!==============================================================================
program test_accuracy_sp
  use iso_fortran_env, only: wp => real32, qp => real128, i8 => int64, output_unit
  use eigen_small
  implicit none

  interface
    subroutine ssyev(jobz, uplo, n, a, lda, w, work, lwork, info)
      import wp
      character, intent(in)    :: jobz, uplo
      integer,   intent(in)    :: n, lda, lwork
      real(wp),  intent(inout) :: a(lda, *)
      real(wp),  intent(out)   :: w(*)
      real(wp),  intent(out)   :: work(*)
      integer,   intent(out)   :: info
    end subroutine ssyev
  end interface

  integer, parameter :: NM = 100000          ! matrici per famiglia e dimensione
  integer, parameter :: NF = 4
  character(len=14), parameter :: fam(NF) = [character(len=14) :: &
    'SPD casuali', 'indefinite', 'cluster gap', 'scale diverse']
  ! Ampiezze (in decadi) degli intervalli casuali; dipendono dalla precisione
  real(wp), parameter :: GAP_DECADI   = 4.0_wp   ! gap relativo del cluster: 1e-1 .. 1e-(1+GR)
  real(wp), parameter :: SCALA_DECADI = 3.0_wp   ! rapporto massimo tra autovalori
  real(wp), parameter :: eps = epsilon(1.0_wp)

  real(wp) :: eO(NF,2:3), eL(NF,2:3), vO(NF,2:3), vL(NF,2:3), rO(NF,2:3), rL(NF,2:3)
  real(wp) :: A(3,3), lam(3), U(3,3), W(3,3), lr(3), work(192), Q(3,3), D(3), r(3)
  real(wp) :: gap, ao, al, lmax
  real(qp) :: Aq(3,3), Vq(3,3), lq(3)
  integer(i8) :: state
  integer :: f, n, m, i, j, k, kk, info, n_ok, n_fail

  n_ok = 0;  n_fail = 0
  eO = 0; eL = 0; vO = 0; vL = 0; rO = 0; rL = 0
  state = 20240607_i8

  write(output_unit,'(/,a)') repeat('=',78)
  write(output_unit,'(a,i0,a)') '  TEST accuratezza sp: eigen_small vs LAPACK vs riferimento quad (', NM, ' matrici/famiglia)'
  write(output_unit,'(a,/)') repeat('=',78)

  do n = 2, 3
    do f = 1, NF
      do m = 1, NM
        A = 0.0_wp
        select case (f)
        case (1)      ! SPD: L L^T con L triangolare inferiore casuale
          Q = 0.0_wp
          do j = 1, n; do i = j, n; Q(i,j) = rnd(); end do; end do
          A(1:n,1:n) = matmul(Q(1:n,1:n), transpose(Q(1:n,1:n)))
        case (2)      ! simmetrica con elementi in [-1,1]
          do j = 1, n; do i = j, n
            A(i,j) = 2.0_wp*rnd() - 1.0_wp;  A(j,i) = A(i,j)
          end do; end do
        case (3)      ! due autovalori vicini [d, d+g] e uno separato
          r = [rnd(), rnd(), rnd()]
          D = [0.1_wp*r(1), 0.0_wp, 1.0_wp]
          D(2) = D(1) + 10.0_wp**(-1.0_wp - GAP_DECADI*r(2))
          A = congruenza(n, D, r)
        case (4)      ! autovalori di modulo 10^-SCALA_DECADI .. 1, segno misto
          r = [rnd(), rnd(), rnd()]
          D = [10.0_wp**(-SCALA_DECADI*r(1)), -10.0_wp**(-SCALA_DECADI*r(2)), 1.0_wp]
          if (r(3) > 0.5_wp) D(2) = -D(2)
          A = congruenza(n, D, r)
        end select

        call eig_safe(A(1:n,1:n), lam(1:n), U(1:n,1:n))
        W = A
        call ssyev('V', 'U', n, W, 3, lr, work, 192, info)
        Aq = real(A, qp)
        call jacobi(n, Aq, Vq, lq)

        lmax = real(max(abs(lq(1)), abs(lq(n))), wp)
        eO(f,n) = max(eO(f,n), maxval(abs(real(real(lam(1:n),qp) - lq(1:n), wp))) / lmax)
        eL(f,n) = max(eL(f,n), maxval(abs(real(real(lr(1:n),qp)  - lq(1:n), wp))) / lmax)
        rO(f,n) = max(rO(f,n), maxval(abs(matmul(A(1:n,1:n), U(1:n,1:n)) &
                                - U(1:n,1:n)*spread(lam(1:n),1,n))) / lmax)
        rL(f,n) = max(rL(f,n), maxval(abs(matmul(A(1:n,1:n), W(1:n,1:n)) &
                                - W(1:n,1:n)*spread(lr(1:n),1,n))) / lmax)
        do k = 1, n
          gap = huge(1.0_wp)
          do kk = 1, n
            if (kk /= k) gap = min(gap, abs(real(lq(k) - lq(kk), wp)))
          end do
          if (gap > 100.0_wp*eps*lmax) then
            ao = angolo(n, U(:,k), Vq(:,k))
            al = angolo(n, W(:,k), Vq(:,k))
            vO(f,n) = max(vO(f,n), ao*gap/(eps*lmax))
            vL(f,n) = max(vL(f,n), al*gap/(eps*lmax))
          end if
        end do
      end do
    end do
  end do

  write(output_unit,'(a)') &
    '  famiglia        n |  autovalori [eps]     |  autovettori [eps/gap] |  residuo [eps]'
  write(output_unit,'(a)') &
    '                    |  eig_small   LAPACK   |  eig_small   LAPACK    |  eig_small   LAPACK'
  write(output_unit,'(2x,a)') repeat('-',76)
  do n = 2, 3
    do f = 1, NF
      write(output_unit,'(2x,a14,i2,a,2f10.1,a,2f10.1,a,2f10.1)') fam(f), n, ' |', &
        eO(f,n)/eps, eL(f,n)/eps, '  |', vO(f,n), vL(f,n), '   |', rO(f,n)/eps, rL(f,n)/eps
    end do
  end do
  write(output_unit,'(/,a)') '  Verifica (eigen_small <= 1.5*LAPACK + 1  e  <= 20):'
  do n = 2, 3
    do f = 1, NF
      call verifica(trim(fam(f))//' '//char(48+n)//'x'//char(48+n)//' autovalori',  eO(f,n)/eps, eL(f,n)/eps)
      call verifica(trim(fam(f))//' '//char(48+n)//'x'//char(48+n)//' autovettori', vO(f,n),     vL(f,n))
      call verifica(trim(fam(f))//' '//char(48+n)//'x'//char(48+n)//' residuo',     rO(f,n)/eps, rL(f,n)/eps)
    end do
  end do

  write(output_unit,'(/,a)') repeat('=',78)
  write(output_unit,'(a,i0,a,i0,a)') '  TOTALE: ', n_ok, ' OK,  ', n_fail, ' FAIL'
  write(output_unit,'(a,/)') repeat('=',78)
  if (n_fail > 0) error stop 1

contains

  !> Generatore LCG deterministico a 64 bit, restituisce un reale in [0,1)
  function rnd() result(x)
    real(wp) :: x
    state = mod(state * 1103515245_i8 + 12345_i8, 2147483648_i8)
    x = real(state, wp) / 2147483648.0_wp
  end function rnd

  !> Aggiorna i contatori: ok se ours <= 1.5*lap + 1 e ours <= 20
  subroutine verifica(label, ours, lap)
    character(len=*), intent(in) :: label
    real(wp), intent(in) :: ours, lap
    if (ours <= 1.5_wp*lap + 1.0_wp .and. ours <= 20.0_wp) then
      n_ok = n_ok + 1
    else
      write(output_unit,'(4x,a,a,2f9.2)') 'FAIL ', trim(label)//':', ours, lap
      n_fail = n_fail + 1
    end if
  end subroutine verifica

  !> A = Q diag(D) Q^T con Q rotazione casuale n x n, resa esattamente simmetrica
  function congruenza(n, D, r) result(B)
    integer, intent(in) :: n
    real(wp), intent(in) :: D(3), r(3)
    real(wp) :: B(3,3), R3(3,3), c(3), s(3), t(3)
    t = 6.2831853_wp * r;  c = cos(t);  s = sin(t)
    R3 = 0.0_wp
    if (n == 2) then
      R3(1:2,1:2) = reshape([c(1), s(1), -s(1), c(1)], [2,2])
    else
      R3 = reshape([ c(1)*c(2),                 s(1)*c(2),                -s(2),    &
                     c(1)*s(2)*s(3)-s(1)*c(3),  s(1)*s(2)*s(3)+c(1)*c(3),  c(2)*s(3), &
                     c(1)*s(2)*c(3)+s(1)*s(3),  s(1)*s(2)*c(3)-c(1)*s(3),  c(2)*c(3) ], [3,3])
    end if
    B = 0.0_wp
    B(1:n,1:n) = matmul(R3(1:n,1:n) * spread(D(1:n), 1, n), transpose(R3(1:n,1:n)))
    B(1:n,1:n) = 0.5_wp * (B(1:n,1:n) + transpose(B(1:n,1:n)))
  end function congruenza

  !> Seno dell'angolo tra v (wp) e w (qp), vettori unitari, indipendente dal segno
  function angolo(n, v, w) result(a)
    integer, intent(in) :: n
    real(wp), intent(in) :: v(3)
    real(qp), intent(in) :: w(3)
    real(wp) :: a
    real(qp) :: d, p(3)
    d = sum(real(v(1:n), qp) * w(1:n))
    p(1:n) = real(v(1:n), qp) - d * w(1:n)
    a = real(sqrt(sum(p(1:n)**2)), wp)
  end function angolo

  !> Jacobi ciclica in precisione quadrupla: autovalori crescenti e autovettori
  subroutine jacobi(n, B, V, l)
    integer, intent(in)    :: n
    real(qp), intent(inout) :: B(3,3)
    real(qp), intent(out)   :: V(3,3), l(3)
    integer  :: sweep, p, q, i
    real(qp) :: th, t, c, s, tmp(3), x
    V = 0.0_qp
    do i = 1, 3; V(i,i) = 1.0_qp; end do
    do sweep = 1, 30
      do p = 1, n-1
        do q = p+1, n
          if (abs(B(p,q)) < 1.0e-60_qp) cycle
          th = (B(q,q) - B(p,p)) / (2.0_qp*B(p,q))
          t  = sign(1.0_qp, th) / (abs(th) + sqrt(th*th + 1.0_qp))
          c  = 1.0_qp / sqrt(t*t + 1.0_qp);  s = t*c
          tmp = B(:,p);  B(:,p) = c*tmp - s*B(:,q);  B(:,q) = s*tmp + c*B(:,q)
          tmp = B(p,:);  B(p,:) = c*tmp - s*B(q,:);  B(q,:) = s*tmp + c*B(q,:)
          tmp = V(:,p);  V(:,p) = c*tmp - s*V(:,q);  V(:,q) = s*tmp + c*V(:,q)
        end do
      end do
    end do
    l = 0.0_qp
    do i = 1, n; l(i) = B(i,i); end do
    do p = 1, n-1
      do q = p+1, n
        if (l(q) < l(p)) then
          x = l(p);  l(p) = l(q);  l(q) = x
          tmp = V(:,p);  V(:,p) = V(:,q);  V(:,q) = tmp
        end if
      end do
    end do
  end subroutine jacobi

end program test_accuracy_sp
