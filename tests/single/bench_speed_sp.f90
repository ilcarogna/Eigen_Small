!==============================================================================
! bench_speed_sp.f90 -- Benchmark di velocita' eigen_small vs LAPACK ssyev (3x3)
!
! (versione singola precisione; timer, statistiche e checksum in doppia)
!
! Misura il tempo per matrice di:
!   - eig_inplace  (eigen_small)
!   - eig_safe     (eigen_small)
!   - ssyev('V')   (LAPACK, workspace preallocato con query ottimale)
!
! La produzione delle matrici di test NON e' temporizzata: tutte le N matrici
! sono generate prima in A(3,3,N). Dentro i cicli temporizzati non c'e' alcun
! I/O. Poiche' eig_inplace e ssyev sovrascrivono l'input, ogni iterazione
! copia A(:,:,k) in una matrice di lavoro: il costo di questa copia viene
! misurato da un ciclo di baseline (stessa copia + stesso checksum, nessun
! calcolo) e sottratto.
!
! Robustezza rispetto a deriva della frequenza della CPU, carico della
! macchina e ordine di esecuzione:
!   - un giro di riscaldamento non misurato prima delle ripetizioni;
!   - i metodi sono INTERCALATI: ogni ripetizione esegue baseline e i tre
!     metodi una volta, in ordine ruotato a ogni ripetizione, cosi' ogni
!     metodo vede in media le stesse condizioni;
!   - il tempo netto di un metodo e' (tempo - baseline) della STESSA
!     ripetizione (misura accoppiata);
!   - si riportano mediana e minimo dei tempi netti, piu' la dispersione
!     (max-min)/mediana: se e' grande le misure non sono affidabili;
!   - lo speedup e' il rapporto delle mediane.
! I checksum sugli autovalori evitano l'eliminazione del calcolo da parte
! del compilatore e verificano che i metodi diano lo stesso risultato.
!==============================================================================
program bench_speed_sp
  use iso_fortran_env, only: sp => real32, dp => real64, int64, output_unit
  use eigen_small
  implicit none

  integer, parameter :: N     = 1000000   ! numero di matrici
  integer, parameter :: NREP  = 15        ! ripetizioni misurate
  integer, parameter :: NMETH = 4         ! 0=baseline 1=inplace 2=safe 3=ssyev

  real(sp), allocatable :: A(:,:,:)
  real(sp) :: U(3,3), lam(3), Uo(3,3), r, work_q(1)
  real(sp), allocatable :: work(:)
  integer  :: lwork, info, i, j, k, rep, q, m
  integer(int64) :: t0, t1, rate
  real(dp) :: tt(0:NMETH-1, NREP)       ! tempi grezzi [s]
  real(dp) :: net(NREP)
  real(dp) :: s(0:NMETH-1)              ! checksum dell'ultima esecuzione
  real(dp) :: med(0:NMETH-1), mn(0:NMETH-1), spread(0:NMETH-1)
  character(len=16), parameter :: nome(0:NMETH-1) = &
       [character(len=16) :: 'baseline (copia)', 'eig_inplace', 'eig_safe', 'LAPACK ssyev']

  ! --- Produzione delle matrici (fuori dal tempo misurato) -------------------
  allocate(A(3,3,N))
  call random_seed(put=[(i, i=1, 33)])
  do k = 1, N
    do j = 1, 3
      do i = 1, 3
        call random_number(r)
        A(i,j,k) = 2.0_sp*r - 1.0_sp
      end do
    end do
    A(:,:,k) = 0.5_sp*(A(:,:,k) + transpose(A(:,:,k)))   ! simmetrica
  end do

  ! --- Workspace ottimale per ssyev -----------------------------------------
  U = A(:,:,1)
  call ssyev('V', 'U', 3, U, 3, lam, work_q, -1, info)
  lwork = max(int(work_q(1)), 3*3 - 1)
  allocate(work(lwork))

  call system_clock(count_rate=rate)

  ! --- Riscaldamento (non misurato): un giro completo di tutti i metodi ------
  do m = 0, NMETH-1
    call run(m)
  end do

  ! --- Misure intercalate, ordine ruotato a ogni ripetizione -----------------
  do rep = 1, NREP
    do q = 0, NMETH-1
      m = mod(q + rep - 1, NMETH)
      call system_clock(t0)
      call run(m)
      call system_clock(t1)
      tt(m, rep) = real(t1-t0, dp)/real(rate, dp)
    end do
  end do

  ! --- Statistiche -----------------------------------------------------------
  do m = 0, NMETH-1
    if (m == 0) then
      net = tt(0,:)
    else
      net = tt(m,:) - tt(0,:)            ! misura accoppiata
    end if
    call ordina(net)
    med(m) = mediana(net)
    mn(m)  = net(1)
    spread(m) = (net(NREP) - net(1)) / med(m)
  end do

  ! --- Report ----------------------------------------------------------------
  write(output_unit,'(a,i0,a,i0,a)') 'Benchmark 3x3 simmetriche: N = ', N, &
        ', ripetizioni = ', NREP, ' (intercalate, ordine ruotato)'
  write(output_unit,'(a)') 'Tempo per matrice, baseline (copia) sottratta:'
  write(output_unit,'(a,t30,a,t44,a,t58,a)') '  metodo', 'mediana [ns]', 'minimo [ns]', 'dispers.'
  do m = 0, NMETH-1
    write(output_unit,'(2x,a,t30,f9.1,t46,f9.1,t58,f7.1,a)') nome(m), &
          1.0e9_dp*med(m)/N, 1.0e9_dp*mn(m)/N, 100.0_dp*spread(m), ' %'
  end do
  write(output_unit,'(/,a,f7.2,a)') 'Speedup eig_inplace / ssyev (mediane): ', med(3)/med(1), ' x'
  write(output_unit,'(a,f7.2,a)')   'Speedup eig_safe    / ssyev (mediane): ', med(3)/med(2), ' x'
  write(output_unit,'(a,f7.2,a)')   'Rapporto eig_safe / eig_inplace      : ', med(2)/med(1), ' x'
  if (maxval(spread(1:)) > 0.10_dp) write(output_unit,'(/,a)') &
    'ATTENZIONE: dispersione > 10 %: macchina rumorosa, misure poco affidabili.'
  write(output_unit,'(/,a)') 'Controllo (somma autovalori, deve coincidere):'
  write(output_unit,'(a,es22.14)') '  eig_inplace : ', s(1)
  write(output_unit,'(a,es22.14)') '  eig_safe    : ', s(2)
  write(output_unit,'(a,es22.14)') '  ssyev       : ', s(3)
  write(output_unit,'(a,es10.2)')  '  |inplace-ssyev|/|ssyev| : ', &
        abs(s(1)-s(3))/max(abs(s(3)), tiny(1.0_dp))
  if (s(0) > huge(1.0_dp)) write(output_unit,*) s(0)   ! anti-eliminazione baseline

contains

  !> Esegue il ciclo su tutte le N matrici per il metodo m; checksum in s(m)
  subroutine run(m)
    integer, intent(in) :: m
    integer :: kk
    real(dp) :: acc
    acc = 0.0_dp
    select case (m)
    case (0)   ! baseline: copia + checksum, nessun calcolo
      do kk = 1, N
        U = A(:,:,kk)
        lam = U(:,1)
        acc = acc + lam(1) + lam(2) + lam(3)
      end do
    case (1)
      do kk = 1, N
        U = A(:,:,kk)
        call eig_inplace(U, lam, 'U')
        acc = acc + lam(1) + lam(2) + lam(3) + U(1,1)*0.0_sp
      end do
    case (2)
      do kk = 1, N
        U = A(:,:,kk)
        call eig_safe(U, lam, Uo, 'U')
        acc = acc + lam(1) + lam(2) + lam(3) + Uo(1,1)*0.0_sp
      end do
    case (3)
      do kk = 1, N
        U = A(:,:,kk)
        call ssyev('V', 'U', 3, U, 3, lam, work, lwork, info)
        acc = acc + lam(1) + lam(2) + lam(3) + U(1,1)*0.0_sp
      end do
    end select
    s(m) = acc
  end subroutine run

  !> Ordinamento crescente per inserzione (NREP piccolo)
  subroutine ordina(v)
    real(dp), intent(inout) :: v(:)
    real(dp) :: x
    integer :: a, b
    do a = 2, size(v)
      x = v(a);  b = a - 1
      do while (b >= 1)
        if (v(b) <= x) exit
        v(b+1) = v(b);  b = b - 1
      end do
      v(b+1) = x
    end do
  end subroutine ordina

  !> Mediana di un vettore gia' ordinato
  pure function mediana(v) result(md)
    real(dp), intent(in) :: v(:)
    real(dp) :: md
    integer :: nn
    nn = size(v)
    if (mod(nn,2) == 1) then
      md = v((nn+1)/2)
    else
      md = 0.5_dp*(v(nn/2) + v(nn/2+1))
    end if
  end function mediana

end program bench_speed_sp
