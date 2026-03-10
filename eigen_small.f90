!==============================================================================
! MODULE eigen_small
!
! Calcolo analitico di autovalori e autovettori per matrici simmetriche
! semidefinite positive 2x2 e 3x3, in singola (sp) e doppia (dp) precisione.
!
! Algoritmi:
!   2x2 -- formule analitiche esatte (O(1), senza iterazione)
!   3x3 -- autovalori: metodo trigonometrico di Smith (1961) / Kopp (2008)
!          autovettori: prodotti vettoriali dei cofattori delle righe
!          raffinamento ortogonale via Gram-Schmidt modificato
!
! Argomento UPLO (character, optional, default 'F'):
!   'U' o 'u' -- solo il triangolo superiore di M e' significativo
!   'L' o 'l' -- solo il triangolo inferiore di M e' significativo
!   'F' o 'f' -- matrice completa (default)
!   In tutti i casi il kernel lavora su una copia W simmetrica completa
!   costruita da symmetrize; i dati originali di M non vengono mai letti
!   fuori dal triangolo selezionato.
!
! Interfacce pubbliche:
!   eig_inplace(M, lambda [,uplo])      M sovrascritto dagli autovettori (colonne)
!   eig_safe   (M, lambda, U [,uplo])  U pre-allocata dal chiamante, M intatta
!
! Convenzioni output:
!   - Autovalori in lambda(:) in ordine CRESCENTE
!   - Colonna k di U = autovettore normalizzato di lambda(k)
!
! Note di implementazione:
!   - Tutti i kernel interni sono PURE: nessun effetto collaterale,
!     nessuna allocazione heap; W e' sempre una variabile locale sullo stack.
!   - Le procedure pubbliche non sono pure per compatibilita' con
!     intent(inout) su M in eig_inplace.
!   - La logica uplo e' incapsulata interamente in symmetrize_dp/sp;
!     dispatch e i kernel non la conoscono.
!
! Riferimenti:
!   Smith, O.K. (1961). Eigenvalues of a symmetric 3x3 matrix.
!     Communications of the ACM, 4(4), 168.
!   Kopp, J. (2008). Efficient numerical diagonalization of hermitian 3x3
!     matrices. Int. J. Mod. Phys. C, 19(03), 523-548. arXiv:physics/0610206
!==============================================================================
module eigen_small
  use iso_fortran_env, only: sp => real32, dp => real64
  implicit none
  private

  ! Pi greco in dp/sp; usato nel metodo trigonometrico 3x3
  real(dp), parameter :: PI_DP = acos(-1.0_dp)
  real(sp), parameter :: PI_SP = acos(-1.0_sp)

  !----------------------------------------------------------------------------
  ! Interfacce pubbliche
  !----------------------------------------------------------------------------

  !> Diagonalizza M in-place: M viene sovrascritta dalla matrice U degli
  !  autovettori (per colonne). uplo indica quale triangolo di M leggere.
  interface eig_inplace
    module procedure :: eig_ip_sp, eig_ip_dp
  end interface eig_inplace

  !> Diagonalizza M senza modificarla; U deve essere pre-allocata dal chiamante
  !  (array statico, automatico o gia' allocato). uplo indica quale triangolo leggere.
  interface eig_safe
    module procedure :: eig_safe_sp, eig_safe_dp
  end interface eig_safe

  public :: eig_inplace, eig_safe

contains

  !=============================================================================
  ! SEZIONE 1 – SIMMETRIZZAZIONE
  !
  ! symmetrize costruisce una copia W simmetrica completa di M leggendo
  ! solo il triangolo indicato da uplo.
  ! E' l'UNICA procedura che conosce uplo; tutto il resto del modulo
  ! lavora su W completa e non deve preoccuparsi di quale triangolo usare.
  !
  ! Caso 'U': W(i,j) = M(i,j) per i<=j;  W(j,i) = M(i,j) per i>j
  ! Caso 'L': W(i,j) = M(i,j) per i>=j;  W(j,i) = M(i,j) per i<j
  ! Caso 'F': W = M  (copia diretta, nessuna riflessione)
  !=============================================================================

  pure subroutine symmetrize_dp(M, W, uplo)
    real(dp), intent(in)  :: M(:,:)
    real(dp), intent(out) :: W(size(M,1), size(M,1))
    character, intent(in) :: uplo
    integer :: i, j, n

    n = size(M, 1)

    select case (uplo)

    case ('U', 'u')
      ! Leggi solo il triangolo superiore (j >= i) e rifletti sotto la diagonale
      do j = 1, n
        do i = 1, j
          W(i,j) = M(i,j)
          W(j,i) = M(i,j)
        end do
      end do

    case ('L', 'l')
      ! Leggi solo il triangolo inferiore (i >= j) e rifletti sopra la diagonale
      do j = 1, n
        do i = j, n
          W(i,j) = M(i,j)
          W(j,i) = M(i,j)
        end do
      end do

    case default
      ! 'F', 'f' o qualsiasi altro: matrice gia' completa, copia diretta
      W = M

    end select
  end subroutine symmetrize_dp

  pure subroutine symmetrize_sp(M, W, uplo)
    real(sp), intent(in)  :: M(:,:)
    real(sp), intent(out) :: W(size(M,1), size(M,1))
    character, intent(in) :: uplo
    integer :: i, j, n

    n = size(M, 1)

    select case (uplo)
    case ('U', 'u')
      do j = 1, n
        do i = 1, j
          W(i,j) = M(i,j)
          W(j,i) = M(i,j)
        end do
      end do
    case ('L', 'l')
      do j = 1, n
        do i = j, n
          W(i,j) = M(i,j)
          W(j,i) = M(i,j)
        end do
      end do
    case default
      W = M
    end select
  end subroutine symmetrize_sp

  !=============================================================================
  ! SEZIONE 2 – UTILITA' INTERNE
  !=============================================================================

  !> Ordinamento crescente di 3 scalari (bubble sort ottimale per N=3)
  pure subroutine sort3_dp(v)
    real(dp), intent(inout) :: v(3)
    real(dp) :: t
    if (v(1) > v(2)) then; t=v(1); v(1)=v(2); v(2)=t; end if
    if (v(2) > v(3)) then; t=v(2); v(2)=v(3); v(3)=t; end if
    if (v(1) > v(2)) then; t=v(1); v(1)=v(2); v(2)=t; end if
  end subroutine sort3_dp

  pure subroutine sort3_sp(v)
    real(sp), intent(inout) :: v(3)
    real(sp) :: t
    if (v(1) > v(2)) then; t=v(1); v(1)=v(2); v(2)=t; end if
    if (v(2) > v(3)) then; t=v(2); v(2)=v(3); v(3)=t; end if
    if (v(1) > v(2)) then; t=v(1); v(1)=v(2); v(2)=t; end if
  end subroutine sort3_sp

  !> Determinante 3x3 (sviluppo di Sarrus sulla prima riga)
  pure function det3_dp(B) result(d)
    real(dp), intent(in) :: B(3,3)
    real(dp) :: d
    d =   B(1,1)*(B(2,2)*B(3,3) - B(2,3)*B(3,2)) &
        - B(1,2)*(B(2,1)*B(3,3) - B(2,3)*B(3,1)) &
        + B(1,3)*(B(2,1)*B(3,2) - B(2,2)*B(3,1))
  end function det3_dp

  pure function det3_sp(B) result(d)
    real(sp), intent(in) :: B(3,3)
    real(sp) :: d
    d =   B(1,1)*(B(2,2)*B(3,3) - B(2,3)*B(3,2)) &
        - B(1,2)*(B(2,1)*B(3,3) - B(2,3)*B(3,1)) &
        + B(1,3)*(B(2,1)*B(3,2) - B(2,2)*B(3,1))
  end function det3_sp

  !> Prodotto vettoriale a x b in R^3
  pure function cross_dp(a, b) result(c)
    real(dp), intent(in) :: a(3), b(3)
    real(dp) :: c(3)
    c(1) = a(2)*b(3) - a(3)*b(2)
    c(2) = a(3)*b(1) - a(1)*b(3)
    c(3) = a(1)*b(2) - a(2)*b(1)
  end function cross_dp

  pure function cross_sp(a, b) result(c)
    real(sp), intent(in) :: a(3), b(3)
    real(sp) :: c(3)
    c(1) = a(2)*b(3) - a(3)*b(2)
    c(2) = a(3)*b(1) - a(1)*b(3)
    c(3) = a(1)*b(2) - a(2)*b(1)
  end function cross_sp

  !> Dato un vettore unitario v, costruisce u1, u2 tali che {v, u1, u2}
  !  sia una terna ortonormale.
  !  Algoritmo: Gram-Schmidt sul vettore canonico meno allineato a v,
  !  poi prodotto vettoriale per garantire ortogonalita' esatta.
  pure subroutine two_ortho_dp(v, u1, u2)
    real(dp), intent(in)  :: v(3)
    real(dp), intent(out) :: u1(3), u2(3)
    real(dp) :: nrm

    ! Sceglie il vettore canonico meno parallelo a v (minimizza cancellazione)
    if (abs(v(1)) <= abs(v(2)) .and. abs(v(1)) <= abs(v(3))) then
      u1 = [1.0_dp, 0.0_dp, 0.0_dp]
    else if (abs(v(2)) <= abs(v(3))) then
      u1 = [0.0_dp, 1.0_dp, 0.0_dp]
    else
      u1 = [0.0_dp, 0.0_dp, 1.0_dp]
    end if
    u1  = u1 - dot_product(u1, v) * v   ! proiezione ortogonale a v
    nrm = norm2(u1)
    u1  = u1 / nrm
    u2  = cross_dp(v, u1)               ! ortogonalita' esatta per costruzione
  end subroutine two_ortho_dp

  pure subroutine two_ortho_sp(v, u1, u2)
    real(sp), intent(in)  :: v(3)
    real(sp), intent(out) :: u1(3), u2(3)
    real(sp) :: nrm

    if (abs(v(1)) <= abs(v(2)) .and. abs(v(1)) <= abs(v(3))) then
      u1 = [1.0_sp, 0.0_sp, 0.0_sp]
    else if (abs(v(2)) <= abs(v(3))) then
      u1 = [0.0_sp, 1.0_sp, 0.0_sp]
    else
      u1 = [0.0_sp, 0.0_sp, 1.0_sp]
    end if
    u1  = u1 - dot_product(u1, v) * v
    nrm = norm2(u1)
    u1  = u1 / nrm
    u2  = cross_sp(v, u1)
  end subroutine two_ortho_sp

  !=============================================================================
  ! SEZIONE 3 – KERNEL 2x2
  !
  ! Per A = [[a, b], [b, d]] gli autovalori sono:
  !   lambda(1,2) = (a+d)/2 -/+ sqrt(((a-d)/2)^2 + b^2)
  !
  ! Autovettori (caso b != 0):
  !   v1 proporzionale a [b, lambda1 - a]
  !   v2 = [-( lambda1-a ), b]  (ortogonale a v1 per costruzione)
  !
  ! Caso diagonale (b ~ 0): vettori canonici ordinati per autovalore.
  !=============================================================================

  pure subroutine core2_dp(a11, a12, a22, lambda, U)
    real(dp), intent(in)  :: a11, a12, a22
    real(dp), intent(out) :: lambda(2), U(2,2)
    real(dp) :: tr2, disc, nrm, thr

    tr2  = 0.5_dp * (a11 + a22)
    ! hypot evita overflow/underflow nel calcolo di sqrt(x^2 + y^2)
    disc = hypot(0.5_dp*(a11 - a22), a12)

    lambda(1) = tr2 - disc
    lambda(2) = tr2 + disc

    ! Soglia relativa per il test di off-diagonale
    thr = epsilon(a12) * max(abs(a11), abs(a22), 1.0_dp)

    if (abs(a12) > thr) then
      U(1,1) = a12
      U(2,1) = lambda(1) - a11
      nrm    = hypot(U(1,1),U(2,1))
      U(:,1) = U(:,1) / nrm
      ! Secondo autovettore ortogonale al primo per costruzione
      U(1,2) = -U(2,1)
      U(2,2) =  U(1,1)
    else
      ! Matrice gia' diagonale: vettori canonici, eventualmente scambiati
      if (a11 <= a22) then
        U(:,1) = [1.0_dp, 0.0_dp]
        U(:,2) = [0.0_dp, 1.0_dp]
      else
        U(:,1) = [0.0_dp, 1.0_dp]
        U(:,2) = [1.0_dp, 0.0_dp]
      end if
    end if
  end subroutine core2_dp

  pure subroutine core2_sp(a11, a12, a22, lambda, U)
    real(sp), intent(in)  :: a11, a12, a22
    real(sp), intent(out) :: lambda(2), U(2,2)
    real(sp) :: tr2, disc, nrm, thr

    tr2  = 0.5_sp * (a11 + a22)
    ! hypot evita overflow/underflow nel calcolo di sqrt(x^2 + y^2)
    disc = hypot(0.5_sp*(a11 - a22), a12)

    lambda(1) = tr2 - disc
    lambda(2) = tr2 + disc

    ! Soglia relativa per il test di off-diagonale
    thr = epsilon(a12) * max(abs(a11), abs(a22), 1.0_sp)

    if (abs(a12) > thr) then
      U(1,1) = a12
      U(2,1) = lambda(1) - a11
      nrm    = hypot(U(1,1),U(2,1))
      U(:,1) = U(:,1) / nrm
      ! Secondo autovettore ortogonale al primo per costruzione
      U(1,2) = -U(2,1)
      U(2,2) =  U(1,1)
    else
      ! Matrice gia' diagonale: vettori canonici, eventualmente scambiati
      if (a11 <= a22) then
        U(:,1) = [1.0_sp, 0.0_sp]
        U(:,2) = [0.0_sp, 1.0_sp]
      else
        U(:,1) = [0.0_sp, 1.0_sp]
        U(:,2) = [1.0_sp, 0.0_sp]
      end if
    end if
  end subroutine core2_sp

  !=============================================================================
  ! SEZIONE 4 – KERNEL 3x3
  !
  ! AUTOVALORI – metodo trigonometrico (Smith 1961, Kopp 2008):
  !   1) Test di diagonalita': off-diag^2 << diag^2 -> restituisce diag(A)
  !   2) q = trace(A)/3
  !   3) p = sqrt(||A-q*I||^2_F / 6)
  !   4) B = (A - q*I)/p  -> autovalori di B in [-2,2] (condizionamento ottimo)
  !   5) r = det(B)/2, clampato in [-1,1]  (protegge acos)
  !   6) phi = acos(r)/3
  !   7) lambda_k = q + 2*p*cos(phi + 2*pi*k/3), k=0,1,2
  !
  ! AUTOVETTORI – prodotti vettoriali dei cofattori:
  !   Le righe di (A - lambda*I) sono ortogonali all'autovettore cercato.
  !   Il prodotto vettoriale di due righe lin. indip. e' proporzionale
  !   all'autovettore. Si sceglie la coppia con norma massima (piu' stabile).
  !
  ! GESTIONE AUTOVALORI RIPETUTI (quattro casi):
  !   1) lambda1 = lambda2 = lambda3  ->  U = I
  !   2) lambda1 = lambda2 < lambda3  ->  v3 ben definito, two_ortho per v1,v2
  !   3) lambda1 < lambda2 = lambda3  ->  v1 ben definito, two_ortho per v2,v3
  !   4) tutti distinti               ->  tre eigvec3 + raffinamento GS
  !=============================================================================

  pure subroutine eigenvalues3_dp(A, lambda)
    real(dp), intent(in)  :: A(3,3)
    real(dp), intent(out) :: lambda(3)
    real(dp) :: p1, p2, q, p, r, phi, B(3,3), p_diag

    ! Somma quadratica degli elementi fuori diagonale
    p1 = A(1,2)**2 + A(1,3)**2 + A(2,3)**2

    ! Test di diagonalita' SCALE-INVARIANTE:
    !   confronta off-diag^2 con diag^2, senza costanti assolute.
    !   Il vecchio "+ 1.0" proteggeva il caso matrice nulla ma rendeva
    !   il test insensibile alla scala: per elementi O(1e-12) la soglia
    !   diventava epsilon_dp * 1.0 >> p1, classificando erroneamente
    !   come diagonale una matrice densa piccola.
    !   Correzione: il caso nulla e' coperto da "0 <= 0" senza addendi.
    p_diag = A(1,1)**2 + A(2,2)**2 + A(3,3)**2
    if (p1 <= epsilon(1.0_dp) * p_diag) then
      lambda = [A(1,1), A(2,2), A(3,3)]
      call sort3_dp(lambda)
      return
    end if

    q  = (A(1,1) + A(2,2) + A(3,3)) / 3.0_dp
    p2 = (A(1,1)-q)**2 + (A(2,2)-q)**2 + (A(3,3)-q)**2 + 2.0_dp*p1
    p  = sqrt(p2 / 6.0_dp)

    ! Matrice normalizzata B = (A - q*I)/p: autovalori in [-2, 2]
    B(1,1) = (A(1,1)-q)/p;  B(2,2) = (A(2,2)-q)/p;  B(3,3) = (A(3,3)-q)/p
    B(1,2) = A(1,2)/p;      B(2,1) = B(1,2)
    B(1,3) = A(1,3)/p;      B(3,1) = B(1,3)
    B(2,3) = A(2,3)/p;      B(3,2) = B(2,3)

    ! r teoricamente in [-1,1]; clamping protegge acos da errori floating-point
    r   = max(-1.0_dp, min(1.0_dp, 0.5_dp * det3_dp(B)))
    phi = acos(r) / 3.0_dp

    ! Tutti e tre calcolati con coseno diretto: evita la sottrazione di Vieta
    ! (stessa logica di eigenvalues3_sp, meno critica in dp ma piu' corretta).
    lambda(1) = q + 2.0_dp*p * cos(phi)
    lambda(2) = q + 2.0_dp*p * cos(phi + 2.0_dp*PI_DP/3.0_dp)
    lambda(3) = q + 2.0_dp*p * cos(phi + 4.0_dp*PI_DP/3.0_dp)
    call sort3_dp(lambda)
  end subroutine eigenvalues3_dp

  pure subroutine eigenvalues3_sp(A, lambda)
    real(sp), intent(in)  :: A(3,3)
    real(sp), intent(out) :: lambda(3)
    real(sp) :: p1, p2, q, p, r, phi, B(3,3), p_diag

    p1 = A(1,2)**2 + A(1,3)**2 + A(2,3)**2

    ! Test di diagonalita' scale-invariante (vedi commento in eigenvalues3_dp)
    p_diag = A(1,1)**2 + A(2,2)**2 + A(3,3)**2
    if (p1 <= epsilon(1.0_sp) * p_diag) then
      lambda = [A(1,1), A(2,2), A(3,3)]
      call sort3_sp(lambda)
      return
    end if

    q  = (A(1,1) + A(2,2) + A(3,3)) / 3.0_sp
    p2 = (A(1,1)-q)**2 + (A(2,2)-q)**2 + (A(3,3)-q)**2 + 2.0_sp*p1
    p  = sqrt(p2 / 6.0_sp)

    B(1,1) = (A(1,1)-q)/p;  B(2,2) = (A(2,2)-q)/p;  B(3,3) = (A(3,3)-q)/p
    B(1,2) = A(1,2)/p;      B(2,1) = B(1,2)
    B(1,3) = A(1,3)/p;      B(3,1) = B(1,3)
    B(2,3) = A(2,3)/p;      B(3,2) = B(2,3)

    r   = max(-1.0_sp, min(1.0_sp, 0.5_sp * det3_sp(B)))
    phi = acos(r) / 3.0_sp

    ! Tutti e tre calcolati con coseno diretto: evita la sottrazione di Vieta
    !   lambda(2) = 3q - lambda(1) - lambda(3)
    ! che perde cifre in singola precisione quando lambda(1) e lambda(3)
    ! sono entrambi vicini a q (cancellazione catastrofica).
    ! I tre angoli 0, 2pi/3, 4pi/3 coprono uniformemente il cerchio.
    lambda(1) = q + 2.0_sp*p * cos(phi)
    lambda(2) = q + 2.0_sp*p * cos(phi + 2.0_sp*PI_SP/3.0_sp)
    lambda(3) = q + 2.0_sp*p * cos(phi + 4.0_sp*PI_SP/3.0_sp)
    call sort3_sp(lambda)
  end subroutine eigenvalues3_sp

  !> Autovettore di A per l'autovalore lam.
  !  Sceglie tra le tre coppie di righe di (A - lam*I) quella il cui
  !  prodotto vettoriale ha norma massima (massima stabilita' numerica).
  pure subroutine eigvec3_dp(A, lam, v)
    real(dp), intent(in)  :: A(3,3), lam
    real(dp), intent(out) :: v(3)
    real(dp) :: r1(3), r2(3), r3(3), c12(3), c13(3), c23(3)
    real(dp) :: n12, n13, n23, nrm

    r1 = [A(1,1)-lam, A(1,2),     A(1,3)    ]
    r2 = [A(2,1),     A(2,2)-lam, A(2,3)    ]
    r3 = [A(3,1),     A(3,2),     A(3,3)-lam]

    c12 = cross_dp(r1, r2);  n12 = dot_product(c12, c12)
    c13 = cross_dp(r1, r3);  n13 = dot_product(c13, c13)
    c23 = cross_dp(r2, r3);  n23 = dot_product(c23, c23)

    if (n12 >= n13 .and. n12 >= n23) then
      v = c12;  nrm = sqrt(n12)
    else if (n13 >= n23) then
      v = c13;  nrm = sqrt(n13)
    else
      v = c23;  nrm = sqrt(n23)
    end if

    if (nrm > epsilon(nrm)) then
      v = v / nrm
    else
      v = [1.0_dp, 0.0_dp, 0.0_dp]   ! fallback: autovalore triplo degenere
    end if
  end subroutine eigvec3_dp

  pure subroutine eigvec3_sp(A, lam, v)
    real(sp), intent(in)  :: A(3,3), lam
    real(sp), intent(out) :: v(3)
    real(sp) :: r1(3), r2(3), r3(3), c12(3), c13(3), c23(3)
    real(sp) :: n12, n13, n23, nrm

    r1 = [A(1,1)-lam, A(1,2),     A(1,3)    ]
    r2 = [A(2,1),     A(2,2)-lam, A(2,3)    ]
    r3 = [A(3,1),     A(3,2),     A(3,3)-lam]

    c12 = cross_sp(r1, r2);  n12 = dot_product(c12, c12)
    c13 = cross_sp(r1, r3);  n13 = dot_product(c13, c13)
    c23 = cross_sp(r2, r3);  n23 = dot_product(c23, c23)

    if (n12 >= n13 .and. n12 >= n23) then
      v = c12;  nrm = sqrt(n12)
    else if (n13 >= n23) then
      v = c13;  nrm = sqrt(n13)
    else
      v = c23;  nrm = sqrt(n23)
    end if

    if (nrm > epsilon(nrm)) then
      v = v / nrm
    else
      v = [1.0_sp, 0.0_sp, 0.0_sp]   ! fallback: autovalore triplo degenere
    end if
  end subroutine eigvec3_sp

  !> Raffinamento di ortogonalita' per autovalori tutti distinti.
  !  Gram-Schmidt parziale: rimuove la componente di col2 lungo col1,
  !  poi usa il prodotto vettoriale per col3 (ortogonalita' esatta di macchina).
  pure subroutine gs3_dp(U)
    real(dp), intent(inout) :: U(3,3)
    real(dp) :: nrm
    U(:,2) = U(:,2) - dot_product(U(:,2), U(:,1)) * U(:,1)
    nrm = norm2(U(:,2))
    if (nrm > epsilon(nrm)) U(:,2) = U(:,2) / nrm
    U(:,3) = cross_dp(U(:,1), U(:,2))
  end subroutine gs3_dp

  pure subroutine gs3_sp(U)
    real(sp), intent(inout) :: U(3,3)
    real(sp) :: nrm
    U(:,2) = U(:,2) - dot_product(U(:,2), U(:,1)) * U(:,1)
    nrm = norm2(U(:,2))
    if (nrm > epsilon(nrm)) U(:,2) = U(:,2) / nrm
    U(:,3) = cross_sp(U(:,1), U(:,2))
  end subroutine gs3_sp

  !> Kernel completo 3x3: autovalori + autovettori con gestione casi ripetuti.
  pure subroutine core3_dp(A, lambda, U)
    real(dp), intent(in)  :: A(3,3)
    real(dp), intent(out) :: lambda(3), U(3,3)
    real(dp) :: tol, scale, As(3,3)

    ! Normalizza A per il calcolo degli autovettori.
    ! Gli autovettori di A e di A/scale sono identici, ma eigvec3 usa un
    ! threshold ASSOLUTO (epsilon) sul modulo del cross-product.
    ! Per matrici con elementi molto piccoli (es. O(1e-12)) i cross-product
    ! hanno norma ~O(scale^2) che puo' scendere sotto epsilon, portando
    ! al fallback [1,0,0] per tutti gli autovettori (|U^TU-I| ~ 1).
    ! Con A normalizzata i cross-product sono O(1) e il threshold funziona.

    scale = maxval(abs(A))
    if (scale > 0.0_dp) then
      As = A / scale
    else
      ! Matrice nulla: qualsiasi base ortonormale e' corretta
      U = 0.0_dp;  U(1,1)=1.0_dp;  U(2,2)=1.0_dp;  U(3,3)=1.0_dp
      lambda = 0.0_dp
      return
    end if

    call eigenvalues3_dp(As, lambda)     ! ← As normalizzata
    lambda = lambda * scale              ! ← ri-scala dopo

    ! Tolleranza per il riconoscimento di autovalori coincidenti.
    ! Usa max(|lambda|, epsilon) invece di max(|lambda|, 1) per restare
    ! scale-invariante anche per matrici con autovalori molto piccoli.
    tol = epsilon(1.0_dp) * 128.0_dp * &
          max(abs(lambda(1)), abs(lambda(3)), epsilon(1.0_dp))

    if (abs(lambda(3) - lambda(1)) < tol) then
      ! Caso 1: lambda1 = lambda2 = lambda3 -> base canonica
      U = 0.0_dp
      U(1,1) = 1.0_dp;  U(2,2) = 1.0_dp;  U(3,3) = 1.0_dp

    else if (abs(lambda(2) - lambda(1)) < tol) then
      ! Caso 2: lambda1 = lambda2 < lambda3 -> v3 ben definito
      ! Passa As e lambda scalato: eigvec3 lavora su matrice O(1)
      call eigvec3_dp(As, lambda(3)/scale, U(:,3))
      call two_ortho_dp(U(:,3), U(:,1), U(:,2))

    else if (abs(lambda(3) - lambda(2)) < tol) then
      ! Caso 3: lambda1 < lambda2 = lambda3 -> v1 ben definito
      call eigvec3_dp(As, lambda(1)/scale, U(:,1))
      call two_ortho_dp(U(:,1), U(:,2), U(:,3))

    else
      ! Caso 4: tutti distinti -> tre autovettori + raffinamento GS
      call eigvec3_dp(As, lambda(1)/scale, U(:,1))
      call eigvec3_dp(As, lambda(2)/scale, U(:,2))
      call eigvec3_dp(As, lambda(3)/scale, U(:,3))
      call gs3_dp(U)
    end if
  end subroutine core3_dp

  pure subroutine core3_sp(A, lambda, U)
    real(sp), intent(in)  :: A(3,3)
    real(sp), intent(out) :: lambda(3), U(3,3)
    real(sp) :: tol, scale, As(3,3)

    ! Normalizza A per il calcolo degli autovettori (vedi commento in core3_dp)
    scale = maxval(abs(A))
    if (scale > 0.0_sp) then
      As = A / scale
    else
      U = 0.0_sp;  U(1,1)=1.0_sp;  U(2,2)=1.0_sp;  U(3,3)=1.0_sp
      lambda = 0.0_sp
      return
    end if

    call eigenvalues3_sp(As, lambda)     ! ← As normalizzata
    lambda = lambda * scale              ! ← ri-scala dopo

    ! Tolleranza per il riconoscimento di autovalori coincidenti.
    ! Usa max(|lambda|, epsilon) invece di max(|lambda|, 1) per restare
    ! scale-invariante anche per matrici con autovalori molto piccoli.
    tol = epsilon(1.0_sp) * 128.0_sp * &
          max(abs(lambda(1)), abs(lambda(3)), epsilon(1.0_sp))

    if (abs(lambda(3) - lambda(1)) < tol) then
      ! Caso 1: lambda1 = lambda2 = lambda3 -> base canonica
      U = 0.0_sp
      U(1,1) = 1.0_sp;  U(2,2) = 1.0_sp;  U(3,3) = 1.0_sp

    else if (abs(lambda(2) - lambda(1)) < tol) then
      ! Caso 2: lambda1 = lambda2 < lambda3 -> v3 ben definito
      call eigvec3_sp(As, lambda(3)/scale, U(:,3))
      call two_ortho_sp(U(:,3), U(:,1), U(:,2))

    else if (abs(lambda(3) - lambda(2)) < tol) then
      ! Caso 3: lambda1 < lambda2 = lambda3 -> v1 ben definito
      call eigvec3_sp(As, lambda(1)/scale, U(:,1))
      call two_ortho_sp(U(:,1), U(:,2), U(:,3))

    else
      ! Caso 4: tutti distinti -> tre autovettori + raffinamento GS
      call eigvec3_sp(As, lambda(1)/scale, U(:,1))
      call eigvec3_sp(As, lambda(2)/scale, U(:,2))
      call eigvec3_sp(As, lambda(3)/scale, U(:,3))
      call gs3_sp(U)
    end if
  end subroutine core3_sp

  !=============================================================================
  ! SEZIONE 5 – DISPATCHER INTERNO
  !
  ! Unico punto di ingresso ai kernel per le procedure pubbliche.
  ! Riceve W gia' simmetrica completa (prodotta da symmetrize) e
  ! chiama core2 o core3 in base a size(W,1).
  !=============================================================================

  pure subroutine dispatch_dp(W, lambda, U)
    real(dp), intent(in)  :: W(:,:)
    real(dp), intent(out) :: lambda(:), U(:,:)

    select case (size(W, 1))
    case (2)
      call core2_dp(W(1,1), W(1,2), W(2,2), lambda(1:2), U(1:2,1:2))
    case (3)
      call core3_dp(W(1:3,1:3), lambda(1:3), U(1:3,1:3))
    case default
      error stop "eigen_small (dp): supportate solo matrici 2x2 o 3x3"
    end select
  end subroutine dispatch_dp

  pure subroutine dispatch_sp(W, lambda, U)
    real(sp), intent(in)  :: W(:,:)
    real(sp), intent(out) :: lambda(:), U(:,:)

    select case (size(W, 1))
    case (2)
      call core2_sp(W(1,1), W(1,2), W(2,2), lambda(1:2), U(1:2,1:2))
    case (3)
      call core3_sp(W(1:3,1:3), lambda(1:3), U(1:3,1:3))
    case default
      error stop "eigen_small (sp): supportate solo matrici 2x2 o 3x3"
    end select
  end subroutine dispatch_sp

  !=============================================================================
  ! SEZIONE 6 – PROCEDURE PUBBLICHE
  !
  ! Schema comune:
  !   1) Leggi uplo (default 'F' se assente)
  !   2) Costruisci W simmetrica completa via symmetrize  [stack, nessun heap]
  !   3) dispatch -> kernel -> risultato
  !
  ! eig_inplace: scrive gli autovettori direttamente in M (M e' inout)
  ! eig_safe:    scrive gli autovettori in U, M e' intent(in) e non cambia
  !=============================================================================

  subroutine eig_ip_dp(M, lambda, uplo)
    real(dp),  intent(inout)        :: M(:,:)   ! [in] matrice; [out] autovettori
    real(dp),  intent(out)          :: lambda(:)
    character, intent(in), optional :: uplo     ! 'U', 'L', 'F' (default 'F')
    real(dp)  :: W(size(M,1), size(M,1))
    character :: uplo_

    uplo_ = 'F';  if (present(uplo)) uplo_ = uplo
    call symmetrize_dp(M, W, uplo_)    ! W = copia simmetrica completa di M
    call dispatch_dp(W, lambda, M)     ! autovettori scritti direttamente in M
  end subroutine eig_ip_dp

  subroutine eig_ip_sp(M, lambda, uplo)
    real(sp),  intent(inout)        :: M(:,:)
    real(sp),  intent(out)          :: lambda(:)
    character, intent(in), optional :: uplo
    real(sp)  :: W(size(M,1), size(M,1))
    character :: uplo_

    uplo_ = 'F';  if (present(uplo)) uplo_ = uplo
    call symmetrize_sp(M, W, uplo_)
    call dispatch_sp(W, lambda, M)
  end subroutine eig_ip_sp

  subroutine eig_safe_dp(M, lambda, U, uplo)
    real(dp),  intent(in)           :: M(:,:)   ! non viene mai modificata
    real(dp),  intent(out)          :: lambda(:)
    real(dp),  intent(out)          :: U(:,:)   ! pre-allocata dal chiamante
    character, intent(in), optional :: uplo
    real(dp)  :: W(size(M,1), size(M,1))
    character :: uplo_

    uplo_ = 'F';  if (present(uplo)) uplo_ = uplo
    call symmetrize_dp(M, W, uplo_)
    call dispatch_dp(W, lambda, U)
  end subroutine eig_safe_dp

  subroutine eig_safe_sp(M, lambda, U, uplo)
    real(sp),  intent(in)           :: M(:,:)
    real(sp),  intent(out)          :: lambda(:)
    real(sp),  intent(out)          :: U(:,:)
    character, intent(in), optional :: uplo
    real(sp)  :: W(size(M,1), size(M,1))
    character :: uplo_

    uplo_ = 'F';  if (present(uplo)) uplo_ = uplo
    call symmetrize_sp(M, W, uplo_)
    call dispatch_sp(W, lambda, U)
  end subroutine eig_safe_sp

end module eigen_small

