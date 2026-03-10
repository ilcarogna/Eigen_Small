# eigen_small

**Fast analytical eigenvalue solver for small symmetric matrices (2×2, 3×3)**

[![Fortran](https://img.shields.io/badge/Fortran-2008-734f96?logo=fortran)](https://fortran-lang.org)
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

## Overview

`eigen_small` computes eigenvalues and eigenvectors of 2×2 and 3×3 symmetric matrices using **analytical formulas** (no iteration). It is an **almost plug-in replacement** for LAPACK's `dsyev`/`ssyev` for small matrices, offering ~15-30× speedup while maintaining comparable accuracy.

**Key features:**
- Analytical solution: O(1) complexity for 2×2, deterministic trigonometric method for 3×3
- Dual precision: generic interfaces for `real32` (sp) and `real64` (dp)
- Pure routines: no heap allocation or side effects
- Packed storage: optional `uplo` parameter for triangular input (`'U'`, `'L'`, or full `'F'`)
- Eigenvalues in ascending order: λ₁ ≤ λ₂ ≤ λ₃
- Orthonormal eigenvectors: column `k` = eigenvector for λₖ

## Algorithms

- **2×2**: Exact analytical formulas (quadratic characteristic equation)
- **3×3**: Smith (1961) / Kopp (2008) trigonometric method + Gram-Schmidt orthogonalization

### References

1. Smith, O.K. (1961). *Eigenvalues of a symmetric 3×3 matrix*. CACM **4**(4), 168. [DOI](https://doi.org/10.1145/355578.366316), [CACM](https://cacm.acm.org/research/eigenvalues-of-a-symmetric-3-x-3-matrix/)
2. Kopp, J. (2008). *Efficient numerical diagonalization of hermitian 3×3 matrices*. IJMPC **19**(03), 523. arXiv: [physics/0610206](https://arxiv.org/abs/physics/0610206)

## Usage

### Public Interfaces

```fortran
call eig_inplace(M, lambda [, uplo])  ! Overwrites M with eigenvectors
call eig_safe(M, lambda, U [, uplo])  ! Preserves M, stores eigenvectors in U
```

**Parameters:**
- `M`: Symmetric matrix (2×2 or 3×3)
- `lambda`: Eigenvalues (ascending order)
- `U`: Eigenvectors (column-major, for `eig_safe` only)
- `uplo`: `'U'` (upper), `'L'` (lower), or `'F'` (full, default)

### Example

```fortran
program example
  use iso_fortran_env, only: dp => real64
  use eigen_small
  implicit none
  
  real(dp) :: M(3,3), lambda(3), U(3,3)
  
  M = reshape([6.0_dp, 2.0_dp, 1.0_dp,  &
               2.0_dp, 5.0_dp, 3.0_dp,  &
               1.0_dp, 3.0_dp, 4.0_dp], [3,3])
  
  call eig_safe(M, lambda, U)
  
  print '(a, 3f10.5)', 'Eigenvalues: ', lambda
  print '(a, 3f10.5)', 'First eigenvector: ', U(:,1)
end program example
```

**Compile and run:**
```bash
gfortran -O3 eigen_small.f90 example.f90 -o example && ./example
```

## Performance

Benchmark speedup: **30× faster** than LAPACK for 3×3 matrices.

For 2×2 matrices, speedup reaches ~**30×**.

## Testing

Test suite in `tests/double/` (dp) and `tests/single/` (sp) validates against LAPACK:
```bash
make test    # Run all tests
make double  # Double precision only
make single  # Single precision only
```

Tests verify eigenvalue accuracy, orthonormality (`U^T U = I`), and reconstruction (`M = U Λ U^T`).

## Limitations

- **Fixed size**: only 2×2 and 3×3 (use LAPACK for larger matrices)
- **Symmetry not checked**: undefined behavior for non-symmetric input
- **Accuracy**: tested for `cond(M) < 10¹²` (dp); extremely ill-conditioned matrices may yield less accurate eigenvectors

## License

MIT License. See `LICENSE` file for details.

## Acknowledgments

This implementation was developed independently from scratch, following the algorithms in:

- **Smith, O.K. (1961)**: *Eigenvalues of a symmetric 3×3 matrix*. Analytical formulas for eigendecomposition ([CACM](https://cacm.acm.org/research/eigenvalues-of-a-symmetric-3-x-3-matrix/), [ACM DOI](https://dl.acm.org/doi/10.1145/355578.366316)).
- **Kopp, J. (2008)**: *Efficient numerical diagonalization of hermitian 3×3 matrices*. Trigonometric method for 3×3 eigendecomposition ([arXiv:physics/0610206](https://arxiv.org/abs/physics/0610206)).

Implementation: Modern Fortran 2008 with dual precision (sp/dp) support.

**Note**: Kopp's original C/C++/Fortran implementations (more comprehensive, includes complex/hermitian matrices) are available under LGPL at: https://www.mpi-hd.mpg.de/personalhomes/globes/3x3/

---

**Author**: Giancarlo de Gasperis  
**Version**: 1.0 (March 2026)
