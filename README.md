# eigen_small

**Fast analytical eigenvalue solver for small symmetric matrices (2×2, 3×3)**

[![Fortran](https://img.shields.io/badge/Fortran-2008-734f96?logo=fortran)](https://fortran-lang.org)
[![License: GPL v3+](https://img.shields.io/badge/License-GPLv3%2B-blue.svg)](LICENSE)

## Overview

`eigen_small` computes eigenvalues and eigenvectors of 2×2 and 3×3 real, symmetric matrices using **analytical formulas** (no iteration). It is an **almost plug-in replacement** for LAPACK's `dsyev`/`ssyev` for small matrices, offering a ~3-4× speedup (see [Performance](#performance)) with accuracy on par with LAPACK, including nearly singular matrices and clustered eigenvalues.

**Key features:**
- Analytical solution: O(1) complexity for 2×2, deterministic trigonometric method for 3×3
- Dual precision: generic interfaces for `real32` (sp) and `real64` (dp)
- Pure routines: no heap allocation or side effects
- Packed storage: optional `uplo` parameter for triangular input (`'U'`, `'L'`, or full `'F'`)
- Eigenvalues in ascending order: λ₁ ≤ λ₂ ≤ λ₃
- Orthonormal eigenvectors: column `k` = eigenvector for λₖ

## Algorithms

- **2×2**: Exact analytical formulas (quadratic characteristic equation)
- **3×3**: eigenvalues by the Smith (1961) / Kopp (2008) trigonometric method; the eigenvector of the most isolated eigenvalue comes from row cofactors, and the remaining pair is solved by a projected 2×2 problem in the orthogonal plane (deflation), with eigenvalues refined by Rayleigh quotients

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

Speedup over LAPACK `dsyev` for 3×3 double-precision matrices, measured with the bundled benchmark (`bin/test_eigen_dp`, 1M random SPD matrices, matrix generation time included in both columns; 20-core x86-64 machine):

| Toolchain | `eigen_small` | LAPACK `dsyev` | Speedup |
|---|---|---|---|
| gfortran `-O3 -march=native`, system LAPACK | ~320 ns | ~1130 ns | **~3.5×** |
| ifx `-O3 -xHost`, MKL (`-qmkl=sequential`) | ~329 ns | ~2290 ns | **~7×** |
| ifx `-O3 -xHost`, MKL (`-qmkl`, threaded) | ~323 ns | ~2650 ns | **~8×** |

`eigen_small` runs at the same speed with both compilers; the difference comes from `dsyev`, which carries a larger fixed overhead in MKL for tiny matrices. Kernel-only timings with gfortran (matrices precomputed) are ~74 ns vs ~265 ns for 2×2 and ~305 ns vs ~1050 ns for 3×3, i.e. also ~3.5×. The dedicated benchmark `bin/bench_speed_dp` / `bin/bench_speed_sp` excludes matrix generation altogether and is robust to CPU noise (interleaved methods in rotating order, paired baseline subtraction, medians over 15 repetitions; see [tests/README.md](tests/README.md)). With gfortran it gives, for 3×3 matrices with eigenvectors, a speedup over `dsyev` of **~3.2×** in dp and **~3.6×** in sp (`eig_inplace` and `eig_safe` share the same kernel and run at the same speed). Results vary with CPU, compiler and LAPACK build; run the benchmark on your own setup. For Intel builds use `make FC=ifx` after sourcing oneAPI's `setvars.sh`.

## Accuracy

Accuracy is compared with LAPACK (`dsyev`/`ssyev`) against a quad-precision Jacobi reference by the test programs `bin/test_accuracy_dp` and `bin/test_accuracy_sp` (about 20 s each; the random generator is deterministic, so the figures below are exactly reproducible). Each run uses 100,000 matrices for each of four families and for both sizes:

- random SPD matrices;
- random indefinite matrices;
- two nearly coincident eigenvalues (relative gap from `1e-1` down to `1e-9` in dp, `1e-5` in sp);
- eigenvalues of mixed sign spanning several orders of magnitude (up to `1e6` in dp, `1e3` in sp).

Worst case over the four families, in units of machine epsilon times `‖A‖` (eigenvectors: units of `ε‖A‖/gap`, the unavoidable error of any backward-stable solver):

| | | eigenvalues | | eigenvectors | | residual `‖AU − UΛ‖` | |
|---|---|---|---|---|---|---|---|
| | | `eigen_small` | LAPACK | `eigen_small` | LAPACK | `eigen_small` | LAPACK |
| dp | 2×2 | 1.3 | 2.1 | 1.3 | 1.6 | 1.7 | 2.4 |
| dp | 3×3 | 3.9 | 7.7 | 2.7 | 5.2 | 3.5 | 6.8 |
| sp | 2×2 | 1.4 | 1.9 | 1.2 | 1.6 | 1.7 | 2.2 |
| sp | 3×3 | 3.8 | 7.8 | 2.5 | 5.0 | 3.8 | 7.1 |

On these classes of matrices `eigen_small` is as accurate as LAPACK or better: about the same for 2×2 and up to about 2× smaller errors for 3×3. All errors are a few units of `ε‖A‖`, as expected for a backward-stable computation. This also holds for clustered eigenvalues, where naive analytical formulas lose accuracy; the deflation scheme described above is what avoids that. The test passes when, for every family and quantity, the error of `eigen_small` is at most `1.5 ×` LAPACK's `+ 1` and at most 20 units, so it also guards against accuracy regressions.

The figures are empirical, measured on the families listed, not a proof of a bound. Not covered: errors measured *relative to each eigenvalue* (very small eigenvalues next to large ones are accurate to `ε‖A‖`, like LAPACK's), matrices larger than 3×3, and non-symmetric input.

## Testing

Test suite in `tests/double/` (dp) and `tests/single/` (sp) validates against LAPACK:
```bash
make help    # Print a help message
make double  # Double precision only
make single  # Single precision only
```

Speed benchmark (matrix generation excluded): `bin/bench_speed_dp` and `bin/bench_speed_sp`.

Tests verify eigenvalue accuracy, orthonormality (`U^T U = I`), the residual `‖M U − U Λ‖`, and reconstruction (`M = U Λ U^T`), including extreme scales and nearly coincident eigenvalues, and compare accuracy with LAPACK against a quad-precision reference.

## Limitations

- **Fixed size**: only 2×2 and 3×3 (use LAPACK for larger matrices)
- **Symmetry not checked**: undefined behavior for non-symmetric input if `uplo`=`'F'`
- **Accuracy**: eigenvectors of eigenvalues separated by a relative gap `g` carry the unavoidable error ~`ε/g`, as with any backward-stable solver; against a quad-precision reference the errors are at or below LAPACK's for all gaps tested (down to `1e-9·λ_max` in dp)

## License

Copyright (C) 2026 Giancarlo de Gasperis

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version (SPDX: `GPL-3.0-or-later`). See the [LICENSE](LICENSE) file for details.

Versions up to and including commit `6db02fc` were released under the MIT License; copies obtained from those versions remain available under the MIT terms.

## Acknowledgments

This implementation was developed independently from scratch, following the algorithms described in:

- **Smith, O.K. (1961)**: *Eigenvalues of a symmetric 3×3 matrix*. Analytical formulas for eigendecomposition ([CACM](https://cacm.acm.org/research/eigenvalues-of-a-symmetric-3-x-3-matrix/), [ACM DOI](https://dl.acm.org/doi/10.1145/355578.366316)).
- **Kopp, J. (2008)**: *Efficient numerical diagonalization of hermitian 3×3 matrices*. Trigonometric method for 3×3 eigendecomposition ([arXiv:physics/0610206](https://arxiv.org/abs/physics/0610206)).

Implementation: Modern Fortran 2008 with dual precision (sp/dp) support.

**Note**: Kopp's original C/C++/Fortran implementations (more comprehensive, includes complex/hermitian matrices) are available under LGPL at: https://www.mpi-hd.mpg.de/personalhomes/globes/3x3/

---

**Author**: Giancarlo de Gasperis  
**Version**: 1.0 (March 2026)
