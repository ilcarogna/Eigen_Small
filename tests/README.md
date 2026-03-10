# Test Suite for eigen_small

This directory contains comprehensive test programs for the `eigen_small` module, organized by precision:

- **`double/`**: Double precision (dp, `real64`) tests
- **`single/`**: Single precision (sp, `real32`) tests

Each test is available in both precision variants with `_dp` and `_sp` suffixes in the compiled binaries.

---

## Test Programs

### 1. `test_eigen` — Comprehensive Validation Suite

**Purpose**: Exhaustive validation of `eigen_small` accuracy and performance vs LAPACK.

**Test coverage (8 tests):**
1. **Symmetrize correctness**: Verifies `uplo='U'`, `'L'`, `'F'` produce identical eigenvalues
2. **2×2 known cases**: Identity, diagonal, dense symmetric, repeated eigenvalues, null matrix
3. **3×3 known cases**: Identity, diagonal (unsorted), dense, repeated eigenvalues, near-singular, scaled (large/small)
4. **Degenerate cases**: Null matrix, scalar multiples (`αI`), block-diagonal, nearly-coincident eigenvalues
5. **Negative eigenvalues**: Tests indefinite matrices (both 2×2 and 3×3)
6. **Interface consistency**: `eig_inplace` vs `eig_safe` vs `uplo` variants
7. **Statistical test**: 2M random SPD matrices, reports max error vs LAPACK
8. **Benchmark**: Performance comparison (eigen_small vs LAPACK on 1M matrices)

**Validation strategy**: Uses LAPACK (`dsyev`/`ssyev`) as reference oracle, not hard-coded values. Checks:
- `|λ_eigen_small - λ_LAPACK| < ε`
- `max|U^T U - I| < ε` (orthonormality)
- `max|M - U Λ U^T| < ε` (reconstruction)

**Build & run:**
```bash
make double  # Compiles bin/test_eigen_dp
bin/test_eigen_dp

make single  # Compiles bin/test_eigen_sp
bin/test_eigen_sp
```

**Expected output**: `TOTALE: X OK, 0 FAIL` (all tests pass), plus benchmark timings.

---

### 2. `test_pseudoinverse` — Pseudo-Inverse via Eigendecomposition (`eig_safe`)

**Purpose**: Demonstrates using eigendecomposition to compute the Moore-Penrose pseudo-inverse and solve linear systems `Ax = b`.

**Algorithm**:
1. Compute `M = U Λ U^T` via `eig_safe`
2. Construct pseudo-inverse: `M^+ = U Λ^+ U^T` where `Λ^+_ii = 1/λ_i` if `λ_i > threshold`, else `0`
3. Solve: `x = M^+ b`

**Validation**: Compares `eigen_small` pseudo-inverse with LAPACK's eigendecomposition on a 3×3 test matrix. Checks residual `|Ax - b|`.

**Build & run:**
```bash
bin/test_pseudoinverse_dp   # Double precision
bin/test_pseudoinverse_sp   # Single precision
```

---

### 3. `test_pseudoinverse_inplace` — Pseudo-Inverse via `eig_inplace`

**Purpose**: Same as `test_pseudoinverse`, but uses `eig_inplace` (overwrites matrix with eigenvectors) instead of `eig_safe`.

**Use case**: When you don't need to preserve the original matrix, `eig_inplace` avoids the copy overhead.

**Build & run:**
```bash
bin/test_pseudoinverse_inplace_dp
bin/test_pseudoinverse_inplace_sp
```

---

### 4. `test_multiple_matrices` — Batch Processing of Linear Systems

**Purpose**: Tests solving **1000 linear systems** `A_i x = b_i` in batch mode using pseudo-inverse.

**Workflow**:
1. Generate 1000 random 3×3 SPD matrices and corresponding `b` vectors
2. Solve each system via eigendecomposition (both `eigen_small` and LAPACK)
3. Compute residuals: `r_i = |A_i x_i - b_i|`
4. Report statistics: max/mean residuals, max discrepancy between methods

**Validation**: Verifies that `eigen_small` and LAPACK produce consistent solutions across multiple systems.

**Build & run:**
```bash
bin/test_multiple_matrices_dp
bin/test_multiple_matrices_sp
```

**Note**: Adjust `N_MATRICES` parameter in source to test with more/fewer systems.

---

## Building All Tests

```bash
make double        # Build all dp tests
make single        # Build all sp tests
make all           # Build both dp and sp tests
make list-double   # List dp test executables
make list-single   # List sp test executables
```

Executables are placed in `bin/` with `_dp` (double) or `_sp` (single) suffixes.

---

## Compilation Requirements

- **Fortran compiler**: `gfortran ≥9.0` or Intel `ifx ≥2023`
- **LAPACK/BLAS**: Required for all tests (comparison with reference implementation)
  - Linux: `-llapack -lblas` (default)
  - Intel MKL: `-mkl` (when using `ifx`)

See main `Makefile` in project root for compiler options.

---

## Expected Test Results

All tests should **pass** (0 FAIL) on standard hardware with:
- **Double precision**: Errors ≲ 10⁻¹³ (eigenvalues), ≲ 10⁻¹⁴ (orthonormality)
- **Single precision**: Errors ≲ 10⁻⁵ (eigenvalues), ≲ 10⁻⁶ (orthonormality)

Benchmark speedup: ~15-30× faster than LAPACK for 2×2/3×3 matrices.
