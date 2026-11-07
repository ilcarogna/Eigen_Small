# Test Suite for eigen_small

This directory contains comprehensive test programs for the `eigen_small` module, organized by precision:

- **`double/`**: Double precision (dp, `real64`) tests
- **`single/`**: Single precision (sp, `real32`) tests

Each test is available in both precision variants with `_dp` and `_sp` suffixes in the compiled binaries.

---

## Test Programs

### 1. `test_eigen` — Comprehensive Validation Suite

**Purpose**: Exhaustive validation of `eigen_small` accuracy and performance vs LAPACK.

**Test coverage (9 tests):**
1. **Symmetrize correctness**: Verifies `uplo='U'`, `'L'`, `'F'` produce identical eigenvalues
2. **2×2 known cases**: Identity, diagonal, dense symmetric, repeated eigenvalues, null matrix
3. **3×3 known cases**: Identity, diagonal (unsorted), dense, repeated eigenvalues, near-singular, scaled (large/small)
4. **Degenerate cases**: Null matrix, scalar multiples (`αI`), block-diagonal, nearly-coincident eigenvalues
5. **Negative eigenvalues**: Tests indefinite matrices (both 2×2 and 3×3)
6. **Interface consistency**: `eig_inplace` vs `eig_safe` vs `uplo` variants
7. **Statistical test**: 2M random SPD matrices, reports max error vs LAPACK
8. **Eigenvector accuracy** (`TEST 8`): residual `‖AU − UΛ‖` and orthonormality on
   - 2×2 at extreme scales (`1e-20`, `1e20`) and with a tiny off-diagonal term (cancellation in `λ₁ − a₁₁`)
   - 3×3 with two nearly coincident eigenvalues (relative gap down to `1e-8` in dp), cluster at the bottom and at the top; the isolated eigenvector must stay accurate against LAPACK
   - 200k random SPD matrices, max residual
9. **Benchmark**: Performance comparison (eigen_small vs LAPACK on 1M matrices; includes matrix generation time — see `bench_speed` for a kernel-only measurement)

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

### 5. `test_accuracy` — Accuracy vs LAPACK against a Quad-Precision Reference

**Purpose**: Check that `eigen_small` is as accurate as LAPACK, not only that it agrees with it. LAPACK is not used as the oracle here: both solvers are compared with a cyclic Jacobi solver in quad precision (`real128`) run on the same matrix.

**Coverage**: 2×2 and 3×3; four families (random SPD, random indefinite, two nearly coincident eigenvalues, mixed-sign eigenvalues on different scales); 100,000 matrices per family and size. Measures the maximum error of eigenvalues (relative to `max|λ|`), eigenvectors (angle × gap, relative to `ε‖A‖`; skipped when `gap < 100 ε‖A‖`) and residual `‖AU − UΛ‖`.

**Pass criterion**: for every family and quantity, error(`eigen_small`) ≤ 1.5 × error(LAPACK) + 1 and ≤ 20 (in units of `ε`). Exits with a non-zero status on failure.

**Determinism**: the random generator is a fixed-seed LCG, so results are reproducible and match the table in the main README.

**Run time**: about 20 s per precision.

**Build & run:**
```bash
bin/test_accuracy_dp   # Double precision
bin/test_accuracy_sp   # Single precision
```

---

### 6. `bench_speed` — Speed Benchmark, Matrix Generation Excluded

**Purpose**: Measure the real speed ratio between `eigen_small` and LAPACK (`dsyev`/`ssyev`) for 3×3 matrices, without counting the cost of producing the test matrices (unlike the benchmark in `test_eigen`).

**Methods timed**: `eig_inplace`, `eig_safe`, LAPACK `dsyev`/`ssyev` with `jobz='V'` (eigenvectors computed) and a preallocated workspace of optimal size (from a `lwork=-1` query).

**Design**:
- All 1M matrices are generated beforehand in `A(3,3,N)`; the timed loops contain no I/O.
- `eig_inplace` and LAPACK overwrite their input, so each iteration copies `A(:,:,k)` into a work matrix. The cost of that copy is timed by a baseline loop (same copy and checksum, no computation) and subtracted.
- Robust to CPU frequency drift, machine load and execution order: one untimed warm-up pass; the four methods are interleaved inside each repetition, in an order rotated at every repetition; each net time is `time − baseline` from the *same* repetition (paired measurement).
- 15 repetitions; reports **median** and minimum net time per matrix and the spread `(max − min)/median`. Speedups are ratios of medians. If any spread exceeds 10 % a warning is printed: the machine is noisy and the figures are unreliable (pin the process to one core, e.g. `taskset -c 2 bin/bench_speed_dp`, and check the CPU governor).
- Eigenvalue checksums prevent the compiler from removing the computation and verify that the methods agree.
- In the `sp` version timers, statistics and checksums are kept in double precision.

**Typical output** (gfortran `-O3 -march=native`, medians; varies with CPU and LAPACK build): a speedup of `eig_inplace` ≈ `eig_safe` over LAPACK of about 3.2× (dp) and 3.6× (sp), i.e. a few hundred ns per matrix against roughly one µs. `eig_inplace` and `eig_safe` must come out equal within noise: they share the same kernel.

**Build & run:**
```bash
bin/bench_speed_dp   # Double precision
bin/bench_speed_sp   # Single precision
```

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
- **Double precision**: Errors ≲ 10⁻¹³ (eigenvalues), ≲ 10⁻¹⁴ (orthonormality and residual `‖AU − UΛ‖`)
- **Single precision**: Errors ≲ 10⁻⁵ (eigenvalues and residual), ≲ 10⁻⁶ (orthonormality)

Benchmark speedup (`bench_speed`, matrix generation excluded): ~3.2× (dp) and ~3.6× (sp) faster than LAPACK for 3×3 matrices with gfortran and system LAPACK, ~7-8× with ifx and MKL (varies with CPU, compiler and LAPACK build; see the main README).
