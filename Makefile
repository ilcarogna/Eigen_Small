# Makefile per progetto Eigen_Small
# Build di test in doppia/singola precisione da tests/double e tests/single

# Compiler selection: gfortran (default) o ifx, es. "make FC=ifx".
# FC e' una variabile predefinita di make (f77): si usa gfortran solo se
# l'utente non l'ha impostata da riga di comando o da ambiente.
ifeq ($(origin FC),default)
    FC := gfortran
endif
COMPILER := $(shell basename $(FC))

# Compiler flags
ifeq ($(COMPILER),gfortran)
    FFLAGS := -O3 -march=native -Wall -Wextra
    MODFLAG := -J
    LAPACK_LIBS := -llapack -lblas
else ifeq ($(COMPILER),ifx)
    FFLAGS := -O3 -xHost -warn all
    MODFLAG := -module
    # -mkl e' deprecato nelle versioni recenti di oneAPI; sequential evita il
    # threading di MKL, inutile (e dannoso) su matrici 2x2/3x3
    LAPACK_LIBS := -qmkl=sequential
else
    FFLAGS := -O3
    MODFLAG := -J
    LAPACK_LIBS := -llapack -lblas
endif

# Add debug flags if DEBUG=1
ifdef DEBUG
    FFLAGS := $(filter-out -O%,$(FFLAGS)) -O0 -g -fbacktrace -fbounds-check
endif

# Directories
SRCDIR := .
TESTDIR := ./tests
TESTDIR_DP := $(TESTDIR)/double
TESTDIR_SP := $(TESTDIR)/single
BINDIR := ./bin
OBJDIR := ./obj
MODDIR := $(OBJDIR)

# Core module
EIGEN_SRC := $(SRCDIR)/eigen_small.f90
EIGEN_OBJ := $(OBJDIR)/eigen_small.o
EIGEN_MOD := $(MODDIR)/eigen_small.mod

# Test sources
TESTS_DP_SRC := $(wildcard $(TESTDIR_DP)/*.f90)
TESTS_SP_SRC := $(wildcard $(TESTDIR_SP)/*_sp.f90)

# Objects
TESTS_DP_OBJ := $(patsubst $(TESTDIR_DP)/%.f90,$(OBJDIR)/%_dp.o,$(TESTS_DP_SRC))
TESTS_SP_OBJ := $(patsubst $(TESTDIR_SP)/%.f90,$(OBJDIR)/%.o,$(TESTS_SP_SRC))

# Executables in bin/ with _dp and _sp suffixes
TESTS_DP_EXE := $(patsubst $(TESTDIR_DP)/%.f90,$(BINDIR)/%_dp,$(TESTS_DP_SRC))
TESTS_SP_EXE := $(patsubst $(TESTDIR_SP)/%.f90,$(BINDIR)/%,$(TESTS_SP_SRC))

# Phony targets
.PHONY: all double single clean help dirs info list list-double list-single

# Default target
all: double single

# Build only double-precision tests
double: dirs $(TESTS_DP_EXE)

# Build only single-precision tests
single: dirs $(TESTS_SP_EXE)

# Create directories if they don't exist
dirs:
	@mkdir -p $(BINDIR) $(OBJDIR) $(MODDIR)

# Compile eigen_small module
$(EIGEN_OBJ): $(EIGEN_SRC) | dirs
	@echo "Compiling eigen_small module..."
	$(FC) $(FFLAGS) $(MODFLAG) $(MODDIR) -c $< -o $@

# Compile dp tests from tests/double
$(OBJDIR)/%_dp.o: $(TESTDIR_DP)/%.f90 $(EIGEN_OBJ) | dirs
	@echo "Compiling DP test $<..."
	$(FC) $(FFLAGS) $(MODFLAG) $(MODDIR) -I $(MODDIR) -c $< -o $@

# Compile sp tests from tests/single
$(OBJDIR)/%_sp.o: $(TESTDIR_SP)/%_sp.f90 $(EIGEN_OBJ) | dirs
	@echo "Compiling SP test $<..."
	$(FC) $(FFLAGS) $(MODFLAG) $(MODDIR) -I $(MODDIR) -c $< -o $@

# Link dp executables
$(BINDIR)/%_dp: $(EIGEN_OBJ) $(OBJDIR)/%_dp.o
	@echo "Linking $@..."
	$(FC) $(FFLAGS) -o $@ $^ $(LAPACK_LIBS)
	@echo "Build successful: $@"

# Link sp executables
$(BINDIR)/%_sp: $(EIGEN_OBJ) $(OBJDIR)/%_sp.o
	@echo "Linking $@..."
	$(FC) $(FFLAGS) -o $@ $^ $(LAPACK_LIBS)
	@echo "Build successful: $@"

# List detected tests
list:
	@echo "DP sources: $(TESTS_DP_SRC)"
	@echo "SP sources: $(TESTS_SP_SRC)"
	@echo "DP executables: $(TESTS_DP_EXE)"
	@echo "SP executables: $(TESTS_SP_EXE)"

# List only double-precision tests
list-double:
	@echo "DP sources: $(TESTS_DP_SRC)"
	@echo "DP executables: $(TESTS_DP_EXE)"

# List only single-precision tests
list-single:
	@echo "SP sources: $(TESTS_SP_SRC)"
	@echo "SP executables: $(TESTS_SP_EXE)"

# Clean build artifacts
clean:
	@echo "Cleaning build artifacts..."
	rm -rf $(OBJDIR) $(BINDIR)
	@echo "Clean complete"

# Help
help:
	@echo "Usage: make [options] [target]"
	@echo ""
	@echo "Targets:"
	@echo "  all          - Build all tests (default)"
	@echo "  double       - Build only tests in tests/double"
	@echo "  single       - Build only tests in tests/single"
	@echo "  list         - Show discovered sources/executables"
	@echo "  list-double  - Show only double-precision sources/executables"
	@echo "  list-single  - Show only single-precision sources/executables"
	@echo "  clean        - Remove build artifacts"
	@echo "  help         - Show this help message"
	@echo ""
	@echo "Structure:"
	@echo "  tests/double/*.f90    -> bin/*_dp"
	@echo "  tests/single/*_sp.f90 -> bin/*_sp"
	@echo ""
	@echo "Variables:"
	@echo "  FC        - Fortran compiler (default: gfortran)"
	@echo "              Example: make FC=ifx  (after sourcing oneAPI setvars.sh)"
	@echo "  DEBUG=1   - Enable debug flags"

# Display compiler info
info:
	@echo "Compiler: $(FC)"
	@echo "Flags: $(FFLAGS)"
	@echo "LAPACK: $(LAPACK_LIBS)"
