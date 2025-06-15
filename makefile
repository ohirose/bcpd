BCPDSRC :=  register/*.c base/*.c
OUTDIR := .

# Detect platform
UNAME_S := $(shell uname -s)

# Platform-specific settings
ifeq ($(UNAME_S),Darwin)
    CC := clang
    CFLAGS += -Xpreprocessor -fopenmp
    LDFLAGS += -lomp -lopenblas

    # Detect Homebrew or MacPorts
    OMP_PREFIX := $(shell brew --prefix libomp 2>/dev/null)
    BLAS_PREFIX := $(shell brew --prefix openblas 2>/dev/null)

    ifeq ($(OMP_PREFIX),)
        OMP_PREFIX := /opt/local
    endif
    CFLAGS += -I$(OMP_PREFIX)/include -I$(BLAS_PREFIX)/include
    LDFLAGS += -L$(OMP_PREFIX)/lib -L$(BLAS_PREFIX)/lib

else ifeq ($(findstring MINGW,$(UNAME_S)),MINGW)
    # MinGW-w64 on Windows (e.g., MSYS2)
    CC := x86_64-w64-mingw32-gcc
    CFLAGS += -fopenmp -DMINGW64
    LDFLAGS += -lopenblas -static
    OUTDIR := ./win
else
    # Assume Linux
    CC := gcc
    CFLAGS += -fopenmp
    LDFLAGS += -lopenblas
endif

# Allow user to override CC, CFLAGS, LDFLAGS
CC ?= $(CC)

# Compilation
all:
	$(CC) -O3 $(CFLAGS) -DUSE_OPENMP $(BCPDSRC) $(LDFLAGS) -o $(OUTDIR)/bcpd
