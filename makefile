CC=gcc
BCPDSRC=  register/*.c base/*.c
OMP_PORT= -Xpreprocessor -fopenmp -I/opt/local/include/libomp /opt/local/lib/libomp/libomp.dylib
OMP_BREW_ITL= -Xpreprocessor -fopenmp -I/usr/local/include/ /usr/local/lib/libomp.dylib
OMP_BREW_ARM= -Xpreprocessor -fopenmp -I/opt/homebrew/include/ /opt/homebrew/lib/libomp.dylib
DEBUG= 

all:
ifeq ($(OPT),-DUSE_OPENMP)
  ifeq ($(ENV),LINUX)
	$(CC) -O3 -fopenmp $(OPT) $(DEBUG) $(BCPDSRC) -o bcpd -lm -llapack
  else ifeq ($(ENV),MINGW32)
	$(CC) -O3 -fopenmp $(OPT) $(DEBUG) $(BCPDSRC) win/*.dll -o win/bcpd -lm -DMINGW32
  else ifeq ($(ENV),HOMEBREW_INTEL)
	clang -O3 $(OPT) $(OMP_BREW_ITL) $(DEBUG) $(BCPDSRC) -o bcpd -lm -llapack -Wuninitialized
  else ifeq ($(ENV),HOMEBREW)
	clang -O3 $(OPT) $(OMP_BREW_ARM) $(DEBUG) $(BCPDSRC) -o bcpd -lm -llapack -Wuninitialized
  else ifeq ($(ENV),MACPORTS)
	clang -O3 $(OPT) $(OMP_PORT) $(DEBUG) $(BCPDSRC) -o bcpd -lm -llapack -Wuninitialized
  endif
else
  ifeq ($(OPT),-DNUSE_OPENMP)
	$(CC) -O3 $(DEBUG) $(BCPDSRC) -o bcpd -lm -llapack
  else ## default ##
	clang -O3 -DUSE_OPENMP $(OMP_BREW_ARM) $(DEBUG) $(BCPDSRC) -o bcpd -lm -llapack
  endif
endif

