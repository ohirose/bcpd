CC=gcc
BCPDSRC=  register/*.c base/*.c
OMP_PORT= -Xpreprocessor -fopenmp -I/opt/local/include/libomp /opt/local/lib/libomp/libomp.dylib
OMP_BREW= -Xpreprocessor -fopenmp -I/usr/local/include/ /usr/local/lib/libomp.dylib
DEBUG= -g -Wall

all:
ifeq ($(OPT),-DUSE_OPENMP)
  ifeq ($(ENV),LINUX)
	$(CC) -O3 -fopenmp $(OPT) $(DEBUG) $(BCPDSRC) -o bcpd -lm -llapack
  else ifeq ($(ENV),MINGW32)
	$(CC) -O3 -fopenmp $(OPT) $(DEBUG) $(BCPDSRC) win/*.dll -o win/bcpd -lm -DMINGW32
  else ifeq ($(ENV),HOMEBREW)
	clang -O3 $(OPT) $(OMP_BREW) $(DEBUG) $(BCPDSRC) -o bcpd -lm -llapack -Wuninitialized
  else ifeq ($(ENV),MACPORTS)
	clang -O3 $(OPT) $(OMP_PORT) $(DEBUG) $(BCPDSRC) -o bcpd -lm -llapack
  endif
else
  ifeq ($(OPT),-DNUSE_OPENMP)
	$(CC) -O3 $(DEBUG) $(BCPDSRC) -o bcpd -lm -llapack
  else ## default ##
	clang -O3 -DUSE_OPENMP $(OMP_PORT) $(DEBUG) $(BCPDSRC) -o bcpd -lm -llapack
  endif
endif

