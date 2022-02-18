TARGET = gwsearch-omp-cascadelake-dev9-newfft-sleef

VER := $(shell git rev-parse HEAD)

# FFTW_DIR - fftw location, remember to compile it with optimizations (e.g. avx)
# SINCOS = {YEPPP or SLEEF or GNUSINCOS or NOSINCOS}
#FFTW_DIR = /opt/fftw/3.3.7-gcc7.2-mvapich2.2
FFTW_DIR = /opt/fftw/3.3.10-gcc9-mvapich2.2
SINCOS =  SLEEF #YEPPP or SLEEF or GNUSINCOS or NOSINCOS

CC = gcc -std=gnu11 -static-libgcc -fopenmp

# pre-processor flags: TIMERS, VERBOSE (rich stdout output)
CFLAGS = -DPREFIX="./candidates" -DCODEVER=$(VER) -DTIMERS=1 -D$(SINCOS) -UUSE_LOCKING #-DVERBOSE 
CFLAGS += -I$(FFTW_DIR)/include 
CFLAGS += -g -Wall -Wno-unused -fno-omit-frame-pointer \
	  -O3 -funsafe-loop-optimizations -funroll-loops \
	  -march=broadwell -mtune=broadwell  \
	  -ffast-math
#	  -march=cascadelake -mtune=cascadelake  \
#	  -march=broadwell -mtune=broadwell -mavx2 -msse4.2 -mfma \
# 	  -march=sandybridge -mtune=sandybridge -mavx
#	  -march=haswell -mtune=haswell -mavx2 -mfma

LDFLAGS = -L$(FFTW_DIR)/lib

LDLIBS = -Wl,--dynamic-linker=/lib64/ld-linux-x86-64.so.2 \
	 -static -lfftw3 -lfftw3_omp -lgsl -lgslcblas \
	 -Wl,-Bdynamic -lc -lm -lrt

ifeq ($(strip $(SINCOS)),SLEEF)
SPATH = lib/sleef-3.5.1
CFLAGS  += -DENABLE_AVX -I$(SPATH)/include
LDFLAGS += -L$(SPATH)/lib64
LDLIBS += -lsleef
VECLIB = $(SPATH)/lib/libsleef.so
endif
ifeq ($(strip $(SINCOS)),YEPPP)
CFLAGS  += -Ilib/yeppp-1.0.0/library/headers
LDFLAGS += -Llib/yeppp-1.0.0/binaries/linux/x86_64
LDLIBS += -lyeppp
VECLIB = lib/yeppp-1.0.0/binaries/linux/x86_64/libyeppp.so
endif


all: veclib $(TARGET)

ifeq ($(strip $(SINCOS)),SLEEF)
veclib:
	@echo -e "\nUsing prebuild SLEEF library from: $(VECLIB)\n"
endif
ifeq ($(strip $(SINCOS)),YEPPP)
veclib:
	@echo -e "\nUsing prebuild YEPPP library from: $(VECLIB)\n"
endif


$(TARGET): main.o jobcore.o timer.o settings.o auxi.o init.o
	$(CC) $(LDFLAGS) -o $@ $^ $(LDLIBS)

jobcore.o: settings.h timer.h auxi.h struct.h

main.o: settings.h struct.h init.h auxi.h

init.o: settings.h struct.h init.h auxi.h

settings.o: settings.h struct.h

auxi.o: auxi.c auxi.h

clean:
	rm -f *.o

uninstall: clean
	rm -f $(TARGET)

.PHONY: all veclib clean uninstall

#----------------------------------------------

sigen: sigen.o settings.o auxi.o
	$(CC) -o $@ $^ -Wl,-Bdynamic -lc -lrt -lm
