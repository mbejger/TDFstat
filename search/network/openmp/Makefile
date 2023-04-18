# min. required for SLEEF is haswell (AVX2+FMA)
ARCH = cascadelake
TARGET = gwsearch-omp-$(ARCH)

VER := $(shell git rev-parse HEAD)

# FFTW_DIR - fftw location
FFTW_DIR = /opt/fftw/3.3.10-gcc9-mvapich2.2
#FFTW_DIR = $(CONDA_PREFIX)
SINCOS = SLEEF #SLEEF or GNUSINCOS

CC = gcc -std=gnu11 -static-libgcc -fopenmp #-flto
#-flto ?? speeds up a little for AVX512 but slows down for AVX2 

# pre-processor flags: TIMERS, VERBOSE (rich stdout output)
CFLAGS = -DPREFIX="./candidates" -DCODEVER=$(VER) -DTIMERS=1 -D$(SINCOS) -UUSE_LOCKING #-DVERBOSE 
CFLAGS += -I$(FFTW_DIR)/include 
CFLAGS += -g -Wall -Wno-unused \
	  -O3 -funsafe-loop-optimizations -funroll-loops -ffast-math \
	  -march=$(ARCH) -mtune=$(ARCH)

LDFLAGS = -L$(FFTW_DIR)/lib

LDLIBS = -Wl,--dynamic-linker=/lib64/ld-linux-x86-64.so.2 \
	 -static -lfftw3 -lfftw3_omp -lgsl -lgslcblas \
	 -Wl,-Bdynamic -lc -lm -lrt

ifeq ($(strip $(SINCOS)),SLEEF)
SPATH = lib/sleef-3.5.1
# supported options are: USE_AVX2 , USE_AVX512 , USE_AVX2_SINCOS
CFLAGS  += -DUSE_AVX2 -I$(SPATH)/include
LDFLAGS += -L$(SPATH)/lib64
LDLIBS += -lsleef
endif


all: veclib $(TARGET)

ifeq ($(strip $(SINCOS)),SLEEF)
veclib:
	@echo -e "\nUsing prebuild SLEEF library from: $(SPATH)\n"
endif

SRC = main.c jobcore.c auxi.c settings.c init.c  timer.c spinmod.c
INC = jobcore.h settings.h auxi.h timer.h init.h struct.h
DEPS = $(SRC) $(INC) Makefile

$(TARGET): $(DEPS)
	$(CC) $(LDFLAGS) -o $@ $(SRC) $(CFLAGS) $(LDLIBS)

clean:
	rm -f *.o 

uninstall: clean
	rm -f $(TARGET)

.PHONY: all $(TARGET) veclib clean uninstall

#----------------------------------------------

sigen: sigen.o settings.o auxi.o
	$(CC) -o $@ $^ -Wl,-Bdynamic -lc -lrt -lm
