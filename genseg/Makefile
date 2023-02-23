CC = gcc
#LAL_FLAGS = -DUSE_LAL -I$(CONDA_PREFIX)/include
LAL_FLAGS = 
#LAL_LIBS = -L$(CONDA_PREFIX)/lib   -llalpulsar -llalsupport -llalmetaio -llal
LAL_LIBS =
INIPARSER_DIR = ../utils/iniparser
INIPARSER_FLAGS = -I$(INIPARSER_DIR)/src
INIPARSER_LIBS = -L$(INIPARSER_DIR) -liniparser
FFTW_DIR = /opt/fftw/3.3.10-gcc9-mvapich2.2
CFLAGS = -Wall $(LAL_FLAGS) $(INIPARSER_FLAGS) -O2 -DDEBUG
# debug flags below
#CFLAGS = -Wall $(LAL_FLAGS) $(INIPARSER_FLAGS) -DDEBUG -g
LOADLIBES = -Wl,--dynamic-linker=/lib64/ld-linux-x86-64.so.2 \
	    -static -lfftw3 -lgsl -lgslcblas  $(INIPARSER_LIBS) $(LAL_LIBS) \
	    -Wl,-Bdynamic -lc -lm -lz -lpthread

.PHONY:	iniparser

all:	iniparser genseg

iniparser:	
	$(MAKE) -C $(INIPARSER_DIR)

genseg:	genseg.o GrubbsOutliers.o get_barycenter.o sid.o
	$(CC) -o genseg  $^ $(LOADLIBES)

laltest: laltest.o get_barycenter.o sid.o
	$(CC) -o laltest $^ $(LOADLIBES)

clean:
	rm -f *.o genseg laltest
	$(MAKE) -C $(INIPARSER_DIR) veryclean
