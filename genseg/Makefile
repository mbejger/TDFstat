CC = gcc
#LAL_FLAGS = -DUSE_LAL -I/usr/share/lalpulsar/include -I/usr/include/lal/include
LAL_FLAGS = 
#LAL_LIBS = -llalpulsar -L/usr/share/lalpulsar/lib
LAL_LIBS =
INIPARSER_DIR = ../utils/iniparser
INIPARSER_FLAGS = -I$(INIPARSER_DIR)/src
INIPARSER_LIBS = -L$(INIPARSER_DIR) -liniparser
FFTW_DIR = /opt/fftw/3.3.10-gcc9-mvapich2.2
CFLAGS = -Wall $(LAL_FLAGS) $(INIPARSER_FLAGS) -O2 -DDEBUG
# debug flags below
#CFLAGS = -Wall $(LAL_FLAGS) $(INIPARSER_FLAGS) -DDEBUG -g
LOADLIBES = -static -lfftw3 -lgsl -lgslcblas $(INIPARSER_LIBS) \
	    -Wl,-Bdynamic -lc -lm  $(LAL_LIBS) 

.PHONY:	iniparser

all:	iniparser genseg

iniparser:	
	$(MAKE) -C $(INIPARSER_DIR)

genseg:	genseg.o GrubbsOutliers.o get_barycenter.o sid.o
	$(CC) -o genseg  $^ $(LOADLIBES)

laltest: laltest.o get_barycenter.o sid.o
	$(CC) -o laltest $^ $(LAL_LIBS)

clean:
	rm -f *.o genseg laltest
	$(MAKE) -C $(INIPARSER_DIR) veryclean
