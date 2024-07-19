#
INST_DIR = $(HOME)/bin
BIN  = aftani_c_pgl_test

fflags = -O2 -Wall -ffixed-line-length-none 
cflags = -O2 -Wall 

LDLIBS =   -L/uufs/chpc.utah.edu/sys/installdir/fftw/3.3.6/lib -lfftw3
FFLAGS = $(DBG) $(fflags) -I/uufs/chpc.utah.edu/sys/installdir/fftw/3.3.6/include

#LDLIBS =    -L/sw/lib -lfftw3
#LDLIBS =   -L/uufs/chpc.utah.edu/common/home/u0919412/Research/lib -lfftw3

FFLAGS = $(DBG) $(fflags)

CFLAGS = $(DBG) $(cflags)

FC = gfortran

CC = gcc

DBG =


FOBJS = driver_c.o swapn.o aftanpg.o aftanipg.o ftfilt.o fmax.o taper.o \
	trigger.o tapers.o tgauss.o mspline.o pred_cur.o misc.o phtovel.o

$(BIN) : $(FOBJS)
	$(FC) $(FFLAGS) $(CFLAGS) $(FOBJS)  -o $(BIN) $(LDLIBS)

install : $(BIN)
	install -s $(BIN) $(INST_DIR)

clean :
	rm -f $(BIN) core $(FOBJS)
