# Minimal Makefile
#CFLAGS=-O
CC=gcc

# Adapt to location of SAC distribution 
SACDIR=/usr/local/sac
SACLIB = -lsacio
CFLAGS = -g -I$(SACDIR)/include  -L$(SACDIR)/lib  

LIBS = 
GSLLIBS = -lgsl -lgslcblas -lm

# Executables will be placed here (NB: no separate installation, final link directly generates executable here)
BIN=$(HOME)/bin/$(ARCH)
#BIN = .

EXEC = multisplit split_cor error_stack

.PHONY:  install clean
default: $(EXEC)

multisplit: multisplit.o  rmeantaper.o gsl_seis.o invfisher.o betai.o sac_help.o
	$(CC) $(CFLAGS) -o $@ $^ $(SACLIB) $(GSLLIBS) $(LIBS)

split_cor: split_cor.o gsl_seis.o sac_help.o
	$(CC) $(CFLAGS) -o $@ $^ $(SACLIB) $(GSLLIBS) $(LIBS)

error_stack: error_stack.o invfisher.o betai.o 
	$(CC) $(CFLAGS) -o $@ $^ $(SACLIB) $(GSLLIBS) $(LIBS)
# betacf.o gammln.o nrutil.o 

gsl_seis.o: gsl_seis.h

multisplit.o: multisplit.h

clean:
	rm *.o $(EXEC)
	
install: $(EXEC)
	install $(EXEC) $(BIN)
