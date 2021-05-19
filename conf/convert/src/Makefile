CC = gcc
#CC = /opt/ibmcmp/vac/7.0/bin/xlc
LD = gcc
#LD = /opt/ibmcmp/xlf/9.1/bin/xlf
default: conv_lat
clean: 
	rm *.o
distclean: 
	rm *.o *~
SRC=comp_su3.c io.c conv_lat.c
OBJ=$(SRC:.c=.o)

CBLASDIR=usr/local/sw/CBLAS
CFLAGS=-I./

LDFLAGS=-L/$(CBLASDIR)/lib -lcblas
#LDFLAGS=-L./

LOADLIBES= -llapack -lblas -lm
conv_lat: $(OBJ)
	$(LD) $(LDFLAGS) $(OBJ) $(LOADLIBES) -o conv_lat