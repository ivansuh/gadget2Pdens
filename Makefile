# Makefile

CC=gcc
DFLAGS=-std=c99 -O2  -DHAVE_STDINT_H -DHAVE_STDBOOL_H -DNONC99=0 -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE #-DVELOCITY -DHEADER 

OBJ = utility.o gadget2Pdens.o 
PRG = gadget2Pdens

all: clean $(OBJ) $(PRG)

$(OBJ): %.o: %.c
	$(CC) -c $(DFLAGS) $<

$(PRG):
	$(CC) $(OBJ) -o $@ $(OPTFLAGS) $(DFLAGS) -lm

clean:
	rm -f *.o $(PRG)
