CC= bspcc
CFLAGS= -std=c99 -Wall -O3
LFLAGS= -lm

OBJMULT= mult.o bspfft.o bspedupack.o

all: mult

mult: $(OBJMULT)
	$(CC) $(CFLAGS) -o mult $(OBJMULT) $(LFLAGS)

mult.o: bspedupack.h
bspfft.o:    bspedupack.h
bspedupack.o: bspedupack.h

.PHONY: clean
clean:
	rm -f $(OBJMULT) mult