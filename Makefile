CC= bspcc
CFLAGS= -std=c99 -Wall -O3
#Replace -O0 -g by -O3 for the final version.code
LFLAGS= -lm

OBJMULT= mult.o
OBJPARMULT= parmult.o bspfft.o bspedupack.o

all: mult parmult

mult: $(OBJMULT)
	$(CC) $(CFLAGS) -o mult $(OBJMULT) $(LFLAGS)
parmult: $(OBJPARMULT)
	$(CC) $(CFLAGS) -o parmult $(OBJPARMULT) $(LFLAGS)

parmult.o: bspedupack.h
mult.o:
bspfft.o:    bspedupack.h
bspedupack.o: bspedupack.h

.PHONY: clean
clean:
	rm -f $(OBJMULT) $(OBJPARMULT) mult parmult