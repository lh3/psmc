CC=			gcc
CFLAGS=		-g -Wall -O2 -Wno-unused-function
DFLAGS=		
OBJS=		
PROG=		cntcpg fq2psmcfa splitfa calD mutDiff
LIBS=		-lm -lz

.SUFFIXES:.c .o

.c.o:
		$(CC) -c -I.. $(CFLAGS) $(DFLAGS) $< -o $@

all:$(PROG)

fq2psmcfa:fq2psmcfa.o
		$(CC) $< -o $@ $(LIBS)

splitfa:splitfa.o
		$(CC) $< -o $@ $(LIBS)

cntcpg:cntcpg.o
		$(CC) $< -o $@ $(LIBS)

calD:calD.o
		$(CC) $< -o $@ $(LIBS)

mutDiff:mutDiff.o
		$(CC) $< -o $@ $(LIBS)

clean:
		rm -f *.o a.out* *~ *.a $(PROG)
