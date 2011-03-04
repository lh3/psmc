CC=			gcc
CFLAGS=		-g -Wall -O2
DFLAGS=		
OBJS=		
PROG=		cntcpg fq2psmcfa
LIBS=		-lm -lz

.SUFFIXES:.c .o

.c.o:
		$(CC) -c -I.. $(CFLAGS) $(DFLAGS) $< -o $@

all:$(PROG)

fq2psmcfa:fq2psmcfa.o
		$(CC) fq2psmcfa.o -o $@ $(LIBS)

cntcpg:cntcpg.o
		$(CC) cntcpg.o -o $@ $(LIBS)

clean:
		rm -f *.o a.out* *~ *.a $(PROG)
