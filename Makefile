CC=			gcc
CFLAGS=		-g -Wall -O2
DFLAGS=		
OBJS=		khmm.o kmin.o cli.o core.o em.o aux.o
PROG=		psmc
LIBS=		-lm -lz

.SUFFIXES:.c .o

.c.o:
		$(CC) -c $(CFLAGS) $(DFLAGS) $< -o $@

all:$(PROG)

psmc:$(OBJS) main.o
		$(CC) $(CCFLAGS) $(DFLAGS) $(OBJS) main.o -o $@ $(LIBS)

khmm.o:khmm.h
kmin.o:kmin.h
cli.o core.o aux.o em.o:psmc.h khmm.h

clean:
		rm -f *.o a.out *~ *.a $(PROG)
