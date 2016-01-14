# Makefile for bankers2 utilities. These utilities use counts of all possible
# subsets in banker's sequence order to generate plotting files, e.g. for
# circos plots or the R package VennDiagram plots.
#
#   Author: David Laehnemann <david.laehnemann@hhu.de>
#
#   license: MIT


PROGS    = bankers2circos


all: $(PROGS)

CC       = gcc
CFLAGS   = -g -Wall -O2
LDFLAGS  =
LIBS     = -lm

OBJS     = bankers2circos.o

bankers2circos.o: bankers2circos.c

bankers2circos: $(OBJS)
	$(CC) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)
