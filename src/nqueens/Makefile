# nqueens/Makefile
#
# (C) 2018 Carlos M. Fonseca <cmfonsec@dei.uc.pt>
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or (at
# your option) any later version.
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

IDIR = ../include
SDIR = ../solver
CC = gcc
OPT = -O3
CFLAGS = -I$(IDIR)
PDEPS = $(IDIR)/problem.h
SDEPS = $(IDIR)/solver.h
LIBS = -lgsl -lgslcblas -lm

.PHONY: all clean debug

all: nqueens-shd nqueens-sga

clean:
	rm -f *.o *~ core

debug:
	@$(MAKE) OPT=-g

nqueens-shd: nqueens-shd.o nqueens.o
	$(CC) -o $@ $(CFLAGS) $(OPT) $^ $(SDIR)/shd.o $(LIBS)

nqueens-sga: nqueens-sga.o nqueens.o
	$(CC) -o $@ $(CFLAGS) $(OPT) $^ $(SDIR)/sga.o $(LIBS)

%.o: %.c
	$(CC) -c -o $@ $(CFLAGS) $(OPT) $<

nqueens.o: nqueens.c nqueens.h $(PDEPS)
nqueens-shd.o: nqueens-shd.c nqueens.h $(PDEPS) $(IDIR)/shd.h $(SDEPS)
nqueens-sga.o: nqueens-sga.c nqueens.h $(PDEPS) $(IDIR)/sga.h $(SDEPS)

