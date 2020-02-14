# solver/Makefile
#
# (C) 2018 Carlos M. Fonseca <cmfonsec@dei.uc.pt>
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License, version 3, as
# published by the Free Software Foundation.
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
CC = gcc
OPT = -O3 -g
CFLAGS = -I$(IDIR)
PDEPS = $(IDIR)/problem.h
SDEPS = $(IDIR)/solver.h
LIBS = -lgsl -lgslcblas -lm

.PHONY: all clean

all: sga.o shd.o ils.o

clean:
	rm -f *.o *~ core

%.o: %.c
	$(CC) -c -o $@ $(CFLAGS) $(OPT) $<

shd.o : shd.c  $(PDEPS) $(IDIR)/shd.h $(SDEPS)
sga.o : sga.c  $(PDEPS) $(IDIR)/sga.h $(SDEPS)
ils.o: ils.c $(PDEPS) $(IDIR)/ils.h $(SDEPS)

