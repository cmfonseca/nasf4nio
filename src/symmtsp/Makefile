# symmtsp/Makefile
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
OPT = -O0 -g
CFLAGS = -I$(IDIR) -Wall
PDEPS = $(IDIR)/problem.h
SDEPS = $(IDIR)/solver.h
LIBS = -lgsl -lgslcblas -lm

.PHONY: all clean debug

all: symmtsp-shd symmtsp-ils symmtsp-sga  symmtsp-test
clean:
	rm -f *.o *~ core

debug:
	@$(MAKE) OPT=-g


symmtsp-shd: symmtsp-shd.o symmtsp.o
	$(CC) -o $@ $(CFLAGS) $(OPT) $^ $(SDIR)/shd.o $(LIBS)

symmtsp-ils: symmtsp-ils.o symmtsp.o
	$(CC) -o $@ $(CFLAGS) $(OPT) $^ $(SDIR)/ils.o $(LIBS)
	
symmtsp-sga: symmtsp-sga.o symmtsp.o
	$(CC) -o $@ $(CFLAGS) $(OPT) $^ $(SDIR)/sga.o $(LIBS)

symmtsp-test: symmtsp-test.o symmtsp.o
	$(CC) -o $@ $(CFLAGS) $(OPT) $^ $(SDIR)/sga.o $(LIBS)

%.o: %.c
	$(CC) -c -o $@ $(CFLAGS) $(OPT) $<

symmtsp.o: symmtsp.c symmtsp.h $(PDEPS)
symmtsp-shd.o: symmtsp-shd.c symmtsp.h $(PDEPS) $(IDIR)/shd.h $(SDEPS)
symmtsp-ils.o: symmtsp-ils.c symmtsp.h $(PDEPS) $(IDIR)/ils.h $(SDEPS)
symmtsp-sga.o: symmtsp-sga.c symmtsp.h $(PDEPS) $(IDIR)/sga.h $(SDEPS)
symmtsp-test.o: symmtsp-test.c symmtsp.h $(PDEPS) $(IDIR)/sga.h $(SDEPS)
