# qap/Makefile
#
# (C) 2018 Carlos M. Fonseca <cmfonsec@dei.uc.pt>
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License, version 3, as
# published bythe Free Software Foundation.
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
OPT = -g -O3
CFLAGS = -I$(IDIR)
PDEPS = $(IDIR)/problem.h
SDEPS = $(IDIR)/solver.h
LIBS = -lgsl -lgslcblas -lm

.PHONY: all clean debug

all: qap-shd qap-ils qap-sga

clean:
	rm -f *.o *~ core

debug:
	@$(MAKE) OPT=-g

qap-shd: qap-shd.o qap.o $(SDIR)/shd.o
	$(CC) -o $@ $(CFLAGS) $(OPT) $^ $(LIBS)

qap-ils: qap-ils.o qap.o $(SDIR)/ils.o
	$(CC) -o $@ $(CFLAGS) $(OPT) $^ $(LIBS)
	
qap-sga: qap-sga.o qap.o $(SDIR)/sga.o
	$(CC) -o $@ $(CFLAGS) $(OPT) $^ $(LIBS)


%.o: %.c
	$(CC) -c -o $@ $(CFLAGS) $(OPT) $<

qap.o: qap.c qap.h $(PDEPS)
qap-shd.o: qap-shd.c qap.h $(PDEPS) $(IDIR)/shd.h $(SDEPS)
qap-ils.o: qap-ils.c qap.h $(PDEPS) $(IDIR)/ils.h $(SDEPS)
qap-sga.o: qap-sga.c qap.h $(PDEPS) $(IDIR)/sga.h $(SDEPS)

