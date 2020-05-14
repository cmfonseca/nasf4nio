# Makefile
#
# (C) 2020 Carlos M. Fonseca <cmfonsec@dei.uc.pt>
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

.PHONY: all clean debug

all:
	@$(MAKE) -C solver
	@$(MAKE) -C nqueens
	@$(MAKE) -C qap
#	@$(MAKE) -C flowshop
#	@$(MAKE) -C knapsack
	@$(MAKE) -C symmtsp
	
clean:
	@$(MAKE) -C solver clean
	@$(MAKE) -C nqueens clean
	@$(MAKE) -C qap clean
	@$(MAKE) -C flowshop clean
	@$(MAKE) -C knapsack clean
	@$(MAKE) -C symmtsp clean
	rm -f *~ core include/*~

debug:
	@$(MAKE) -C solver debug
	@$(MAKE) -C nqueens debug
	@$(MAKE) -C qap debug
#	@$(MAKE) -C flowshop debug
#	@$(MAKE) -C knapsack debug
	@$(MAKE) -C symmtsp


