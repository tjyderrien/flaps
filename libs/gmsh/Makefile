#----------------------------------------------------------------------------#
# Copyright (c) 2014 Alexandre MOUTON alexandre.mouton@math.univ-lille1.fr   #
#                                                                            #
# This file is part of MESH2VF.                                              #
#                                                                            #
# MESH2VF is free software: you can redistribute it and/or modify it under   #
# the terms of the GNU General Public License as published by the Free       #
# Software Foundation, either version 3 of the License, or at your option)   #
# any later version.                                                         #
#                                                                            #
# MESH2VF is distributed in the hope that it will be useful, but WITHOUT ANY #
# WARRANTY; without even the implied warranty of MERCHANTABILITY or          #
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for   #
# more details.                                                              #
#                                                                            #
# You should have received a copy of the GNU General Public License along    #
# with MESH2VF. If not, see <http://www.gnu.org/licenses/>.                  #
#----------------------------------------------------------------------------#
LIB		= 
INC		=

F90FLAGS	= -O2 -C -g -pg #-w -backtrace
F90			= gfortran #ifort 

LD		= $(F90)


# Declaration of the objects
#---------------------------
EXEC     = msh2vf
OBJ_EXEC = $(addsuffix .o,$(EXEC))
SRC_EXEC = $(addsuffix .f90,$(EXEC))
MODULE   = libmsh2vf msh2vf
OBJ_MOD  = $(addsuffix .o,$(MODULE))
SRC_MOD  = $(addsuffix .f90,$(MODULE))
OBJ = $(OBJ_MOD)
SRC = $(SRC_MOD)


# Compilation
#------------
all: clean $(OBJ) $(EXEC)

$(EXEC) : $(OBJS)
	@echo "Loading $@ ..."
	@$(F90) $(INC) $(F90FLAGS) $(OBJ) $(LIB)  -o $(EXEC) 
	@echo "done $@"

$(OBJ): %.o:	%.f90
	$(F90) $(F90FLAGS) -c $< -o $@ $(INC)

clean:
	@rm -f *.o
	@rm -f $(EXEC)
	@rm -f *.mod

