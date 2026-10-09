#     MPI Makefile
#
#     Copyright (C) 2012 Rafael Martinez-Gordillo
#
#   This file is distributed under the terms of the
#   GNU General Public License. See the file `LICENSE'
#   in the root directory of the present distribution,
#   or http://www.gnu.org/copyleft/gpl.txt
#

default: libmpi.a

VPATH+=$(SRC_DIR)/MIO/MPI

%.o:%.f90
	$(F90) -c $(FCFLAGS) $(INCFLAGS) $< -o $@
%.o:%.F90
	$(F90) -c $(FCFLAGS) $(INCFLAGS) $(CPPFLAGS) $< -o $@
%.o:%.mod

.PHONY: clean

ifndef ROOT_DIR
$(error Run make from root directory)
endif

# LIST OF OBJECT FILES
MPI_OBJS := mpi.o mpitime.o mpi_inc.o mpikinds.o mpi_types.o

libmpi.a: mpi_type.f90 $(MPI_OBJS)
	@$(AR) $(ARFLAGS) libmpi.a $(MPI_OBJS)
	@cp *.a ../libMIO.a
	@cp *.mod $(BUILD_DIR)
	@cp *.mod $(BUILD_DIR)/MIO

mpi_types.f90: $(SRC_DIR)/MIO/MPI/type_gen.f90 $(SRC_DIR)/MIO/MPI/mpi_type.f90
	$(F90_SERIAL) -o type_gen $(SRC_DIR)/MIO/MPI/type_gen.f90
	@./type_gen $(SRC_DIR)/MIO/MPI/mpi_type.f90

clean:
	$(RM) *.mod *.o *.a
	$(RM) type_gen mpi_types.f90

# PREREQUISITES OF OBJECT FILES
mpi.o: mpitime.o mpi_inc.o mpi_types.o
mpi_types.o: mpitime.o mpi_inc.o mpikinds.o
mpitime.o: mpikinds.o