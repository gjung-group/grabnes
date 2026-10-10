#     MIO Makefile
#
#     Copyright (C) 2012 Rafael Martinez-Gordillo
#
#   This file is distributed under the terms of the
#   GNU General Public License. See the file `LICENSE'
#   in the root directory of the present distribution,
#   or http://www.gnu.org/copyleft/gpl.txt
#

default: libMIO.a

VPATH+=$(SRC_DIR)/MIO

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
MIO_OBJS := io.o sys.o format.o precision.o timer.o debug.o \
            file.o mem.o input.o options.o string.o parser.o \
            omp.o

MPI_COMP := $(findstring MPI,$(findstring -DMPI,$(FCFLAGS)))
ifeq ($(MPI_COMP),MPI)
MPI:= libmpi.a
MPI_LIB:=$(LIB_DIR)/libmpi.a
else
MPI:=
MPI_LIB:=
endif

libMIO.a: printMPI $(MPI) $(MIO_OBJS) mio.o | $(LIB_DIR)
	$(AR) $(ARFLAGS) libMIO.a mio.o $(MIO_OBJS)
	cp *.a $(LIB_DIR)
	cp *.mod $(BUILD_DIR)

$(LIB_DIR):
	@(mkdir $(LIB_DIR))

clean:
	$(RM) *.mod *.o *.a
	$(RM) mem_gen mem.f90
	@if [ -d MPI ] ; then \
	cd MPI ; $(MAKE) clean ; \
	else \
	echo 'Nothing to clean in MPI' ; \
	fi

printMPI:
ifeq ($(MPI_COMP),MPI)
	@echo "MPI compilation of library MIO"
else
	@echo "Serial compilation of library MIO"
endif

MPI_DIR:=$(BUILD_DIR)/MIO/MPI
libmpi.a: | $(MPI_DIR)
	@echo
	@echo "    > MPI library compilation <"
	@echo
	@cp -f $(SRC_DIR)/MIO/MPI/Makefile $(MPI_DIR)
	@(cd $(MPI_DIR) ; $(MAKE) libmpi.a)
	@echo "    ---------------------------"

$(MPI_DIR):
	@(mkdir $(MPI_DIR))

mem.f90: $(SRC_DIR)/MIO/mem_gen.f90 $(SRC_DIR)/MIO/mem_temp.f90
mem.f90: $(SRC_DIR)/MIO/mem_in.f90
	$(F90_SERIAL) -o mem_gen $(SRC_DIR)/MIO/mem_gen.f90
	@./mem_gen $(SRC_DIR)/MIO/mem_in.f90 $(SRC_DIR)/MIO/mem_temp.f90 mem.f90

# PREREQUISITES OF OBJECT FILES
debug.o: file.o io.o
file.o: io.o sys.o format.o
format.o: precision.o io.o
input.o: file.o io.o sys.o options.o mem.o $(SRC_DIR)/MIO/input_aliases.inc
input.o: string.o parser.o
mem.o: format.o sys.o
parser.o: sys.o precision.o format.o io.o string.o
sys.o: io.o format.o ioerror.f90
timer.o: sys.o precision.o mem.o

