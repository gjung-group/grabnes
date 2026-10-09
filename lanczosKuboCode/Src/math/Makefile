#     Math library Makefile
#
#     Copyright (C) 2012 Rafael Martinez-Gordillo
#
#   This file is distributed under the terms of the
#   GNU General Public License. See the file `LICENSE'
#   in the root directory of the present distribution,
#   or http://www.gnu.org/copyleft/gpl.txt
#

default: libmath.a

VPATH+=$(SRC_DIR)/math

%.o:%.f90
	$(F90) -c $(FCFLAGS) $(INCFLAGS) $< -o $@
%.o:%.F90
	$(F90) -c $(FCFLAGS) $(INCFLAGS) $(CPPFLAGS) $< -o $@
%.o:%.mod

.PHONY: clean

# LIST OF OBJECT FILES
MATH_OBJS := math.o adams.o lagrangediff.o trapezoidal.o numerov.o vectors.o \
             qsort.o

libmath.a: $(MATH_OBJS) | $(LIB_DIR)
	@$(AR) $(ARFLAGS) libmath.a $(MATH_OBJS)
	@echo $(LIB_DIR)
	@cp *.a $(LIB_DIR)
	@cp *.mod $(BUILD_DIR)

$(LIB_DIR):
	@(mkdir $(LIB_DIR))

clean:
	$(RM) *.mod *.o *.a

