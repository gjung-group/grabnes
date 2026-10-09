# Kubo3 Makefile
#
# Copyright (C) 2012 Rafael Martinez-Gordillo
#

# The build configuration (make.sys) is read once by the root Makefile, which
# exports the compiler, flag, and directory variables used here.

# Optional module compiled only when FCFLAGS contains -DSEMICL and its sources are present
SEMICL_OBJ := $(if $(findstring -DSEMICL,$(FCFLAGS)),semicl.o,)

default: grabnes

# Compilation rules
%.o: %.f90
	$(F90) -c $(FCFLAGS) $(INCFLAGS) $< -o $@
%.o: %.F90
	$(F90) -c $(FCFLAGS) $(INCFLAGS) $(CPPFLAGS) $< -o $@
%.o: %.mod

ifndef ROOT_DIR
$(error Run make from root directory)
endif

.PHONY: clean libMIO.a

MIO_DIR := $(BUILD_DIR)/MIO
MATH_DIR := $(BUILD_DIR)/math

libMIO.a: | $(MIO_DIR)
	@echo
	@echo "  >>> MIO library compilation <<<"
	@echo
	@cp -f $(SRC_DIR)/MIO/Makefile $(MIO_DIR)
	@(cd $(MIO_DIR); $(MAKE) libMIO.a)
	@echo " ==============================="
	@echo

libmath.a: | $(MATH_DIR)
	@echo
	@echo "  >>> Math library compilation <<<"
	@echo
	@cp -f $(SRC_DIR)/math/Makefile $(MATH_DIR)
	@(cd $(MATH_DIR); $(MAKE) libmath.a)
	@echo " ==============================="
	@echo

$(MIO_DIR):
	@(mkdir $(MIO_DIR))

$(MATH_DIR):
	@(mkdir $(MATH_DIR))

LIBS += $(LIB_DIR)/libMIO.a $(LIB_DIR)/libmath.a

objprint:
	@echo
	@echo " >>> Object files compilation <<<"
	@echo

OBJS := grabnes_start.o version.o grabnes_end.o name.o constants.o \
        atoms.o cell.o calc.o ham.o parallel.o tbpar.o neigh.o random.o \
        kuboarrays.o kubo.o kubosubs.o gauss.o magf.o diag.o interface.o hybrid.o \
        scf.o moireBLShift.o $(SEMICL_OBJ)

# Linking step
grabnes: $(SRC_DIR)/grabnes.f90 libMIO.a libmath.a objprint $(OBJS)
	$(F90) -o $(BIN_DIR)/grabnes $(SRC_DIR)/grabnes.f90 $(OBJS) $(LDFLAGS) $(LIBS)

version.o: version.info $(SRC_DIR)/version.f90
	@echo ">> Updating version in file version.f90 <<"; \
	exec 3<$(SRC_DIR)/version.f90; \
	exec 4>version.f90; \
	while read line <&3; do \
	echo "$$line" | sed 's/KUBO_VERSION/$(VERSION)/g;s/KUBO_UPDATE/$(DATE)/g' >&4; \
	done;
	$(F90) -c $(FCFLAGS) $(INCFLAGS) version.f90

clean:
	$(RM) *.mod *.o *.a
	$(RM) $(BUILD_DIR)/version.f90
	@if [ -f MIO/Makefile ]; then \
	cd MIO; $(MAKE) clean; \
	else \
	echo 'Nothing to clean in MIO'; \
	fi
	@if [ -f math/Makefile ]; then \
	cd math; $(MAKE) clean; \
	else \
	echo 'Nothing to clean in math'; \
	fi

# LIST OF OBJECT FILES
grabnes_start.o: version.o name.o random.o
atoms.o: constants.o cell.o parallel.o hybrid.o
calc.o: ham.o tbpar.o kuboarrays.o kubo.o kubosubs.o magf.o moireBLShift.o diag.o $(SEMICL_OBJ)
semicl.o: constants.o cell.o name.o
cell.o: constants.o
diag.o: cell.o atoms.o ham.o neigh.o constants.o name.o tbpar.o scf.o
gauss.o: atoms.o cell.o random.o
ham.o: atoms.o tbpar.o cell.o constants.o neigh.o gauss.o magf.o moireBLShift.o
ham.o: interface.o random.o
hybrid.o: random.o cell.o constants.o
interface.o: atoms.o cell.o
kubo.o: constants.o atoms.o kuboarrays.o
kubo.o: kubosubs.o random.o
kubosubs.o: kuboarrays.o constants.o atoms.o neigh.o
magf.o: cell.o constants.o
moireBLShift.o: cell.o constants.o
random.o: parallel.o name.o
scf.o: cell.o atoms.o neigh.o ham.o
