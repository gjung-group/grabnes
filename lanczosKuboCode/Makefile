# Locate the solver root from this Makefile rather than from $PWD, so that
# 'make -C <dir>' and 'make -f <dir>/Makefile' work from any directory.
ROOT_DIR:=$(patsubst %/,%,$(dir $(abspath $(lastword $(MAKEFILE_LIST)))))

# Include make.sys (override MAKE_SYS to use a configuration kept elsewhere)
MAKE_SYS?=$(ROOT_DIR)/make.sys
ifeq ($(wildcard $(MAKE_SYS)),)
$(error Build configuration '$(MAKE_SYS)' not found. Copy make.sys.example to make.sys and adjust it)
endif
include $(MAKE_SYS)

# Export LDFLAGS
export LDFLAGS

TIMER:=$(findstring DNOTIMER,,$(FCFLAGS))
ifneq ($(TIMER),DNOTIMER)
TIMER:=$(findstring DTIMER,,$(FCFLAGS))
ifneq ($(TIMER),DTIMER)
FCFLAGS+=-DTIMER
endif
endif

F90_SERIAL?=$(F90)

# BUILD_DIR and BIN_DIR may be overridden for an out-of-tree build
BUILD_DIR?=$(ROOT_DIR)/build
SRC_DIR=$(ROOT_DIR)/Src
LIB_DIR:=$(BUILD_DIR)/lib
VPATH:=$(ROOT_DIR):$(SRC_DIR)
BIN_DIR?=$(ROOT_DIR)/bin

VERSION:=$(word 2,$(shell cat $(ROOT_DIR)/version.info))-rev$(word 4,$(shell cat $(ROOT_DIR)/version.info))
ODATE:=$(subst -, ,$(word 8,$(shell cat $(ROOT_DIR)/version.info)))
YEAR:=$(word 1,$(ODATE))
MONTH:=$(word 2,$(ODATE))
DAY:=$(word 3,$(ODATE))
DATE:=$(MONTH)'-'$(DAY)'-'$(YEAR)

export ROOT_DIR BUILD_DIR SRC_DIR BIN_DIR F90
export FCFLAGS LIBS LIB_DIR VPATH ARFLAGS CPPFLAGS
export VERSION DATE INCFLAGS F90_SERIAL

default: grabnes

all: grabnes

grabnes : | $(BUILD_DIR) $(BIN_DIR)

grabnes: flags
	@echo
	@echo ">> GraBNES compilation <<"
	@echo "Version: "$(VERSION)
	@echo "Updated: "$(DATE)
	@echo
	@cp -f $(SRC_DIR)/Makefile $(BUILD_DIR)
	@(cd $(BUILD_DIR) ; $(MAKE) grabnes)

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)
$(BIN_DIR):
	mkdir -p $(BIN_DIR)

.PHONY: clean flags
clean:
	@if [ -d $(BUILD_DIR) ] ; then \
	cd $(BUILD_DIR) ; $(MAKE) clean ; \
	else \
	echo 'Nothing to clean in BUILD_DIR' ; \
	fi
	@if [ -d $(LIB_DIR) ]; then \
	cd $(LIB_DIR) ; $(RM) *.a ; \
	fi
	@$(RM) $(BIN_DIR)/grabnes

