MODULE_TOPDIR = ../..

PGM = i.hyper.smac

include $(MODULE_TOPDIR)/include/Make/Script.make
include $(MODULE_TOPDIR)/include/Make/Html.make
include $(MODULE_TOPDIR)/include/Make/Other.make
include $(MODULE_TOPDIR)/include/Make/Python.make

MODULES = smac aod wvc radtran o3
ETCDIR = $(ETC)/i_hyper_lib
PYFILES := $(patsubst %,$(ETCDIR)/%.py,$(MODULES))
MODULE_SRCS := $(patsubst %,lib/%.py,$(MODULES))

default: script html $(PYFILES)

# Ensure HTML manual is installed
$(HTMLDIR)/$(PGM).html: $(PGM).html
	$(INSTALL_DATA) $(PGM).html $(HTMLDIR)/$(PGM).html

$(ETCDIR):
	$(MKDIR) $@

$(ETCDIR)/%.py: lib/%.py | $(ETCDIR)
	$(INSTALL_DATA) $< $@

install:
	cp -r $(ETCDIR) $(INST_DIR)/etc
