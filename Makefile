MODULE_TOPDIR = ../..

PGM = i.hyper.smac

include $(MODULE_TOPDIR)/include/Make/Script.make
include $(MODULE_TOPDIR)/include/Make/Html.make
include $(MODULE_TOPDIR)/include/Make/Other.make
include $(MODULE_TOPDIR)/include/Make/Python.make

MODULES = smac aod wvc radtran o3 utils smac_coef_generator
ETCDIR = $(ETC)/i_hyper_lib
COEFDIR = $(ETCDIR)/COEFS
PYFILES := $(patsubst %,$(ETCDIR)/%.py,$(MODULES))
MODULE_SRCS := $(patsubst %,lib/%.py,$(MODULES))

default: script html $(PYFILES) install_coefs

# Ensure HTML manual is installed
$(HTMLDIR)/$(PGM).html: $(PGM).html
	$(INSTALL_DATA) $(PGM).html $(HTMLDIR)/$(PGM).html

$(ETCDIR):
	$(MKDIR) $@

$(ETCDIR)/%.py: lib/%.py | $(ETCDIR)
	$(INSTALL_DATA) $< $@

# Install SMAC coefficient files
install_coefs: | $(ETCDIR)
	cp -r COEFS $(ETCDIR)/

install:
	cp -r $(ETCDIR) $(INST_DIR)/etc
