# Makefile for i.hyper.smac GRASS GIS module
# Supports multiple GRASS versions (8.5, 8.6, etc.)

# Auto-detect GRASS installation directory
ifeq ($(wildcard $(HOME)/dev/grass),)
    ifeq ($(wildcard /usr/local/grass86),)
        ifeq ($(wildcard /usr/local/grass85),)
            $(error "GRASS installation not found. Please install GRASS or set MODULE_TOPDIR")
        else
            MODULE_TOPDIR = /usr/local/grass85
        endif
    else
        MODULE_TOPDIR = $(HOME)/dev/grass
    endif
else
    MODULE_TOPDIR = $(HOME)/dev/grass
endif

# Module information
PGM = i.hyper.smac

# Include GRASS makefiles
include $(MODULE_TOPDIR)/include/Make/Script.make
include $(MODULE_TOPDIR)/include/Make/Html.make
include $(MODULE_TOPDIR)/include/Make/Other.make
include $(MODULE_TOPDIR)/include/Make/Python.make

# Module components
MODULES = smac aod wvc radtran o3 utils smac_coef_generator lut spectral_polish adjacency gas_absorption parallel_lut opencl_accelerator
ETCDIR = $(ETC)/i_hyper_lib
COEFDIR = $(ETCDIR)/COEFS
PYFILES := $(patsubst %,$(ETCDIR)/%.py,$(MODULES))
MODULE_SRCS := $(patsubst %,lib/%.py,$(MODULES))

# Default target - install modules first, then generate docs
default: $(PYFILES) install_coefs script html

# Ensure HTML manual is installed
$(HTMLDIR)/$(PGM).html: $(PGM).html
	$(INSTALL_DATA) $(PGM).html $(HTMLDIR)/$(PGM).html

# Ensure MD manual is installed
$(HTMLDIR)/$(PGM).md: $(PGM).md
	$(INSTALL_DATA) $(PGM).md $(HTMLDIR)/$(PGM).md

# Create library directory
$(ETCDIR):
	$(MKDIR) $@

# Install Python modules
$(ETCDIR)/%.py: lib/%.py | $(ETCDIR)
	$(INSTALL_DATA) $< $@

# Install SMAC coefficient files
install_coefs: | $(ETCDIR)
	cp -r COEFS $(ETCDIR)/

# Install all components
install: script html $(PYFILES) install_coefs
	# Copy library to installation directory
	cp -r $(ETCDIR) $(INST_DIR)/etc/
	# Install documentation
	$(INSTALL_DATA) $(PGM).html $(HTMLDIR)/$(PGM).html
	$(INSTALL_DATA) $(PGM).md $(HTMLDIR)/$(PGM).md
	# Install Python script without .py extension (GRASS policy)
	$(INSTALL) $(PGM).py $(SCRIPTDIR)/$(PGM)

# Clean build artifacts
clean:
	rm -rf $(ETCDIR)
	rm -f $(HTMLDIR)/$(PGM).html
	rm -f $(HTMLDIR)/$(PGM).md

# Uninstall from system
uninstall:
	rm -rf $(INST_DIR)/etc/i_hyper_lib
	rm -f $(HTMLDIR)/$(PGM).html
	rm -f $(HTMLDIR)/$(PGM).md
	rm -f $(BIN)/$(PGM)

# Install to multiple GRASS versions
install-all: install
	@if [ -d "/usr/local/grass85" ] && [ "$(MODULE_TOPDIR)" != "/usr/local/grass85" ]; then \
		echo "Installing to GRASS 8.5..."; \
		cp -r $(ETCDIR) /usr/local/grass85/etc/; \
		$(INSTALL_DATA) $(PGM).html /usr/local/grass85/docs/html/$(PGM).html; \
		$(INSTALL_DATA) $(PGM).md /usr/local/grass85/docs/html/$(PGM).md; \
	fi
	@if [ -d "/usr/local/grass86" ] && [ "$(MODULE_TOPDIR)" != "/usr/local/grass86" ]; then \
		echo "Installing to GRASS 8.6..."; \
		cp -r $(ETCDIR) /usr/local/grass86/etc/; \
		$(INSTALL_DATA) $(PGM).html /usr/local/grass86/docs/html/$(PGM).html; \
		$(INSTALL_DATA) $(PGM).md /usr/local/grass86/docs/html/$(PGM).md; \
	fi

# Show configuration
show-config:
	@echo "GRASS Module: $(PGM)"
	@echo "GRASS TopDir: $(MODULE_TOPDIR)"
	@echo "Install Dir: $(INST_DIR)"
	@echo "Library Dir: $(ETCDIR)"
	@echo "Modules: $(MODULES)"
	@echo "Python Files: $(PYFILES)"

# Test installation
test: install
	@echo "Testing module installation..."
	@if [ -f "$(SCRIPTDIR)/$(PGM)" ]; then \
		echo "✅ Script installed: $(SCRIPTDIR)/$(PGM)"; \
	else \
		echo "❌ Script not found: $(SCRIPTDIR)/$(PGM)"; \
	fi
	@if [ -d "$(INST_DIR)/etc/i_hyper_lib" ]; then \
		echo "✅ Library installed: $(INST_DIR)/etc/i_hyper_lib"; \
		ls -la $(INST_DIR)/etc/i_hyper_lib/*.py | wc -l | xargs echo "  Python modules:"; \
	else \
		echo "❌ Library not found: $(INST_DIR)/etc/i_hyper_lib"; \
	fi
	@if [ -f "$(HTMLDIR)/$(PGM).html" ]; then \
		echo "✅ HTML documentation: $(HTMLDIR)/$(PGM).html"; \
	else \
		echo "❌ HTML docs not found: $(HTMLDIR)/$(PGM).html"; \
	fi
	@if [ -f "$(HTMLDIR)/$(PGM).md" ]; then \
		echo "✅ MD documentation: $(HTMLDIR)/$(PGM).md"; \
	else \
		echo "❌ MD docs not found: $(HTMLDIR)/$(PGM).md"; \
	fi

.PHONY: default install clean uninstall install-all show-config test
