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
MODULES = smac aod wvc radtran o3 utils smac_coef_generator lut spectral_polish adjacency gas_absorption parallel_lut smart_lut opencl_accelerator
ETCDIR = $(ETC)/i_hyper_lib
COEFDIR = $(ETCDIR)/COEFS
PYFILES := $(patsubst %,$(ETCDIR)/%.py,$(MODULES))
MODULE_SRCS := $(patsubst %,lib/%.py,$(MODULES))

# Default target - install modules first, then generate docs
default: $(PYFILES) install_coefs script html docs

# Documentation targets
docs: html-docs programmer-docs

# Generate programmer documentation using Doxygen
programmer-docs: $(HTMLDIR)/$(PGM)_programmer.html

# Ensure programmer manual HTML is installed
$(HTMLDIR)/$(PGM)_programmer.html: i_hyper_smac.dox
	@echo "Generating programmer documentation..."
	@if command -v doxygen >/dev/null 2>&1; then \
		doxygen i_hyper_smac.dox; \
		if [ -d "html" ]; then \
			cp -r html/* $(HTMLDIR)/; \
			echo "✓ Programmer documentation installed to $(HTMLDIR)"; \
		else \
			echo "⚠ Doxygen documentation generation completed, but no output found"; \
		fi; \
	else \
		echo "⚠ Doxygen not found. Installing programmer manual markdown..."; \
		$(INSTALL_DATA) PROGRAMMER_MANUAL.md $(HTMLDIR)/$(PGM)_programmer.md; \
	fi

# Ensure programmer manual MD is installed
$(HTMLDIR)/$(PGM)_programmer.md: PROGRAMMER_MANUAL.md
	$(INSTALL_DATA) PROGRAMMER_MANUAL.md $(HTMLDIR)/$(PGM)_programmer.md

# Generate both HTML and markdown programmer docs
programmer-docs-full: $(HTMLDIR)/$(PGM)_programmer.html $(HTMLDIR)/$(PGM)_programmer.md
	@echo "✓ Programmer documentation (HTML + MD) generated"

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
install: script html $(PYFILES) install_coefs programmer-docs
	# Copy library to installation directory
	cp -r $(ETCDIR) $(INST_DIR)/etc/
	# Install documentation
	$(INSTALL_DATA) $(PGM).html $(HTMLDIR)/$(PGM).html
	$(INSTALL_DATA) $(PGM).md $(HTMLDIR)/$(PGM).md
	# Install programmer manual
	@if [ -f "$(HTMLDIR)/$(PGM)_programmer.html" ]; then \
		echo "✓ Programmer HTML documentation installed"; \
	fi
	@if [ -f "$(HTMLDIR)/$(PGM)_programmer.md" ]; then \
		echo "✓ Programmer MD documentation installed"; \
	fi
	# Install Python script without .py extension (GRASS policy)
	$(INSTALL) $(PGM).py $(SCRIPTDIR)/$(PGM)

# Clean build artifacts
clean:
	rm -rf $(ETCDIR)
	rm -f $(HTMLDIR)/$(PGM).html
	rm -f $(HTMLDIR)/$(PGM).md
	rm -f $(HTMLDIR)/$(PGM)_programmer.html
	rm -f $(HTMLDIR)/$(PGM)_programmer.md
	rm -rf html docs/html docs/latex
	@echo "✓ Build artifacts cleaned"

# Uninstall from system
uninstall:
	rm -rf $(INST_DIR)/etc/i_hyper_lib
	rm -f $(HTMLDIR)/$(PGM).html
	rm -f $(HTMLDIR)/$(PGM).md
	rm -f $(HTMLDIR)/$(PGM)_programmer.html
	rm -f $(HTMLDIR)/$(PGM)_programmer.md
	rm -f $(BIN)/$(PGM)
	@echo "✓ Module uninstalled"

# Install to multiple GRASS versions
install-all: install
	@if [ -d "/usr/local/grass85" ] && [ "$(MODULE_TOPDIR)" != "/usr/local/grass85" ]; then \
		echo "Installing to GRASS 8.5..."; \
		cp -r $(ETCDIR) /usr/local/grass85/etc/; \
		$(INSTALL_DATA) $(PGM).html /usr/local/grass85/docs/html/$(PGM).html; \
		$(INSTALL_DATA) $(PGM).md /usr/local/grass85/docs/html/$(PGM).md; \
		if [ -f "$(HTMLDIR)/$(PGM)_programmer.html" ]; then \
			$(INSTALL_DATA) $(HTMLDIR)/$(PGM)_programmer.html /usr/local/grass85/docs/html/; \
		fi; \
		if [ -f "$(HTMLDIR)/$(PGM)_programmer.md" ]; then \
			$(INSTALL_DATA) $(HTMLDIR)/$(PGM)_programmer.md /usr/local/grass85/docs/html/; \
		fi; \
	fi
	@if [ -d "/usr/local/grass86" ] && [ "$(MODULE_TOPDIR)" != "/usr/local/grass86" ]; then \
		echo "Installing to GRASS 8.6..."; \
		cp -r $(ETCDIR) /usr/local/grass86/etc/; \
		$(INSTALL_DATA) $(PGM).html /usr/local/grass86/docs/html/$(PGM).html; \
		$(INSTALL_DATA) $(PGM).md /usr/local/grass86/docs/html/$(PGM).md; \
		if [ -f "$(HTMLDIR)/$(PGM)_programmer.html" ]; then \
			$(INSTALL_DATA) $(HTMLDIR)/$(PGM)_programmer.html /usr/local/grass86/docs/html/; \
		fi; \
		if [ -f "$(HTMLDIR)/$(PGM)_programmer.md" ]; then \
			$(INSTALL_DATA) $(HTMLDIR)/$(PGM)_programmer.md /usr/local/grass86/docs/html/; \
		fi; \
	fi

# Show configuration
show-config:
	@echo "GRASS Module: $(PGM)"
	@echo "GRASS TopDir: $(MODULE_TOPDIR)"
	@echo "Install Dir: $(INST_DIR)"
	@echo "Library Dir: $(ETCDIR)"
	@echo "HTML Dir: $(HTMLDIR)"
	@echo "Modules: $(MODULES)"
	@echo "Python Files: $(PYFILES)"
	@echo "Documentation:"
	@echo "  - User Manual: $(HTMLDIR)/$(PGM).html"
	@echo "  - User MD: $(HTMLDIR)/$(PGM).md"
	@echo "  - Programmer Manual: $(HTMLDIR)/$(PGM)_programmer.html"
	@echo "  - Programmer MD: $(HTMLDIR)/$(PGM)_programmer.md"

# Check documentation dependencies
doc-deps:
	@echo "Checking documentation dependencies..."
	@if command -v doxygen >/dev/null 2>&1; then \
		echo "✓ Doxygen found: $(shell doxygen --version)"; \
	else \
		echo "⚠ Doxygen not found - HTML programmer docs will not be generated"; \
		echo "  Install with: sudo apt-get install doxygen graphviz"; \
	fi
	@if [ -f "i_hyper_smac.dox" ]; then \
		echo "✓ Doxygen configuration found"; \
	else \
		echo "✗ Doxygen configuration not found"; \
	fi
	@if [ -f "PROGRAMMER_MANUAL.md" ]; then \
		echo "✓ Programmer manual markdown found"; \
	else \
		echo "✗ Programmer manual markdown not found"; \
	fi

# Generate documentation only
docs-only: doc-deps programmer-docs-full
	@echo "✓ Documentation generation complete"

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

# Run unit tests
test-unit: install
	@echo "Running unit tests..."
	@cd testsuite && python run_tests.py --unit-only

# Run integration tests
test-integration: install
	@echo "Running integration tests..."
	@cd testsuite && python run_tests.py --integration-only

# Run performance tests
test-performance: install
	@echo "Running performance tests..."
	@cd testsuite && python run_tests.py --performance-only

# Run all tests
test-all: install
	@echo "Running complete test suite..."
	@cd testsuite && python run_tests.py --verbose

# Generate test coverage
test-coverage: install
	@echo "Generating test coverage..."
	@cd testsuite && python -m coverage run run_tests.py && python -m coverage report

# Continuous integration test suite
test-ci: install
	@echo "Running CI test suite..."
	@cd testsuite && python run_tests.py --unit-only --integration-only

.PHONY: default install clean uninstall install-all show-config test test-unit test-integration test-performance test-all test-coverage test-ci docs html-docs programmer-docs programmer-docs-full doc-deps docs-only
