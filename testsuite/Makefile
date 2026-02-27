MODULE_TOPDIR = ../..

include $(MODULE_TOPDIR)/include/Make/Dir.make

default: parsertest

TESTS = test_i_hyper_smac.py

parsertest:
	@for test in $(TESTS) ; do \
		echo "Running $$test..." ; \
		$(MODULE_TOPDIR)/tools/gunittest/testrunner.py $$test ; \
	done

.PHONY: default parsertest
