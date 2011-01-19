# __BEGIN_LICENSE__
# Copyright (C) 2006-2010 United States Government as represented by
# the Administrator of the National Aeronautics and Space Administration.
# All Rights Reserved.
# __END_LICENSE__


########################################################################
# sources
########################################################################
ACLOCAL_AMFLAGS = -I m4 -I thirdparty/m4
DISTCHECK_CONFIGURE_FLAGS = --enable-ccache
LCOV_ARGS = --ignore-errors graph,gcov --follow --compat-libtool

SUBDIRS = src

EXTRA_DIST = \
  autogen \
  config.options.example \
  config/rules.mak \
  docs \
  INSTALLGUIDE \
  README.md \
  scripts/include-graph.py \
  scripts/.lcovrc \
  scripts/python/vw/gdb/__init__.py \
  scripts/python/vw/gdb/types.py \
  scripts/python/vw/__init__.py \
  scripts/relicense.pl \
  thirdparty/cxxtest \
  thirdparty/gtest/CHANGES \
  thirdparty/gtest/CONTRIBUTORS \
  thirdparty/gtest/COPYING \
  thirdparty/gtest/include/gtest/gtest.h \
  thirdparty/gtest/README \
  thirdparty/gtest/src/gtest-all.cc \
  thirdparty/gtest/src/gtest_main.cc


########################################################################
# special top-level rules
########################################################################

dist-hook:
	rm -rf `find $(distdir) -name '.svn' -o -name '.git'`
	make -C $(distdir)/docs/workbook clean
if DIST_LICENSE
	cp $(DIST_LICENSE) $(distdir)/COPYING
endif
if DIST_LICENSE_SUMMARY
	find $(distdir) -type f | scripts/relicense.pl $(DIST_LICENSE_SUMMARY)
endif
if DIST_CONFIG_OPTIONS_DEFAULT
	cp $(DIST_CONFIG_OPTIONS_DEFAULT) $(distdir)/config/config.options.default
endif
if DIST_REMOVE_MODULES
	for module in $(DIST_REMOVE_MODULES) ; do \
	  rm -rf $(distdir)/src/vw/$${module} $(distdir)/src/vw/$${module}.h ; \
	  grep -v "$${module}.*Makefile" $(distdir)/configure.ac > configure.new ; \
	  rebuild=0 ; \
	  if ! diff -q $(distdir)/configure.ac configure.new ; then \
		  mv configure.new $(distdir)/configure.ac ; \
		  rebuild=1 ; \
	  fi ; \
	  grep -v "= $${module}" $(distdir)/src/vw/Makefile.am > Makefile.new ; \
	  if ! diff -q $(distdir)/src/vw/Makefile.am Makefile.new ; then \
		  mv Makefile.new $(distdir)/src/vw/Makefile.am ; \
		  rebuild=1 ; \
	  fi ; \
	  rm -f configure.new Makefile.new ; \
	done ; \
	if test x"$$rebuild" = "x1"; then \
		( cd $(distdir) && autoreconf && touch src/vw/config.h.pre.in && rm -rf autom4te.cache ) ; \
	fi
endif

check-local:
	@echo ==============================
	@echo All tests passed successfully!
	@echo ==============================

.gcno.stamp: .FORCE
	$(MAKE) $(AM_MAKEFLAGS) STAMPFILE=$@ STAMPPATTERN='*.gcno' copy-timestamp

.gcda.stamp: .FORCE
	$(MAKE) $(AM_MAKEFLAGS) STAMPFILE=$@ STAMPPATTERN='*.gcda' copy-timestamp

#This strange HOME bit is so lcov finds the config file we have in the scripts dir.
vw-base.trace: .gcno.stamp
	rm -f $@
	HOME=$(top_srcdir)/scripts $(LCOV) $(LCOV_ARGS) --capture --test-name 'Vision-Workbench' --directory $(top_srcdir) --output-file vw-base.trace --initial
	$(MAKE) $(AM_MAKEFLAGS) TRACEFILE=vw-base.trace coverage-filter

vw-test.trace: .gcda.stamp
	rm -f $@
	HOME=$(top_srcdir)/scripts $(LCOV) $(LCOV_ARGS) --capture --test-name 'Vision-Workbench' --directory $(top_srcdir) --output-file vw-test.trace
	$(MAKE) $(AM_MAKEFLAGS) TRACEFILE=vw-test.trace coverage-filter

vw.trace: vw-base.trace vw-test.trace
	rm -f $@
	HOME=$(top_srcdir)/scripts $(LCOV) $(LCOV_ARGS) -a vw-base.trace -a vw-test.trace -o vw.trace

$(top_srcdir)/coverage/index.html: vw.trace
	rm -rf $(top_srcdir)/coverage
	HOME=$(top_srcdir)/scripts $(GENHTML) --demangle-cpp -o $(top_srcdir)/coverage vw.trace
	@echo "Generated code coverage report in $(top_srcdir)/coverage"

coverage: $(top_srcdir)/coverage/index.html

# The tracefiles can get really big. This removes unnecessary data from them.
coverage-filter:
	rm -f $(TRACEFILE).tmp
	HOME=$(top_srcdir)/scripts $(LCOV) $(LCOV_ARGS) -e $(TRACEFILE)     '*/src/vw/*'              -o $(TRACEFILE).tmp
	HOME=$(top_srcdir)/scripts $(LCOV) $(LCOV_ARGS) -r $(TRACEFILE).tmp '*/tests/*'               -o $(TRACEFILE).tmp
	HOME=$(top_srcdir)/scripts $(LCOV) $(LCOV_ARGS) -r $(TRACEFILE).tmp '*/tools/*'               -o $(TRACEFILE).tmp
	HOME=$(top_srcdir)/scripts $(LCOV) $(LCOV_ARGS) -r $(TRACEFILE).tmp '*/gui/*'                 -o $(TRACEFILE).tmp
	HOME=$(top_srcdir)/scripts $(LCOV) $(LCOV_ARGS) -r $(TRACEFILE).tmp '*/src/vw/Plate/google/*' -o $(TRACEFILE).tmp
	mv $(TRACEFILE).tmp $(TRACEFILE)
	rm -f $(TRACEFILE).tmp

# This copies the timestamp of the newest file that matches STAMPPATERN to the
# timestamp of STAMPFILE
copy-timestamp:
	@test -f $(STAMPFILE) || touch -t 197001010001 $(STAMPFILE)
	@find $(top_srcdir)/src/vw \( -name '$(STAMPPATTERN)' -a -newer $(STAMPFILE) \) \
		-a -exec test {} -nt $(STAMPFILE) \; \
	    -a -exec touch -r {} $(STAMPFILE) \;

distclean-local:
	cd $(srcdir)/docs/workbook && $(MAKE) $(AM_MAKEFLAGS) distclean

workbook: all
	cd $(srcdir)/docs/workbook && $(MAKE) $(AM_MAKEFLAGS) workbook.pdf

doxygen:
	cd src && doxygen

doc: workbook doxygen

install-with-doc: doc install
	mkdir -p $(DESTDIR)usr/share/doc/visionworkbench
	cp docs/workbook/workbook.pdf $(DESTDIR)usr/share/doc/visionworkbench/
	cp -r docs/doxygen $(DESTDIR)usr/share/doc/visionworkbench/apiReference

push-coverage:
	if test -n "$(LCOV_REMOTE_PATH)" ; then rsync --delete -P -e ssh -r coverage/ $(LCOV_REMOTE_PATH); fi

DEHYDRA_MODULES = \
  $(top_srcdir)/scripts/templates.js

DEHYDRA_ARGS = \
  $(top_srcdir)/thirdparty/static-checking.js \
  --output-path=$(top_srcdir)/src/vw \
  --no-backup \
  --dehydra-modules=$(DEHYDRA_MODULES)

# Generate instantiation files for all templates
instantiate-all:
	echo "#include <vw/vw.h>" | \
		$(DEHYDRA_GCC) -x c++ - -S -o /dev/null -fplugin=$(DEHYDRA_PLUGIN) -fplugin-arg='$(DEHYDRA_ARGS)' $(VW_CPPFLAGS)

# List the instantiation percentages
instantiate-skipped:
	@find -type f -name TestInstantiate\*List.hh -print0 | sort -z | xargs -0 -n1 awk '/^\/\//{n++}; END { printf "%5.1f%% of %4i lines skipped: %s\n", ((100*n)/NR), NR, FILENAME }'

# This is a list of every header (eventually) included when you include vw/vw.h
all-headers:
	echo "#include <vw/vw.h>" | \
	$(CXX) -x c++ -E - $(VW_CPPFLAGS) | \
	grep '^#.*vw/' | sed -e 's#.*src/\(vw/[^"]\+\)".*#\1#' | sort -u

one-ctags-file:
	ctags -R \
		--regex-c++='/VW_DEFINE_EXCEPTION\(\s*(.*)\s*,\s*(.*)\s*\)/\1/s/' \
		--regex-c++='/VW_DEFINE_EXCEPTION_EXT\(\s*(.*)\s*,\s*(.*)\s*\)/\1/s/' \
		--regex-c++='/__VW_UNARY_MATH_FUNCTOR\(\s*([A-Z].*)\s*,\s*(.*)\s*\)/Arg\1Functor/s/' \
		--regex-c++='/__VW_BINARY_MATH_FUNCTOR\(\s*([A-Z].*)\s*,\s*(.*)\s*\)/ArgArg\1Functor/s/' \
		--regex-c++='/VW_DEFINE_ENUM\(\s*(.*)\s*,\s*(.*)\s*,\s*\(\s*(.*)\s*\)\)/\1/g/' \
		--languages=C,C++ --exclude='*thirdparty/*' --exclude='*tests/*' --exclude='*vw/GPU/*' . --exclude='' thirdparty/cxxtest/cxxtest/TestSuite.h thirdparty/gtest

.PHONY: push-coverage instantiate-all instantiate-skipped all-headers one-ctags-file copy-timestamp coverage-filter

# This rule is like .PHONY except it doesn't apply recursively to rules that depend on phony rules
.FORCE:

########################################################################
# general
########################################################################
CLEANFILES = vw.trace vw-base.trace vw-test.trace

include $(top_srcdir)/config/rules.mak
