# __BEGIN_LICENSE__
# Copyright (C) 2006, 2007 United States Government as represented by
# the Administrator of the National Aeronautics and Space Administration.
# All Rights Reserved.
# __END_LICENSE__


########################################################################
# sources
########################################################################
ACLOCAL_AMFLAGS = -I m4 -I thirdparty/m4

SUBDIRS = src

EXTRA_DIST = config/rules.mak config.options.example NOTES docs \
			 thirdparty/cxxtest INSTALLGUIDE scripts/include-graph.py \
			 scripts/relicense.pl scripts/templates.js thirdparty/m4/.keep \
			 scripts/create-test.sh
########################################################################
# special top-level rules
########################################################################

dist-hook:
	rm -rf `find $(distdir) -name '.svn'`
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
	  mv configure.new $(distdir)/configure.ac ; \
	  grep -v "= $${module}" $(distdir)/src/vw/Makefile.am > Makefile.new ; \
	  mv Makefile.new $(distdir)/src/vw/Makefile.am ; \
	done
	( cd $(distdir) && autoreconf && touch src/vw/config.h.pre.in && rm -rf autom4te.cache )
endif

check-local:
	@echo ==============================
	@echo All tests passed successfully!
	@echo ==============================

distclean-local:
	cd docs/workbook && $(MAKE) $(AM_MAKEFLAGS) distclean

workbook: all
	cd docs/workbook && $(MAKE) $(AM_MAKEFLAGS) workbook.pdf

doxygen: 
	cd src && doxygen

doc: workbook doxygen

coverage:
	$(LCOV) --capture --test-name 'Vision-Workbench-Base' --directory $(top_srcdir) --output-file vw-base.trace --initial
	$(LCOV) -e vw-base.trace '*/src/vw/*' -o vw-base.trace
	$(LCOV) -r vw-base.trace '*/tests/*'  -o vw-base.trace
	$(LCOV) -r vw-base.trace '*/tools/*'  -o vw-base.trace
	$(LCOV) -r vw-base.trace '*/gui/*'  -o vw-base.trace
	$(LCOV) --capture --test-name 'Vision-Workbench'      --directory $(top_srcdir) --output-file vw-test.trace
	$(LCOV) -e vw-test.trace '*/src/vw/*' -o vw-test.trace
	$(LCOV) -r vw-test.trace '*/tests/*'  -o vw-test.trace
	$(LCOV) -r vw-test.trace '*/tools/*'  -o vw-test.trace
	$(LCOV) -r vw-test.trace '*/gui/*'  -o vw-test.trace
	$(LCOV) -a vw-base.trace -a vw-test.trace -o vw.trace
	$(GENHTML) -o $(top_srcdir)/coverage vw.trace
	@echo "Generated code coverage report in $(top_srcdir)/coverage"

push-coverage:
	if test -n "$(LCOV_REMOTE_PATH)" ; then rsync -P -e ssh -r coverage/ $(LCOV_REMOTE_PATH); fi

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
	ctags -R --languages=C,C++ --exclude='*thirdparty/*' --exclude='*tests/*'

.PHONY: coverage push-coverage instantiate-all instantiate-skipped all-headers one-ctags-file

########################################################################
# general
########################################################################
CLEANFILES = vw.trace vw-base.trace vw-test.trace

include $(top_srcdir)/config/rules.mak
