# __BEGIN_LICENSE__
# 
# Copyright (C) 2006 United States Government as represented by the
# Administrator of the National Aeronautics and Space Administration
# (NASA).  All Rights Reserved.
# 
# Copyright 2006 Carnegie Mellon University. All rights reserved.
# 
# This software is distributed under the NASA Open Source Agreement
# (NOSA), version 1.3.  The NOSA has been approved by the Open Source
# Initiative.  See the file COPYING at the top of the distribution
# directory tree for the complete NOSA document.
# 
# THE SUBJECT SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY WARRANTY OF ANY
# KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT
# LIMITED TO, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL CONFORM TO
# SPECIFICATIONS, ANY IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR
# A PARTICULAR PURPOSE, OR FREEDOM FROM INFRINGEMENT, ANY WARRANTY THAT
# THE SUBJECT SOFTWARE WILL BE ERROR FREE, OR ANY WARRANTY THAT
# DOCUMENTATION, IF PROVIDED, WILL CONFORM TO THE SUBJECT SOFTWARE.
# 
# __END_LICENSE__

########################################################################
# sources
########################################################################
ACLOCAL_AMFLAGS = -I m4 -I thirdparty/m4

SUBDIRS = src

EXTRA_DIST = config/rules.mak config.options.example NOTES docs \
			 thirdparty/cxxtest INSTALLGUIDE scripts/include-graph.py \
			 scripts/relicense.pl scripts/static-checking.js scripts/templates.js
########################################################################
# special top-level rules
########################################################################

install-data-hook:
	mkdir -p $(DESTDIR)$(prefix)/include/vw/
	$(INSTALL) src/vw/config.h $(DESTDIR)$(prefix)/include/vw/

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

lcov:
	$(LCOV) --capture --test-name 'Vision-Workbench-Base' --directory $(top_srcdir) --output-file vw-base.trace --initial
	$(LCOV) --capture --test-name 'Vision-Workbench'      --directory $(top_srcdir) --output-file vw-test.trace
	$(LCOV) -a vw-base.trace -a vw-test.trace -o vw.trace
	$(LCOV) -e vw.trace '*/src/vw/*' -o vw.trace
	$(LCOV) -r vw.trace '*/tests/*'  -o vw.trace
	$(GENHTML) -o $(top_srcdir)/coverage vw.trace
	@echo "Generated code coverage report in $(top_srcdir)/coverage"

push-lcov:
	if test -n "$(LCOV_REMOTE_PATH)" ; then rsync -P -e ssh -r coverage/ $(LCOV_REMOTE_PATH); fi

DEHYDRA_MODULES = \
  $(top_srcdir)/scripts/templates.js

DEHYDRA_ARGS = \
  $(top_srcdir)/scripts/static-checking.js \
  --output-path=$(top_srcdir)/src/vw \
  --dehydra-modules=$(DEHYDRA_MODULES)

instantiate-all:
	echo "#include <vw/vw.h>" | \
		$(DEHYDRA_GCC) -x c++ - -S -o /dev/null -fplugin=$(DEHYDRA_PLUGIN) -fplugin-arg='$(DEHYDRA_ARGS)' $(VW_CPPFLAGS)

# This is a list of every header (eventually) included when you include vw/vw.h
all-headers:
	echo "#include <vw/vw.h>" | \
	$(CXX) -x c++ -E - $(VW_CPPFLAGS) | \
	grep '^#.*vw/' | sed -e 's#.*src/\(vw/[^"]\+\)".*#\1#' | sort -u


########################################################################
# general
########################################################################
CLEANFILES = vw.trace vw-base.trace vw-test.trace

include $(top_srcdir)/config/rules.mak
