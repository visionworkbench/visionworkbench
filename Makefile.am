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

SUBDIRS = src

EXTRA_DIST = config/rules.mak config/relicense config.options.example \
	NOTES docs thirdparty 

########################################################################
# special top-level rules
########################################################################

install-data-hook:
	$(INSTALL) src/vw/config.h $(prefix)/include/vw/

dist-hook:
	rm -rf `find $(distdir) -name '.svn'`
if DIST_LICENSE
	cp $(DIST_LICENSE) $(distdir)/COPYING
endif
if DIST_LICENSE_SUMMARY
	find $(distdir) -type f | config/relicense $(DIST_LICENSE_SUMMARY) -
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
	cd src/vw && doxygen

########################################################################
# general
########################################################################

include $(top_srcdir)/config/rules.mak
