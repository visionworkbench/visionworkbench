# __BEGIN_LICENSE__
#
# Copyright (C) 2006 United States Government as represented by the
# Administrator of the National Aeronautics and Space Administration
# (NASA).  All Rights Reserved.
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

if MAKE_MODULE_FILEIO

TESTS = TestDiskImageResource

nodist_TestDiskImageResource_SOURCES = TestDiskImageResource.cxx

endif

########################################################################
# general
########################################################################

AM_CPPFLAGS = @MODULE_FILEIO_CPPFLAGS@ -I$(top_srcdir)/thirdparty/cxxtest
AM_LDFLAGS  = @PKG_FILEIO_LDFLAGS@

check_PROGRAMS = $(TESTS)

EXTRA_DIST = $(patsubst %,%.h,$(TESTS)) rgb2x2.png rgb2x2.tif rgb2x2.jpg

CLEANFILES = $(patsubst %,%.cxx,$(TESTS)) tmp.png tmp.tif

include $(top_srcdir)/config/rules.mak
