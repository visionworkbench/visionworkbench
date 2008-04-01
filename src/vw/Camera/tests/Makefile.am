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

if MAKE_MODULE_CAMERA

TESTS = TestCAHVORModel TestLinearPushbroomModel TestPinholeModel TestBundleAdjust TestExifData # TestLensDistortion TestCameraTransform

nodist_TestCAHVORModel_SOURCES = TestCAHVORModel.cxx
nodist_TestLinearPushbroomModel_SOURCES = TestLinearPushbroomModel.cxx
nodist_TestPinholeModel_SOURCES = TestPinholeModel.cxx
nodist_TestBundleAdjust_SOURCES = TestBundleAdjust.cxx
nodist_TestExifData_SOURCES = TestExifData.cxx
#nodist_TestLensDistortion_SOURCES = TestLensDistortion.cxx
#nodist_TestCameraTransform_SOURCES = TestCameraTransform.cxx

endif

########################################################################
# general
########################################################################

AM_CPPFLAGS = @VW_CPPFLAGS@ -I$(top_srcdir)/thirdparty/cxxtest
AM_LDFLAGS  = @VW_LDFLAGS@ @PKG_CAMERA_LIBS@ @PKG_VW_LIBS@

check_PROGRAMS = $(TESTS)
EXTRA_DIST = $(patsubst %,%.h,$(TESTS))
CLEANFILES = $(patsubst %,%.cxx,$(TESTS))

include $(top_srcdir)/config/rules.mak
