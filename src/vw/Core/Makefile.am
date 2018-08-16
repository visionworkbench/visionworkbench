# __BEGIN_LICENSE__
#  Copyright (c) 2006-2013, United States Government as represented by the
#  Administrator of the National Aeronautics and Space Administration. All
#  rights reserved.
#
#  The NASA Vision Workbench is licensed under the Apache License,
#  Version 2.0 (the "License"); you may not use this file except in
#  compliance with the License. You may obtain a copy of the License at
#  http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.
# __END_LICENSE__


########################################################################
# sources
########################################################################

if MAKE_MODULE_CORE

include_HEADERS = \
  Cache.h Cache.tcc \
  CompoundTypes.h \
  Condition.h \
  ConfigParser.h \
  Debugging.h \
  Exception.h \
  Features.h \
  Functors.h \
  FundamentalTypes.h \
  Log.h \
  ProgressCallback.h \
  RunOnce.h \
  Settings.h \
  StringUtils.h \
  Stopwatch.h \
  System.h \
  Thread.h \
  ThreadPool.h \
  ThreadQueue.h \
  TypeDeduction.h \
  VarArray.h \
  CmdUtils.h

libvwCore_la_SOURCES = \
  Cache.cc \
  ConfigParser.cc \
  Debugging.cc \
  Exception.cc \
  Log.cc \
  ProgressCallback.cc \
  Settings.cc \
  StringUtils.cc \
  Stopwatch.cc \
  System.cc \
  Thread.cc \
  ThreadPool.cc \
  CmdUtils.cc

libvwCore_la_LIBADD = @MODULE_CORE_LIBS@

lib_LTLIBRARIES = libvwCore.la

endif

########################################################################
# general
########################################################################

AM_CPPFLAGS = @VW_CPPFLAGS@ @PKG_BOOST_CPPFLAGS@
AM_LDFLAGS = @VW_LDFLAGS@ $(LIBTOOL_LIB_FLAGS)

SUBDIRS = . tests

includedir = $(prefix)/include/vw/Core

include $(top_srcdir)/config/rules.mak
