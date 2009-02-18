# __BEGIN_LICENSE__
# Copyright (C) 2006, 2007 United States Government as represented by
# the Administrator of the National Aeronautics and Space Administration.
# All Rights Reserved.
# __END_LICENSE__


########################################################################
# tests (using cxxtest)
########################################################################

SUFFIXES = .cxx

CXXTEST_DIR =   $(top_srcdir)/thirdparty/cxxtest
CXXTEST_GEN =   $(CXXTEST_DIR)/cxxtestgen.pl
CXXTEST_ARGS =  --error-printer

.h.cxx:
	$(CXXTEST_GEN) $(CXXTEST_ARGS) -o $@ $<

newtest:
	@if test -z "$(NAME)"; then echo "run make NAME=TestName [MODULE=ModuleName] newtest"; else $(top_srcdir)/scripts/create-test.sh $(NAME) $(MODULE); fi

.PHONY: newtest
