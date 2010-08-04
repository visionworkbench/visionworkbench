// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


// TestBlockFileIO.h
#include <gtest/gtest.h>
#include <string>
#include <vw/Core/Log.h>
#include <vw/FileIO/DiskImageResource_internal.h>

using namespace vw;
using namespace vw::internal;

static void test_extension(std::string const& fn)
{
  std::cout << vw::stringify("Testing ") + fn << "\n";
}

TEST( ExtTester, Test ) {
  foreach_ext("rwtest", test_extension);
}
