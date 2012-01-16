// __BEGIN_LICENSE__
// Copyright (C) 2006-2011 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <gtest/gtest_VW.h>
#include <test/Helpers.h>
#include <vw/Image/ImageView.h>

using namespace vw;
using namespace vw::test;
using namespace std;

// spike operator new... every nothrow new fails [this lets us ask for
// unlimited memory without the system actually trying to give it]
void* operator new[](size_t, const nothrow_t&) throw() {
  return 0;
}

void operator delete[](void*, const nothrow_t&) throw() {
  /* noop */
}

TEST(ImageView, TooMuchMemory) {
  int32 cols, rows, planes;
  cols = rows = 1<<15;
  planes = 1 << 9;

  ImageView<uint8> img;

  // We can probably only trip this error on 32-bit platforms
  if (sizeof(size_t) <= sizeof(int32))
    EXPECT_THROW_MSG(img.set_size(cols, rows, planes), ArgumentErr, "too many pixels");

  rows = planes = 1;
  EXPECT_THROW_MSG(img.set_size(cols, rows, planes), ArgumentErr, "too many bytes");
}
