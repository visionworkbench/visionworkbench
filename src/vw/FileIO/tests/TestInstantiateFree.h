#include <vw/FileIO.h>

using namespace vw;

typedef PixelRGBA<int> PixelT;
typedef char CharT;
typedef char _CharT;
typedef ImageView<PixelT> ImageT;
typedef ImageT ElemT;

#include "TestInstantiateFreeList.hh"

class TestInstantiateFileIOFree : public CxxTest::TestSuite
{
  public: void test_inst() {}
};
