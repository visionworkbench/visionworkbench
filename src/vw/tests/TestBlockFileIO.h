// __BEGIN_LICENSE__
//
// Copyright (C) 2006 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration
// (NASA).  All Rights Reserved.
// 
// This software is distributed under the NASA Open Source Agreement
// (NOSA), version 1.3.  The NOSA has been approved by the Open Source
// Initiative.  See the file COPYING at the top of the distribution
// directory tree for the complete NOSA document.
// 
// THE SUBJECT SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY WARRANTY OF ANY
// KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT
// LIMITED TO, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL CONFORM TO
// SPECIFICATIONS, ANY IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR
// A PARTICULAR PURPOSE, OR FREEDOM FROM INFRINGEMENT, ANY WARRANTY THAT
// THE SUBJECT SOFTWARE WILL BE ERROR FREE, OR ANY WARRANTY THAT
// DOCUMENTATION, IF PROVIDED, WILL CONFORM TO THE SUBJECT SOFTWARE.
//
// __END_LICENSE__

// TestBlockFileIO.h
#define CXXTEST_ABORT_TEST_ON_FAIL
#include <cxxtest/TestSuite.h>
#include <vw/FileIO/DiskImageResource.h>

#include <vw/Image.h>
#include <vw/FileIO.h>

using namespace vw;

class TestDiskImageResource : public CxxTest::TestSuite
{
public:

  void test_read_image() {
    std::cout << std::endl;
    ImageView<PixelRGB<uint8> > image;
    DiskImageResource *dim;
    TS_ASSERT_THROWS_NOTHING( dim = DiskImageResource::open( "mural.png" ) );
    std::cout << "Block size: " << dim->native_read_block_size() << std::endl;
    TS_ASSERT_THROWS_NOTHING( dim->read(image, BBox2i(100,100,100,100) ) );
    std::cout << "Image: " << image.cols() << "x" << image.rows() << std::endl;
    TS_ASSERT_THROWS_NOTHING( write_image( "mural.tif", image ) );
  }

}; // class TestBlockFileIO
