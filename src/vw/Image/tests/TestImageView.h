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

// TestImageView.h
#include <cxxtest/TestSuite.h>

#include <vw/Image/ImageView.h>
#include <vw/Image/PixelTypes.h>

using namespace std;
using namespace vw;

class TestImageView : public CxxTest::TestSuite
{
public:

  void testDefaultConstructor()
  {
    ImageView<double> test_double;
    TS_ASSERT_EQUALS(test_double.cols(), 0);
    TS_ASSERT_EQUALS(test_double.rows(), 0);
    TS_ASSERT_EQUALS(test_double.planes(), 0);
    TS_ASSERT_EQUALS(test_double.data(), (double*)0);

    ImageView<PixelRGBA<uint8> > test_rgba;
    TS_ASSERT_EQUALS(test_rgba.cols(), 0);
    TS_ASSERT_EQUALS(test_rgba.rows(), 0);
    TS_ASSERT_EQUALS(test_rgba.planes(), 0);
    TS_ASSERT_EQUALS(test_rgba.data(), (PixelRGBA<uint8>*)0);
  }

  void testColsRowsConstructor()
  {
    ImageView<double> test_double(3,4);
    TS_ASSERT_EQUALS(test_double.cols(), 3);
    TS_ASSERT_EQUALS(test_double.rows(), 4);
    TS_ASSERT_EQUALS(test_double.planes(), 1);
    TS_ASSERT_DIFFERS(test_double.data(), (double*)0);

    ImageView<PixelRGBA<uint8> > test_rgba(3,4);
    TS_ASSERT_EQUALS(test_rgba.cols(), 3);
    TS_ASSERT_EQUALS(test_rgba.rows(), 4);
    TS_ASSERT_EQUALS(test_rgba.planes(), 1);
    TS_ASSERT_DIFFERS(test_rgba.data(), (PixelRGBA<uint8>*)0);
  }

  void testColsRowsPlanesConstructor()
  {
    ImageView<double> test_double(4,3,2);
    TS_ASSERT_EQUALS(test_double.cols(), 4);
    TS_ASSERT_EQUALS(test_double.rows(), 3);
    TS_ASSERT_EQUALS(test_double.planes(), 2);
    TS_ASSERT_DIFFERS(test_double.data(), (double*)0);

    ImageView<PixelRGBA<uint8> > test_rgba(4,3,2);
    TS_ASSERT_EQUALS(test_rgba.cols(), 4);
    TS_ASSERT_EQUALS(test_rgba.rows(), 3);
    TS_ASSERT_EQUALS(test_rgba.planes(), 2);
    TS_ASSERT_DIFFERS(test_rgba.data(), (PixelRGBA<uint8>*)0);
  }

  void testCopyConstructor()
  {
    ImageView<double> test_double(3,4);
    ImageView<double> test2_double( test_double );
    TS_ASSERT_EQUALS(test2_double.cols(), 3);
    TS_ASSERT_EQUALS(test2_double.rows(), 4);
    TS_ASSERT_EQUALS(test2_double.planes(), 1);
    TS_ASSERT_EQUALS(test2_double.data(),test_double.data());

    ImageView<PixelRGBA<uint8> > test_rgba(3,4);
    ImageView<PixelRGBA<uint8> > test2_rgba( test_rgba );
    TS_ASSERT_EQUALS(test2_rgba.cols(), 3);
    TS_ASSERT_EQUALS(test2_rgba.rows(), 4);
    TS_ASSERT_EQUALS(test2_rgba.planes(), 1);
    TS_ASSERT_EQUALS(test2_rgba.data(), test_rgba.data());
  }

  void testSetSize()
  {
    ImageView<double> test_double(3,4);
    test_double.set_size(2,3);
    TS_ASSERT_EQUALS(test_double.cols(), 2);
    TS_ASSERT_EQUALS(test_double.rows(), 3);
    TS_ASSERT_EQUALS(test_double.planes(), 1);
    test_double.set_size(0,0);
    TS_ASSERT_EQUALS(test_double.cols(), 0);
    TS_ASSERT_EQUALS(test_double.rows(), 0);
    TS_ASSERT_EQUALS(test_double.data(), (double*)0);

    ImageView<PixelRGBA<uint8> > test_rgba(3,4);
    test_rgba.set_size(2,3);
    TS_ASSERT_EQUALS(test_rgba.cols(), 2);
    TS_ASSERT_EQUALS(test_rgba.rows(), 3);
    TS_ASSERT_EQUALS(test_rgba.planes(), 1);
    test_rgba.set_size(0,0);
    TS_ASSERT_EQUALS(test_rgba.cols(), 0);
    TS_ASSERT_EQUALS(test_rgba.rows(), 0);
    TS_ASSERT_EQUALS(test_rgba.data(), (PixelRGBA<uint8>*)0);
  }

  void testReset()
  {
    ImageView<double> test_double(3,4);
    test_double.reset();
    TS_ASSERT_EQUALS(test_double.cols(), 0);
    TS_ASSERT_EQUALS(test_double.rows(), 0);
    TS_ASSERT_EQUALS(test_double.planes(), 0);
    TS_ASSERT_EQUALS(test_double.data(), (double*)0);

    ImageView<PixelRGBA<uint8> > test_rgba(3,4);
    test_rgba.reset();
    TS_ASSERT_EQUALS(test_rgba.cols(), 0);
    TS_ASSERT_EQUALS(test_rgba.rows(), 0);
    TS_ASSERT_EQUALS(test_rgba.planes(), 0);
    TS_ASSERT_EQUALS(test_rgba.data(), (PixelRGBA<uint8>*)0);
  }

};


