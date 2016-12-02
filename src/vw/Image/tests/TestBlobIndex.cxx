// __BEGIN_LICENSE__
//  Copyright (c) 2009-2013, United States Government as represented by the
//  Administrator of the National Aeronautics and Space Administration. All
//  rights reserved.
//
//  The NGT platform is licensed under the Apache License, Version 2.0 (the
//  "License"); you may not use this file except in compliance with the
//  License. You may obtain a copy of the License at
//  http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
// __END_LICENSE__


#include <test/Helpers.h>
#include <vw/Math/BBox.h>
#include <vw/Math/Vector.h>
#include <vw/Image/BlobIndex.h>
#include <vw/Image/Manipulation.h>
#include <vw/Image/MaskViews.h>
#include <vw/Image/PerPixelViews.h>
#include <vw/Image/PixelMask.h>
#include <vw/Image/PixelMath.h>
#include <vw/Image/Statistics.h>

// The BlobIndex tests are currently disabled.
// - To enable them a dependency needs to be added to FileIO.

//#include <vw/FileIO/DiskImageView.h>

#include <boost/assign/std/vector.hpp>
#include <boost/assign/list_of.hpp>

using namespace vw;
using namespace boost::assign;



TEST( BlobIndex, BlobIndex ) {
  typedef PixelMask<uint8> MPx;
  ImageView<MPx> im(7,5);
  fill( crop(im,1,1,3,1), MPx(255) ); // Blob 1
  fill( crop(im,2,3,2,1), MPx(255) ); // Blob 2
  fill( crop(im,5,1,1,3), MPx(255) );
  im(6,2) = MPx(5); // Blob 3

  ImageView<uint8> idx =  blob_index( im );
  EXPECT_EQ( 3, int(max_channel_value(idx)) );
  EXPECT_EQ( idx(1,1), idx(2,1) ); // Verify blob 1
  EXPECT_NE( idx(0,0), idx(1,1) );
  EXPECT_EQ( idx(2,4), idx(3,4) ); // Verify blob 2
  EXPECT_NE( idx(1,1), idx(2,4) );
  EXPECT_EQ( idx(5,1), idx(6,2) ); // Verify blob 3
  EXPECT_NE( idx(6,2), idx(1,1) );
}
/*
TEST(BlobIndexThreaded, TestImage1) {
  DiskImageView<PixelGray<uint8> > input("ThreadTest1.tif");
  EXPECT_EQ( 20, input.cols() );
  BlobIndexThreaded bindex( create_mask(input,255), 1000, 10 );
  EXPECT_EQ( 2u, bindex.num_blobs() );
}

TEST(BlobIndexThreaded, TestImage2) {
  DiskImageView<PixelGray<uint8> > input("ThreadTest2.tif");
  EXPECT_EQ( 10, input.cols() );
  BlobIndexThreaded bindex( create_mask(input,255), 1000, 5 );
  EXPECT_EQ( 1u, bindex.num_blobs() );
}

TEST(BlobIndexThreaded, TestImage3) {
  DiskImageView<PixelGray<uint8> > input("ThreadTest3.tif");
  EXPECT_EQ( 10, input.cols() );
  BlobIndexThreaded bindex( create_mask(input,255), 1000, 5 );
  EXPECT_EQ( 2u, bindex.num_blobs() );
}
*/


TEST(BlobIndex, BlobCompressedIntersect) {
  std::vector<std::list<int32> > starts, ends;
  starts += list_of(0), list_of(0), list_of(0), list_of(0), list_of(0);
  ends += list_of(5), list_of(1), list_of(1), list_of(1), list_of(1);
  blob::BlobCompressed test_blob( Vector2i(5,5), starts, ends );

  test_blob.print();

  EXPECT_FALSE( test_blob.intersects( BBox2i(6,6,2,2) ) );
  EXPECT_FALSE( test_blob.intersects( BBox2i(4,4,1,1) ) );
  EXPECT_FALSE( test_blob.intersects( BBox2i(6,4,3,1) ) );
  EXPECT_FALSE( test_blob.intersects( BBox2i(3,5,2,8) ) );
  EXPECT_TRUE( test_blob.intersects( BBox2i(6,2,4,10) ) );
  EXPECT_TRUE( test_blob.intersects( BBox2i(3,4,6,2) ) );
  EXPECT_TRUE( test_blob.intersects( BBox2i(4,7,2,2) ) );
}


TEST( BlobIndex, BlobSizesView ) {

  // Create a test image (zero pixels will be masked)
  ImageView<float> image(6, 6);
  image(0,0) = 0.0;  image(1,0) = 0.0;  image(2,0) = 0.0;  image(3,0) = 1.0;  image(4,0) = 8.0;  image(5,0) = 0.0;
  image(0,1) = 0.0;  image(1,1) = 1.0;  image(2,1) = 0.0;  image(3,1) = 0.0;  image(4,1) = 0.0;  image(5,1) = 0.0;
  image(0,2) = 0.0;  image(1,2) = 1.0;  image(2,2) = 0.0;  image(3,2) = 0.0;  image(4,2) = 1.0;  image(5,2) = 1.0;
  image(0,3) = 1.0;  image(1,3) = 1.0;  image(2,3) = 0.0;  image(3,3) = 0.0;  image(4,3) = 0.0;  image(5,3) = 1.0;
  image(0,4) = 1.0;  image(1,4) = 0.0;  image(2,4) = 1.0;  image(3,4) = 0.0;  image(4,4) = 0.0;  image(5,4) = 0.0;
  image(0,5) = 0.0;  image(1,5) = 1.0;  image(2,5) = 0.0;  image(3,5) = 0.0;  image(4,5) = 0.0;  image(5,5) = 0.0;

  // Run the algorithm
  ImageView<uint32> output = get_blob_sizes(create_mask(image), 256, 256); // Tile sizes don't matter for tiny image

  // Check results
  EXPECT_EQ(0,  output(3,3));
  EXPECT_EQ(7,  output(0,3));
  EXPECT_EQ(2,  output(3,0));
  EXPECT_EQ(3,  output(5,3));
  EXPECT_EQ(0,  output(2,1));
}



