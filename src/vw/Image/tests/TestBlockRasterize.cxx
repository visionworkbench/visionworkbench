// __BEGIN_LICENSE__
//  Copyright (c) 2006-2013, United States Government as represented by the
//  Administrator of the National Aeronautics and Space Administration. All
//  rights reserved.
//
//  The NASA Vision Workbench is licensed under the Apache License,
//  Version 2.0 (the "License"); you may not use this file except in
//  compliance with the License. You may obtain a copy of the License at
//  http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
// __END_LICENSE__


#include <test/Helpers.h>
#include <vw/Image/BlockRasterize.h>
#include <vw/Image/BlockImageOperator.h>

using namespace vw;
using namespace std;

TEST(BlockRasterize, Basic) {
  typedef ImageView<uint32> Image;
  typedef BlockRasterizeView<Image> Block;
  Vector2i block(1,1);

  Image img1(2,2), img2;
  img1(0,0) = 1; img1(0,1) = 2; img1(1,0) = 3; img1(1,1) = 4;

  Block b1(img1, block, 4);
  Block b2 = block_rasterize(img1, block, 4);
  Block b3 = block_cache(img1, block, 4);
  Block b4 = block_cache(img1, block, 4, vw_system_cache());

  img2 = b1;
  EXPECT_RANGE_EQ(img1.begin(), img1.end(), img2.begin(), img2.end());
  img2 = b2;
  EXPECT_RANGE_EQ(img1.begin(), img1.end(), img2.begin(), img2.end());
  img2 = b3;
  EXPECT_RANGE_EQ(img1.begin(), img1.end(), img2.begin(), img2.end());
  img2 = b4;
  EXPECT_RANGE_EQ(img1.begin(), img1.end(), img2.begin(), img2.end());
}

/// Count the number of pixels above a threshold on a per-block basis.
class ImageBlockThresholdFunctor {
  
  double m_threshold;
  size_t m_count;
  Mutex  m_mutex;
  
public:
  
  ImageBlockThresholdFunctor(double threshold)
    : m_threshold(threshold), m_count(0) {}

  size_t get_count() const {return m_count;}

  void reset_count() { m_count = 0; }

  // This is the only thread safe function!
  template <class T>
  void operator()(ImageView<T> const& image, BBox2i const& bbox) {
    // Count values in this region.
    size_t local_count = 0;
    for (int r=0; r<bbox.height(); ++r) {
      for (int c=0; c<bbox.width(); ++c) {
        if (image(c,r) >= m_threshold) {
          ++local_count;
        }
      }
    }
    // Safely add to the shared total.
    m_mutex.lock();
    m_count += local_count;
    m_mutex.unlock();
  }
};


TEST(BlockImageOperator, Basic) {

  const double   threshold = 17.0;
  const Vector2i block_size(15, 15);

  ImageBlockThresholdFunctor threshold_functor(threshold);

  // Generate a test image and compute the correct answer.
  size_t real_count = 0;
  ImageView<uint8> image(100,100);
  for (int i=0; i<100; ++i) {
      for (int j=0; j<100; ++j) {
      uint8 value = i/7;
      image(i,j) = value;
      if (value >= threshold)
        ++real_count;
    }
  }

  // Use BlockImageOperator to count thresholds in multiple threads.
  block_op_cache(image, threshold_functor, block_size);
  size_t result = threshold_functor.get_count();
  EXPECT_EQ(real_count, result);

  threshold_functor.reset_count();
  EXPECT_EQ(0, threshold_functor.get_count());

  // Try it again without using the cache.
  block_op(image, threshold_functor, block_size);
  result = threshold_functor.get_count();
  EXPECT_EQ(real_count, result);
}
