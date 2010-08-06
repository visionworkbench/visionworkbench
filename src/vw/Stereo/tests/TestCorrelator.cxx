// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


// TestCorrelator.h
#include <gtest/gtest.h>

#include <vw/Image/UtilityViews.h>
#include <vw/Stereo/CorrelatorView.h>
#include <vw/Image/Transform.h>

#include <boost/random/linear_congruential.hpp>

using namespace vw;
using namespace vw::stereo;

//   void test_image_sum_functions() {
//     ImageView<uint8> image = ConstantView<uint8>(1, 100, 100);
//     int kern_width = 10;
//     int kern_height = 10;

//     // Computes the complete sum at every pixel
//     for (int j = kern_height/2; j < image.rows()-kern_height/2; ++j)
//       for (int i = kern_width/2; i < image.cols()-kern_width/2; ++i) {
//         int16 sum = compute_sum(image, BBox2i(i-kern_width/2,j-kern_height/2,kern_width,kern_height));
//         TS_ASSERT_EQUALS(sum,kern_width*kern_height);
//       }

//     // Computes the complete sum at the beginning of the scanline, and
//     // updates as it moves across the image.
//     for (int j = kern_height/2; j < image.rows()-kern_height/2; ++j) {
//       bool need_to_compute_full_sum = true;
//       for (int i = kern_width/2; i < image.cols()-kern_width/2; ++i) {
//         int16 sum;
//         if (need_to_compute_full_sum) {
//           sum = compute_sum(image, BBox2i(i-kern_width/2,j-kern_height/2,kern_width,kern_height));
//           need_to_compute_full_sum = false;
//         } else {
//           sum = update_sum(image, sum, BBox2i(i-kern_width/2,j-kern_height/2,kern_width,kern_height));
//         }
//         TS_ASSERT_EQUALS(sum,kern_width*kern_height);
//       }
//     }
//   }

class BasicCorrelationTest : public ::testing::Test {
protected:
  BasicCorrelationTest() {}

  virtual void SetUp() {
    boost::rand48 gen(10);
    image1 = 255*uniform_noise_view( gen, 50, 50 );
    image2 = transform(image1, TranslateTransform(3,3),
                       ZeroEdgeExtension(), NearestPixelInterpolation());
    mask.set_size(50,50);
    fill(mask,PixelMask<uint8>(255));
  }

  template <class ViewT, class MViewT, class PreProcT>
  CorrelatorView<typename ViewT::pixel_type,
                 typename MViewT::pixel_type,
                 PreProcT>
  correlate( ImageViewBase<ViewT> const& input1,
             ImageViewBase<ViewT> const& input2,
             ImageViewBase<MViewT> const& mask,
             PreProcT const& proc) {
    typedef typename ViewT::pixel_type pixel_type;
    typedef typename MViewT::pixel_type mask_type;
    CorrelatorView<pixel_type,mask_type,PreProcT> corr( input1, input2, mask, mask, proc, false );
    corr.set_search_range( BBox2i(0,0,6,6) );
    corr.set_kernel_size(  Vector2i(7,7)   );
    return corr;
  }

  template <class ViewT>
  void check_error( ImageViewBase<ViewT> const& input,
                    float success = 0.9 ) {
    ViewT const& disparity_map = input.impl();
    int count_correct = 0;
    int count_valid = 0;
    for (int j = 0; j < disparity_map.rows(); ++j)
      for (int i = 0; i < disparity_map.cols(); ++i)
        if ( is_valid( disparity_map(i,j) ) ) {
          count_valid++;
          if ( disparity_map(i,j).child() == Vector2f(3,3) )
            count_correct++;
        }
    std::cout << "Correct Amount: "
              << float(count_correct)/float(count_valid) << "\n";
    EXPECT_LT( success, float(count_correct)/float(count_valid) );
  }

  ImageView<uint8> image1, image2;
  ImageView<PixelMask<uint8> > mask;
};

TEST_F( BasicCorrelationTest, NullPreprocess ) {
  ImageView<PixelMask<Vector2f> > disparity_map =
    correlate( image1, image2, mask,
               NullStereoPreprocessingFilter() );
  check_error( disparity_map, 0.95 );
}

TEST_F( BasicCorrelationTest, SlogPreprocess ) {
  ImageView<PixelMask<Vector2f> > disparity_map =
    correlate( image1, image2, mask,
               SlogStereoPreprocessingFilter() );
  check_error( disparity_map, 0.80 );
}

TEST_F( BasicCorrelationTest, LogPreprocess ) {
  ImageView<PixelMask<Vector2f> > disparity_map =
    correlate( image1, image2, mask,
               LogStereoPreprocessingFilter() );
  check_error( disparity_map, 0.80 );
}

TEST_F( BasicCorrelationTest, BlurPreprocess ) {
  ImageView<PixelMask<Vector2f> > disparity_map =
    correlate( image1, image2, mask,
               BlurStereoPreprocessingFilter() );
  check_error( disparity_map, 0.75 );
}
