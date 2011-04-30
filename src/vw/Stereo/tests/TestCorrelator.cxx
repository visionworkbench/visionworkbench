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
             PreProcT const& proc,
             stereo::CorrelatorType type = stereo::ABS_DIFF_CORRELATOR ) {
    typedef typename ViewT::pixel_type pixel_type;
    typedef typename MViewT::pixel_type mask_type;
    CorrelatorView<pixel_type,mask_type,PreProcT> corr( input1, input2, mask, mask, proc, false );
    corr.set_search_range( BBox2i(0,0,6,6) );
    corr.set_kernel_size(  Vector2i(7,7)   );
    corr.set_correlator_options( 1, type );
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
    EXPECT_GT( float(count_correct)/float(count_valid), success );
  }

  ImageView<uint8> image1, image2;
  ImageView<PixelMask<uint8> > mask;
};

TEST_F( BasicCorrelationTest, NullPreprocess ) {
  typedef NullStereoPreprocessingFilter FilterT;

  ImageView<PixelMask<Vector2f> > disparity_map =
    correlate( image1, image2, mask, FilterT(),
               stereo::ABS_DIFF_CORRELATOR );
  check_error( disparity_map, 0.95 );

  disparity_map =
    correlate( image1, image2, mask, FilterT(),
               stereo::SQR_DIFF_CORRELATOR );
  check_error( disparity_map, 0.93 );

  disparity_map =
    correlate( image1, image2, mask, FilterT(),
               stereo::NORM_XCORR_CORRELATOR );
  check_error( disparity_map, 0.95 );
}

TEST_F( BasicCorrelationTest, SlogPreprocess ) {
  typedef SlogStereoPreprocessingFilter FilterT;

  ImageView<PixelMask<Vector2f> > disparity_map =
    correlate( image1, image2, mask, FilterT(),
               stereo::ABS_DIFF_CORRELATOR );
  check_error( disparity_map, 0.84 );

  disparity_map =
    correlate( image1, image2, mask, FilterT(),
               stereo::SQR_DIFF_CORRELATOR );
  check_error( disparity_map, 0.84 );

  disparity_map =
    correlate( image1, image2, mask, FilterT(),
               stereo::NORM_XCORR_CORRELATOR );
  check_error( disparity_map, 0.915 );
}

TEST_F( BasicCorrelationTest, LogPreprocess ) {
  typedef LogStereoPreprocessingFilter FilterT;

  ImageView<PixelMask<Vector2f> > disparity_map =
    correlate( image1, image2, mask, FilterT(),
               stereo::ABS_DIFF_CORRELATOR );
  check_error( disparity_map, 0.85 );

  disparity_map =
    correlate( image1, image2, mask, FilterT(),
               stereo::SQR_DIFF_CORRELATOR );
  check_error( disparity_map, 0.81 );

  disparity_map =
    correlate( image1, image2, mask, FilterT(),
               stereo::NORM_XCORR_CORRELATOR );
  check_error( disparity_map, 0.89 );
}

TEST_F( BasicCorrelationTest, BlurPreprocess ) {
  typedef BlurStereoPreprocessingFilter FilterT;

  ImageView<PixelMask<Vector2f> > disparity_map =
    correlate( image1, image2, mask, FilterT(),
               stereo::ABS_DIFF_CORRELATOR );
  check_error( disparity_map, 0.79 );

  disparity_map =
    correlate( image1, image2, mask, FilterT(),
               stereo::SQR_DIFF_CORRELATOR );
  check_error( disparity_map, 0.75 );

  disparity_map =
    correlate( image1, image2, mask, FilterT(),
               stereo::NORM_XCORR_CORRELATOR );
  check_error( disparity_map, 0.79 );
}
