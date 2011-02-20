// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifndef __VW_STEREO_OPTIMIZED_CORRELATOR__
#define __VW_STEREO_OPTIMIZED_CORRELATOR__

#include <vw/Image/ImageMath.h>
#include <vw/Stereo/DisparityMap.h>
#include <vw/Stereo/Correlate.h>

// Boost
#include <boost/thread/xtime.hpp>

namespace vw {
namespace stereo {

  // ---------------------------------------------------------------------------
  //                           COST FUNCTIONS
  // ---------------------------------------------------------------------------

  class StereoCostFunction {

  protected:
    BBox2i m_left_bbox;
    ImageView<float> m_dst;
    int32 m_kernel_size;
    float m_kernel_size_i2;
    int32 m_half_kernel;

  public:
    // The constructor allocates the buffers for the box filter once
    // and only once.
    StereoCostFunction(int32 width, int32 height, BBox2i const& search_window, int32 kernel_size) :
    m_left_bbox(BBox2i((search_window.max().x() < 0) ? (-search_window.max().x()) : 0,
                       (search_window.max().y() < 0) ? (-search_window.max().y()) : 0,
                       (search_window.min().x() < 0) ? width - abs(search_window.max().x()) : width - abs(search_window.min().x()),
                       (search_window.min().y() < 0) ? height - abs(search_window.max().y()) : height - abs(search_window.min().y()) )),
      m_dst(m_left_bbox.width(), m_left_bbox.height()),
      m_kernel_size(kernel_size), m_kernel_size_i2(1.0/float(kernel_size*kernel_size)),
      m_half_kernel(kernel_size/2) {}

    virtual ~StereoCostFunction() {}

    BBox2i const& bbox() const { return m_left_bbox; }
    int32 kernel_size() const { return m_kernel_size; }

    virtual ImageView<float> calculate(int32 dx, int32 dy) = 0;
    virtual int32 cols() const = 0;
    virtual int32 rows() const = 0;
    virtual int32 sample_size() const = 0; // What is the side length of
                                         // the square of surrounding
                                         // pixels needed to calculate
                                         // the cost for a single
                                         // pixel?
  protected:

    // Efficient box filter implemenation.  This filter is called
    // repeatedly, but we allocate the image buffers only once (in the
    // constructor, above).
    template <class BoxViewT>
    ImageView<float> box_filter(ImageViewBase<BoxViewT> const& img) {
      VW_ASSERT(img.impl().cols() == m_dst.cols() && img.impl().rows() == m_dst.rows(),
                ArgumentErr() << "StereoCostFunction::box_filter() : image size (" << img.impl().cols() << " " << img.impl().rows() << ") does not match box filter size (" << m_dst.cols() << " " << m_dst.rows() << ").");

      Vector<float> cSum(img.impl().cols());

      // Seed the column sum buffer
      for (int32 x = 0; x < img.impl().cols(); x++) {
        cSum(x) = 0;
        for (int32 ky = 0; ky < m_kernel_size; ky++) {
          cSum(x) += img.impl()(x, ky);
        }
      }

      for (int32 y = 0; y < img.impl().rows() - m_kernel_size; y++) {
        // Seed the row sum
        float rsum = 0;
        for (int32 i = 0; i < m_kernel_size; i++) {
          rsum += cSum(i);
        }

        for (int32 x = 0; x < img.impl().cols() - m_kernel_size; x++) {
          m_dst(x + m_half_kernel, y + m_half_kernel) = rsum;
          // Update the row sum
          rsum += cSum(x + m_kernel_size) - cSum(x);
        }

        // Update the column sum
        for (int32 i = 0; i < img.impl().cols(); i++) {
          cSum(i) += img.impl()(i, y + m_kernel_size) - img.impl()(i, y);
        }
      }

      return m_dst * m_kernel_size_i2;
    }
  };

  class AbsDifferenceCost : public StereoCostFunction {
    ImageView<float> m_left, m_right;

  public:
    template <class ViewT>
    AbsDifferenceCost(ImageViewBase<ViewT> const& left,
                      ImageViewBase<ViewT> const& right,
                      BBox2i const& search_window,
                      int32 kern_size) : StereoCostFunction(left.impl().cols(), left.impl().rows(),
                                                                 search_window, kern_size),
                                       m_left(copy(left.impl())),
                                       m_right(copy(right.impl())) {
      VW_ASSERT(m_left.impl().cols() == m_right.impl().cols(), ArgumentErr() << "Left and right images not the same width");
      VW_ASSERT(m_left.impl().rows() == m_right.impl().rows(), ArgumentErr() << "Left and right images not the same height");
    }

    virtual ImageView<float> calculate(int32 dx, int32 dy);

    virtual int32 cols() const { return m_left.cols(); }
    virtual int32 rows() const { return m_left.rows(); }
    virtual int32 sample_size() const { return this->kernel_size(); }
  };


  class SqDifferenceCost : public StereoCostFunction {
    ImageView<float> m_left, m_right;
  public:
    template <class ViewT>
    SqDifferenceCost(ImageViewBase<ViewT> const& left,
                     ImageViewBase<ViewT> const& right,
                     BBox2i const& search_window,
                     int32 kern_size) : StereoCostFunction(left.impl().cols(),
                                                         left.impl().rows(),
                                                         search_window, kern_size),
                                      m_left(copy(left.impl())),
                                      m_right(copy(right.impl()) ) {
      VW_ASSERT(m_left.cols() == m_right.cols(), ArgumentErr() << "Left and right images not the same width");
      VW_ASSERT(m_left.rows() == m_right.rows(), ArgumentErr() << "Left and right images not the same height");
    }

    virtual ImageView<float> calculate(int32 dx, int32 dy);

    virtual int32 cols() const { return m_left.cols(); }
    virtual int32 rows() const { return m_left.rows(); }
    virtual int32 sample_size() const { return this->kernel_size(); }
  };

  class NormXCorrCost : public StereoCostFunction {
    ImageView<float> m_left, m_left_mean, m_left_variance;
    ImageView<float> m_right, m_right_mean, m_right_variance;

  public:
    template<class ViewT>
    NormXCorrCost(ImageViewBase<ViewT> const& left,
                  ImageViewBase<ViewT> const& right,
                  BBox2i const& search_window,
                  int32 kern_size) : StereoCostFunction(left.impl().cols(), left.impl().rows(),
                                                      search_window, kern_size),
                                   m_left(copy(left.impl())),
                                   m_right(copy(right.impl())) {
      VW_ASSERT(m_left.cols() == m_right.cols(), ArgumentErr() << "Left and right images not the same width");
      VW_ASSERT(m_left.rows() == m_right.rows(), ArgumentErr() << "Left and right images not the same height");

      vw::ImageView<float> left_mean_sq = this->box_filter(m_left * m_left);
      m_left_mean = this->box_filter(m_left);
      m_left_variance = left_mean_sq - m_left_mean * m_left_mean;

      vw::ImageView<float> right_mean_sq = this->box_filter(m_right * m_right);
      m_right_mean = this->box_filter(m_right);
      m_right_variance = right_mean_sq - m_right_mean * m_right_mean;
    }

    virtual ImageView<float> calculate(int32 dx, int32 dy);

    virtual int32 cols() const { return m_left.cols(); }
    virtual int32 rows() const { return m_left.rows(); }
    virtual int32 sample_size() const { return this->kernel_size(); }
  };

  class BlurCost : public StereoCostFunction {
    boost::shared_ptr<StereoCostFunction> m_base_cost;
    int32 m_blur_size;
  public:
    BlurCost(boost::shared_ptr<StereoCostFunction> base_cost,
             BBox2i const& search_window,
             int32 blur_size) : StereoCostFunction(base_cost->cols(), base_cost->rows(),
                                                 search_window, blur_size),
                              m_base_cost(base_cost),
                              m_blur_size(blur_size) {}

    virtual ImageView<float> calculate(int32 dx, int32 dy) {
      return this->box_filter(m_base_cost->calculate(dx, dy));
    }

    int32 cols() const { return m_base_cost->cols(); }
    int32 rows() const { return m_base_cost->rows(); }
    int32 sample_size() const {
      return m_blur_size + m_base_cost->sample_size();
    }
  };

  ImageView<PixelMask<Vector2f> > correlate(boost::shared_ptr<StereoCostFunction> const& cost_function,
                                            BBox2i const& search_window,
                                            ProgressCallback const& progress = ProgressCallback::dummy_instance() );

  class OptimizedCorrelator {

    BBox2i m_search_window;
    int32 m_kern_size;
    float m_cross_correlation_threshold;
    float m_corrscore_rejection_threshold;
    int32 m_cost_blur;
    stereo::CorrelatorType m_correlator_type;

  public:

    // See Correlate.h for CorrelatorType options.
    OptimizedCorrelator(BBox2i const& search_window,
                        int32 const& kernel_size,
                        float const& cross_correlation_threshold,
                        float const& corrscore_rejection_threshold,
                        int32 const& cost_blur = 1,
                        stereo::CorrelatorType correlator_type = ABS_DIFF_CORRELATOR ) :
      m_search_window(search_window),
      m_kern_size(kernel_size),
      m_cross_correlation_threshold(cross_correlation_threshold),
      m_corrscore_rejection_threshold(corrscore_rejection_threshold),
      m_cost_blur(cost_blur),
      m_correlator_type(correlator_type) {}

    template <class ViewT, class PreProcFilterT>
    ImageView<PixelMask<Vector2f> > operator()(ImageViewBase<ViewT> const& image0,
                                               ImageViewBase<ViewT> const& image1,
                                               PreProcFilterT const& preproc_filter) {

      // Check to make sure that image0 and image1 have equal dimensions
      if ((image0.impl().cols() != image1.impl().cols()) ||
          (image0.impl().rows() != image1.impl().rows())) {
        vw_throw( ArgumentErr() << "Primary and secondary image dimensions do not agree!" );
      }

      // Check to make sure that the images are single channel/single plane
      if (!(image0.channels() == 1 && image0.impl().planes() == 1 &&
            image1.channels() == 1 && image1.impl().planes() == 1)) {
        vw_throw( ArgumentErr() << "Both images must be single channel/single plane images!" );
      }

      typedef typename PreProcFilterT::result_type preproc_type;
      preproc_type left_image = preproc_filter(image0);
      preproc_type right_image = preproc_filter(image1);

      BBox2i r2l_window(-m_search_window.max().x(), -m_search_window.max().y(),
                        m_search_window.width(), m_search_window.height());

      boost::shared_ptr<StereoCostFunction> l2r_cost, r2l_cost;

      if (m_correlator_type == ABS_DIFF_CORRELATOR) {
        l2r_cost.reset(new AbsDifferenceCost(left_image, right_image, m_search_window, m_kern_size));
        r2l_cost.reset(new AbsDifferenceCost(right_image, left_image, r2l_window, m_kern_size));
      } else if (m_correlator_type == SQR_DIFF_CORRELATOR) {
        l2r_cost.reset(new SqDifferenceCost(left_image, right_image, m_search_window, m_kern_size));
        r2l_cost.reset(new SqDifferenceCost(right_image, left_image, r2l_window, m_kern_size));
      } else if (m_correlator_type == NORM_XCORR_CORRELATOR) {
        l2r_cost.reset(new NormXCorrCost(left_image, right_image, m_search_window, m_kern_size));
        r2l_cost.reset(new NormXCorrCost(right_image, left_image, r2l_window, m_kern_size));
      } else {
        vw_throw(ArgumentErr() << "OptimizedCorrelator: unknown correlator type " << m_correlator_type << ".");
      }

      boost::shared_ptr<StereoCostFunction> l2r_cost_and_blur = l2r_cost;
      boost::shared_ptr<StereoCostFunction> r2l_cost_and_blur = r2l_cost;
      if (m_cost_blur > 1) {
        l2r_cost_and_blur.reset(new BlurCost(l2r_cost, m_search_window, m_cost_blur));
        r2l_cost_and_blur.reset(new BlurCost(r2l_cost, r2l_window, m_cost_blur));
      }

      ImageView<PixelMask<Vector2f> > result_l2r = stereo::correlate(l2r_cost_and_blur, m_search_window);
      ImageView<PixelMask<Vector2f> > result_r2l = stereo::correlate(r2l_cost_and_blur, r2l_window);

      // Cross check the left and right disparity maps
      cross_corr_consistency_check(result_l2r, result_r2l, m_cross_correlation_threshold, false);

      return result_l2r;
    }
  };



}}   // namespace vw::stereo

#endif // __VW_STEREO_OPTIMIZED_CORRELATOR__
