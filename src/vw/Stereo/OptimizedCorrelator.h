#ifndef __VW_STEREO_OPTIMIZED_CORRELATOR__
#define __VW_STEREO_OPTIMIZED_CORRELATOR__

#include <vw/Core/Stopwatch.h>
#include <vw/Image/ImageView.h>
#include <vw/Image/Manipulation.h>
#include <vw/Stereo/DisparityMap.h>
#include <vw/Stereo/Correlate.h>

// Boost
#include <boost/thread/xtime.hpp>

namespace vw { 
namespace stereo {
  // ---------------------------------------------------------------------------
  //                        STANDALONE BOX FILTER
  // ---------------------------------------------------------------------------

  // Efficient box filter implemenation.  
  template <class ViewT>
  static ImageView<float> standalone_box_filter(ImageViewBase<ViewT> const& img, Vector2i kernel_size) {
    int kern_width = kernel_size.x();
    int kern_height = kernel_size.y();

    ImageView<float> src = img.impl();
    ImageView<float> dst(src.cols(), src.rows());
    std::vector<float> rSum(src.rows());

    typedef typename ImageView<float>::pixel_accessor src_accessor;
    typedef typename ImageView<float>::pixel_accessor dst_accessor;
    
    // Seed the row sum buffer
    src_accessor row_acc = src.origin();
    for (int y = 0; y < src.rows(); ++y) {
      rSum[y] = 0;
      src_accessor col_acc = row_acc;
      for (int kx = 0; kx < kern_width; ++kx) {
        rSum[y] += *col_acc; // src(kx, y);
        col_acc.next_col();
      }
      row_acc.next_row();
    }
      
    dst_accessor dst_col_acc = dst.origin();
    dst_col_acc.advance(kern_width / 2, kern_height / 2);
    for (int x = 0; x < src.cols() - kern_width; ++x) {
      // Seed the col sum
      float cSum = 0;
      for (int j = 0; j < kern_height; ++j) {
        cSum += rSum[j];
      }
      
      dst_accessor dst_row_acc = dst_col_acc;
      for (int y = 0; y < src.rows() - kern_height; ++y) {
        *dst_row_acc = cSum;  // dst(x + kern_width / 2, y + kern_height / 2) = cSum
        dst_row_acc.next_row();
        
        // Update the col sum
        cSum += rSum[y + kern_height] - rSum[y];
      }
      
      // Update the row sum
      for (int j = 0; j < src.rows(); ++j) {
        rSum[j] += src(x + kern_width, j) - src(x, j);
      }
      dst_col_acc.next_col();
    }
    
    return dst;
  }

  // ---------------------------------------------------------------------------
  //                           COST FUNCTIONS
  // ---------------------------------------------------------------------------

  class StereoCostFunction {

  protected:
    BBox2i m_left_bbox;
    ImageView<float> m_src;
    ImageView<float> m_dst;
    std::vector<float> m_rSum;
    int m_kernel_size;

  public:
    // The constructor allocates the buffers for the box filter once
    // and only once.
    StereoCostFunction(int width, int height, BBox2i const& search_window, int kernel_size) :     
      m_left_bbox(BBox2i((search_window.max().x() < 0) ? (-search_window.max().x()) : 0,
                         (search_window.max().y() < 0) ? (-search_window.max().y()) : 0,
                         (search_window.min().x() < 0) ? width - abs(search_window.max().x()) : width - abs(search_window.min().x()),
                         (search_window.min().y() < 0) ? height - abs(search_window.max().y()) : height - abs(search_window.min().y()) )),
      m_src(m_left_bbox.width(), m_left_bbox.height()), 
      m_dst(m_left_bbox.width(), m_left_bbox.height()), 
      m_rSum(m_left_bbox.height()), 
      m_kernel_size(kernel_size) {}

    BBox2i const& bbox() const { return m_left_bbox; }
    int kernel_size() const { return m_kernel_size; }

    virtual ImageView<float> calculate(int dx, int dy) = 0;
    virtual int cols() const = 0;
    virtual int rows() const = 0;
    virtual int sample_size() const = 0; // What is the side length of
                                         // the square of surrounding
                                         // pixels needed to calculate
                                         // the cost for a single
                                         // pixel?
    virtual ~StereoCostFunction() {}

    // Efficient box filter implemenation.  This filter is called
    // repeatedly, but we allocate the image buffers only once (in the
    // constructor, above).
    template <class BoxViewT>
    ImageView<float> box_filter(ImageViewBase<BoxViewT> const& img) {
      int kern_width = m_kernel_size;
      int kern_height = m_kernel_size;

      VW_ASSERT(img.impl().cols() == m_dst.cols() && img.impl().rows() == m_dst.rows(),
                ArgumentErr() << "StereoCostFunction::box_filter() : image size (" << img.impl().cols() << " " << img.impl().rows() << ") does not match box filter size (" << m_dst.cols() << " " << m_dst.rows() << ").");
      m_src = img.impl();

      typedef typename ImageView<float>::pixel_accessor src_accessor;
      typedef typename ImageView<float>::pixel_accessor dst_accessor;
      
      // Seed the row sum buffer
      src_accessor row_acc = m_src.origin();
      for (int y = 0; y < m_src.rows(); ++y) {
        m_rSum[y] = 0;
        src_accessor col_acc = row_acc;
        for (int kx = 0; kx < kern_width; ++kx) {
          m_rSum[y] += *col_acc; // src(kx, y);
          col_acc.next_col();
        }
        row_acc.next_row();
      }
      
      dst_accessor dst_col_acc = m_dst.origin();
      dst_col_acc.advance(kern_width / 2, kern_height / 2);
      for (int x = 0; x < m_src.cols() - kern_width; ++x) {
        // Seed the col sum
        float cSum = 0;
        for (int j = 0; j < kern_height; ++j) {
          cSum += m_rSum[j];
        }
        
        dst_accessor dst_row_acc = dst_col_acc;
        for (int y = 0; y < m_src.rows() - kern_height; ++y) {
          *dst_row_acc = cSum;  // dst(x + kern_width / 2, y + kern_height / 2) = cSum
          dst_row_acc.next_row();
          
          // Update the col sum
          cSum += m_rSum[y + kern_height] - m_rSum[y];
        }
        
        // Update the row sum
        for (int j = 0; j < m_src.rows(); ++j) {
          m_rSum[j] += m_src(x + kern_width, j) - m_src(x, j);
        }
        dst_col_acc.next_col();
      }
      
      return m_dst;
    }
  };

  class AbsDifferenceCost : public StereoCostFunction {
    ImageView<float> m_left, m_right;

  public:
    template <class ViewT>
    AbsDifferenceCost(ImageViewBase<ViewT> const& left, 
                      ImageViewBase<ViewT> const& right, 
                      BBox2i const& search_window,
                      int kern_size) : StereoCostFunction(left.impl().cols(), left.impl().rows(), 
                                                                 search_window, kern_size), 
                                       m_left(copy(left.impl())),
                                       m_right(copy(right.impl())) {
      VW_ASSERT(m_left.impl().cols() == m_right.impl().cols(), ArgumentErr() << "Left and right images not the same width");
      VW_ASSERT(m_left.impl().rows() == m_right.impl().rows(), ArgumentErr() << "Left and right images not the same height");
    }
    
    virtual ImageView<float> calculate(int dx, int dy);
    
    virtual int cols() const { return m_left.cols(); }
    virtual int rows() const { return m_left.rows(); }
    virtual int sample_size() const { return this->kernel_size(); }
  };

  
  class SqDifferenceCost : public StereoCostFunction {
    ImageView<float> m_left, m_right;
  public:
    template <class ViewT>
    SqDifferenceCost(ImageViewBase<ViewT> const& left, 
                     ImageViewBase<ViewT> const& right, 
                     BBox2i const& search_window,
                     int kern_size) : StereoCostFunction(left.impl().cols(), 
                                                         left.impl().rows(), 
                                                         search_window, kern_size), 
                                      m_left(copy(left.impl())),
                                      m_right(copy(right.impl()) ) {
      VW_ASSERT(m_left.cols() == m_right.cols(), ArgumentErr() << "Left and right images not the same width");
      VW_ASSERT(m_left.rows() == m_right.rows(), ArgumentErr() << "Left and right images not the same height");
    }
    
    virtual ImageView<float> calculate(int dx, int dy);
    
    virtual int cols() const { return m_left.cols(); }
    virtual int rows() const { return m_left.rows(); }
    virtual int sample_size() const { return this->kernel_size(); }
  };
  
  class NormXCorrCost : public StereoCostFunction {
    ImageView<float> m_left, m_left_mean, m_left_variance;
    ImageView<float> m_right, m_right_mean, m_right_variance;
    
  public:
    template<class ViewT>
    NormXCorrCost(ImageViewBase<ViewT> const& left, 
                  ImageViewBase<ViewT> const& right, 
                  BBox2i const& search_window,
                  int kern_size) : StereoCostFunction(left.impl().cols(), left.impl().rows(), 
                                                      search_window, kern_size), 
                                   m_left(copy(left.impl())), 
                                   m_right(copy(right.impl())) {
      VW_ASSERT(m_left.cols() == m_right.cols(), ArgumentErr() << "Left and right images not the same width");
      VW_ASSERT(m_left.rows() == m_right.rows(), ArgumentErr() << "Left and right images not the same height");
      
      vw::ImageView<float> left_mean_sq = standalone_box_filter(m_left * m_left, Vector2i(kern_size,kern_size));
      m_left_mean = standalone_box_filter(m_left, Vector2i(kern_size,kern_size));
      m_left_variance = left_mean_sq - m_left_mean * m_left_mean;
      
      vw::ImageView<float> right_mean_sq = standalone_box_filter(m_right * m_right, Vector2i(kern_size,kern_size));
      m_right_mean = standalone_box_filter(m_right, Vector2i(kern_size,kern_size));
      m_right_variance = right_mean_sq - m_right_mean * m_right_mean;
    }
    
    virtual ImageView<float> calculate(int dx, int dy);
    
    virtual int cols() const { return m_left.cols(); }
    virtual int rows() const { return m_left.rows(); }
    virtual int sample_size() const { return this->kernel_size(); }
  };
  
  class BlurCost : public StereoCostFunction {
    boost::shared_ptr<StereoCostFunction> m_base_cost;
    int m_blur_size;
  public:
    BlurCost(boost::shared_ptr<StereoCostFunction> base_cost,
             BBox2i const& search_window,
             int blur_size) : StereoCostFunction(base_cost->cols(), base_cost->rows(), 
                                                 search_window, blur_size), 
                              m_base_cost(base_cost), 
                              m_blur_size(blur_size) {}
    
    virtual ImageView<float> calculate(int dx, int dy) {
      return this->box_filter(m_base_cost->calculate(dx, dy));
    }
    
    int cols() const { return m_base_cost->cols(); }
    int rows() const { return m_base_cost->rows(); }
    int sample_size() const { 
      return m_blur_size + m_base_cost->sample_size();
    }
  };

  ImageView<PixelDisparity<float> > correlate(boost::shared_ptr<StereoCostFunction> const& cost_function,
                                              BBox2i const& search_window,
                                              ProgressCallback const& progress = ProgressCallback::dummy_instance() );  



  class OptimizedCorrelator {
    
    BBox2i m_search_window;
    int m_kern_size;
    double m_cross_correlation_threshold;
    double m_corrscore_rejection_threshold;
    int m_cost_blur;
    stereo::CorrelatorType m_correlator_type;

  public:
    
    // See Correlate.h for CorrelatorType options.
    OptimizedCorrelator(BBox2i const& search_window,
                        int const kernel_size,
                        int cross_correlation_threshold,
                        float corrscore_rejection_threshold,
                        int cost_blur = 1,
                        stereo::CorrelatorType correlator_type = ABS_DIFF_CORRELATOR ) : 
      m_search_window(search_window),
      m_kern_size(kernel_size),
      m_cross_correlation_threshold(cross_correlation_threshold),
      m_corrscore_rejection_threshold(corrscore_rejection_threshold),
      m_cost_blur(cost_blur),
      m_correlator_type(correlator_type) {}
      
    template <class ViewT, class PreProcFilterT>
    ImageView<PixelDisparity<float> > operator()(ImageViewBase<ViewT> const& image0,
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

      //Run the correlator and record how long it takes to run.
      Stopwatch timer;
      timer.start();
      
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
      
      ImageView<PixelDisparity<float> > result_l2r = stereo::correlate(l2r_cost_and_blur, m_search_window);
      ImageView<PixelDisparity<float> > result_r2l = stereo::correlate(r2l_cost_and_blur, r2l_window);

      timer.stop();
      //      double lapse__ = timer.elapsed_seconds();
      
      // Cross check the left and right disparity maps
      cross_corr_consistency_check(result_l2r, result_r2l, m_cross_correlation_threshold, false);

      int matched = 0;
      int total = 0;
      int nn = 0;
      for (int j = 0; j < result_l2r.rows(); j++) {
        for (int i = 0; i < result_l2r.cols(); i++) {
          total++;
          if (!(result_l2r(i,j).missing())) {
            matched++;
          }
          nn++;
        }
      } 

      //      vw_out(InfoMessage, "stereo")
      //      vw_out(0) << "\tCorrelation took " << lapse__ << " sec";
      //      double nTries = (m_search_window.max().y() - m_search_window.min().y() + 1) * (m_search_window.max().x() - m_search_window.min().x() + 1);
      //      double rate = nTries * left_image.cols() * left_image.rows() / lapse__ / 1.0e6;
      //      vw_out(0) << "\t(" << rate << " M disparities/second)\n";
      
      //      double score = (100.0 * matched) / total;
      //      vw_out(0) << "\tCorrelation rate: " << score << "\n\n";
    
      return result_l2r;
    }      
  };



}}   // namespace vw::stereo

#endif // __VW_STEREO_OPTIMIZED_CORRELATOR__
