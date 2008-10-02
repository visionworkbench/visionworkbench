#include <vw/Stereo/OptimizedCorrelator.h>

using namespace vw;
using namespace stereo;

template <class ScoreT>
struct DisparityScore {
  ScoreT best, worst;
  int32 hdisp, vdisp;
  
  DisparityScore() {
    best = ScalarTypeLimits<ScoreT>::highest();
    worst = ScalarTypeLimits<ScoreT>::lowest();
    hdisp = vdisp = 0;
  }
};

struct AbsDifferenceFunctor : BinaryReturnTemplateType<DifferenceType> {
  template <class Arg1T, class Arg2T>
  typename result<AbsDifferenceFunctor(Arg1T, Arg2T)>::type
  inline operator()(Arg1T const& arg1, Arg2T const& arg2) const { return fabs(arg1-arg2); }
};

template <class Image1T, class Image2T>
inline BinaryPerPixelView<Image1T, Image2T, AbsDifferenceFunctor> abs_difference(ImageViewBase<Image1T> const& image1, ImageViewBase<Image2T> const& image2) {
  return BinaryPerPixelView<Image1T, Image2T, AbsDifferenceFunctor>(image1.impl(), image2.impl(), AbsDifferenceFunctor());
}

struct SqDifferenceFunctor : BinaryReturnTemplateType<DifferenceType> {
  template <class Arg1T, class Arg2T>
  typename result<SqDifferenceFunctor(Arg1T, Arg2T)>::type
  inline operator()(Arg1T const& arg1, Arg2T const& arg2) const {
    typename result<SqDifferenceFunctor(Arg1T, Arg2T)>::type diff = arg1 - arg2;
    return diff * diff;
  }
};

template <class Image1T, class Image2T>
inline BinaryPerPixelView<Image1T, Image2T, SqDifferenceFunctor> sq_difference(ImageViewBase<Image1T> const& image1, ImageViewBase<Image2T> const& image2) {
  return BinaryPerPixelView<Image1T, Image2T, SqDifferenceFunctor>(image1.impl(), image2.impl(), SqDifferenceFunctor());
}

// ---------------------------------------------------------------------------
//                           COST FUNCTIONS
// ---------------------------------------------------------------------------

ImageView<float> AbsDifferenceCost::calculate(int dx, int dy) {
  BBox2i right_bbox = this->bbox() + Vector2i(dx, dy);
  CropView<EdgeExtensionView<ImageView<float>, ZeroEdgeExtension> > left_window(edge_extend(m_left,ZeroEdgeExtension()), this->bbox());
  CropView<EdgeExtensionView<ImageView<float>, ZeroEdgeExtension> > right_window(edge_extend(m_right,ZeroEdgeExtension()), right_bbox);
  return this->box_filter(abs_difference(left_window, right_window));
}


ImageView<float> SqDifferenceCost::calculate(int dx, int dy) {
  BBox2i right_bbox = this->bbox() + Vector2i(dx, dy);
  CropView<EdgeExtensionView<ImageView<float>, ZeroEdgeExtension> > left_window(edge_extend(m_left,ZeroEdgeExtension()), this->bbox());
  CropView<EdgeExtensionView<ImageView<float>, ZeroEdgeExtension> > right_window(edge_extend(m_right,ZeroEdgeExtension()), right_bbox);
  return this->box_filter(abs_difference(left_window, right_window));
}


ImageView<float> NormXCorrCost::calculate(int dx, int dy) {
  BBox2i right_bbox = this->bbox() + Vector2i(dx, dy);
  CropView<EdgeExtensionView<ImageView<float>, ZeroEdgeExtension> > left_window(edge_extend(m_left, ZeroEdgeExtension()), this->bbox());
  CropView<EdgeExtensionView<ImageView<float>, ZeroEdgeExtension> > left_mean_window(edge_extend(m_left_mean, ZeroEdgeExtension()), this->bbox());
  CropView<EdgeExtensionView<ImageView<float>, ZeroEdgeExtension> > left_variance_window(edge_extend(m_left_variance, ZeroEdgeExtension()), this->bbox());
  
  CropView<EdgeExtensionView<ImageView<float>, ZeroEdgeExtension> > right_window(edge_extend(m_right, ZeroEdgeExtension()), right_bbox);
  CropView<EdgeExtensionView<ImageView<float>, ZeroEdgeExtension> > right_mean_window(edge_extend(m_right_mean, ZeroEdgeExtension()), right_bbox);
  CropView<EdgeExtensionView<ImageView<float>, ZeroEdgeExtension> > right_variance_window(edge_extend(m_right_variance, ZeroEdgeExtension()), right_bbox);
  
  ImageView<float> left_right_mean = box_filter(left_window * right_window);
  
  
  // We take the absolute value and subtract 1 here so that
  // invalide (empty) regions of the image (which would normally
  // evaluate to zero) end up evaluating to 1.0, and areas that
  // are highly correlated are closer to zero.
  return abs((pow(left_right_mean - left_mean_window * right_mean_window, 2) /
              left_variance_window / right_variance_window)-1);
}




// ---------------------------------------------------------------------------
//                           CORRELATE()
// ---------------------------------------------------------------------------
ImageView<PixelDisparity<float> > vw::stereo::correlate(boost::shared_ptr<StereoCostFunction> const& cost_function,
                                                        BBox2i const& search_window,
                                                        ProgressCallback const& progress) {
  
  int width = cost_function->cols();
  int height = cost_function->rows();
  
  ImageView<DisparityScore<float> > result_buf(width, height);
  ImageView<float> cost_buf(width, height);
  
  int current_iteration = 0;
  int total_iterations = (search_window.width() + 1) * (search_window.height() + 1);
  
  BBox2i left_bbox = cost_function->bbox();
  for (int dy = search_window.min().y(); dy <= search_window.max().y(); dy++) {
    for (int dx = search_window.min().x(); dx <= search_window.max().x(); dx++) {
      BBox2i right_bbox = cost_function->bbox() + Vector2i (dx,dy);
      //        vw_out(0) << "Correlating: [" << dx << "  " << dy << "]         " << left_bbox << "        " << width << " " << height << "  \n";
      CropView<ImageView<DisparityScore<float> > > result_buf_window(result_buf, left_bbox);
      CropView<ImageView<float> > cost_buf_window(cost_buf, left_bbox);
        
      // Calculate cost function
      cost_buf_window = cost_function->calculate(dx,dy);

      // For debugging
      //
      //         std::ostringstream ostr;
      //         ostr << "costimage_" << dy << "_" << dx << ".tif";
      //         write_image(ostr.str(), cost_buf_window);

      CropView<ImageView<DisparityScore<float> > >::pixel_accessor result_row_acc = result_buf_window.origin();
      CropView<ImageView<float> >::pixel_accessor cost_buf_row_acc = cost_buf_window.origin();
      for (int y = 0; y < left_bbox.height(); y++) {
        CropView<ImageView<DisparityScore<float> > >::pixel_accessor result_col_acc = result_row_acc;
        CropView<ImageView<float> >::pixel_accessor cost_buf_col_acc = cost_buf_row_acc;
        for (int x = 0; x < left_bbox.width(); x++) {
          if (*cost_buf_col_acc < (*result_col_acc).best) {
            (*result_col_acc).best = *cost_buf_col_acc;
            (*result_col_acc).hdisp = dx;
            (*result_col_acc).vdisp = dy;
          }
          if (*cost_buf_col_acc > (*result_col_acc).worst) {
            (*result_col_acc).worst = *cost_buf_col_acc;
          }
          result_col_acc.next_col();
          cost_buf_col_acc.next_col();
        }
        result_row_acc.next_row();
        cost_buf_row_acc.next_row();
      }

      progress.report_fractional_progress(++current_iteration, total_iterations);
      progress.abort_if_requested();
    }
  } 
  //    vw_out(0) << "Correlating: done.                            \n";
    

  // convert from the local result buffer to the return format
  ImageView<PixelDisparity<float> > result(width, height);
    
  for (int x = 0; x < width; x++) {
    for (int y = 0; y < height; y++) {
      if (result_buf(x, y).best == ScalarTypeLimits<float>::highest() || 
          result_buf(x, y).best == result_buf(x, y).worst) {
        result(x, y) = PixelDisparity<float>();
      } else {
        result(x, y) = PixelDisparity<float>(result_buf(x, y).hdisp, result_buf(x, y).vdisp);
      }
    }
  }
  progress.report_finished();
  return result; 
}
  
