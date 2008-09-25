#ifndef __VW_STEREO_AFFINE_SUBPIXEL_VIEW__
#define __VW_STEREO_AFFINE_SUBPIXEL_VIEW__

#include <vw/Image/ImageView.h>
#include <vw/Image/ImageViewRef.h>
#include <vw/Stereo/DisparityMap.h>
#include <vw/Stereo/Correlate.h>


namespace vw {
namespace stereo {

  /// An image view for performing image correlation
  class AffineSubpixelView : public ImageViewBase<AffineSubpixelView> {

    ImageViewRef<PixelDisparity<float> > m_disparity_map;
    ImageViewRef<float> m_left_image; 
    ImageViewRef<float> m_right_image;
    ImageViewRef<float> m_left_log_image; 
    ImageViewRef<float> m_right_log_image;

    boost::shared_ptr< CropView<ImageView<float> > > m_left_cached_log_image;
    boost::shared_ptr< CropView<ImageView<float> > > m_right_cached_log_image;
    BBox2i m_left_bbox, m_right_bbox;

    // General Settings
    int m_kern_width, m_kern_height;
    bool m_do_h_subpixel, m_do_v_subpixel;
    bool m_verbose;

    // Affine subpixel specific settings
    int m_weight_threshold;
    
    boost::shared_ptr< CropView<ImageView<float> > > m_x_deriv;
    boost::shared_ptr< CropView<ImageView<float> > > m_y_deriv;
    ImageView<float> m_weight_template;

    mutable ImageView<float> m_weight;

    // Private methods
    ImageView<float> compute_gaussian_weight_image(int kern_width, int kern_height) const;      
    void cache(BBox2i bbox);

    template <class DisparityPatchViewT>
    inline int adjust_weight_image(ImageView<float> &weight,
                                   ImageViewBase<DisparityPatchViewT> const& disparity_map_patch,
                                   ImageView<float> const& weight_template) const {
      
      //    const float continuity_threshold_squared = 64;  // T = 8
      int center_pix_x = weight_template.cols()/2;
      int center_pix_y = weight_template.rows()/2;
      PixelDisparity<float> center_pix = disparity_map_patch.impl()(center_pix_x, center_pix_y);
      
      float sum = 0;
      int num_good_pix = 0;
      ImageView<float>::pixel_accessor weight_row_acc = weight.origin();
      ImageView<float>::pixel_accessor template_row_acc = weight_template.origin();
      typename DisparityPatchViewT::pixel_accessor disp_row_acc = disparity_map_patch.impl().origin();
      for (int j = 0; j < weight_template.rows(); ++j) {
        ImageView<float>::pixel_accessor weight_col_acc = weight_row_acc;
        ImageView<float>::pixel_accessor template_col_acc = template_row_acc;
        typename DisparityPatchViewT::pixel_accessor disp_col_acc = disp_row_acc;
        for (int i = 0; i < weight_template.cols(); ++i ) {
          
          // Mask is zero if the disparity map's pixel is missing...
          if ( (*disp_col_acc).missing()) 
            *weight_col_acc = 0;
          
          //         // ... or if there is a large discontinuity ...
          //         if (pow( (*disp_col_acc).h()-center_pix.h(),2) + pow( (*disp_col_acc).v()-center_pix.v(),2) >= continuity_threshold_squared)
          //           *weight_col_acc = 0;
          
          // ... otherwise we use the weight from the weight template
          else {
            *weight_col_acc = *template_col_acc;
            sum += *weight_col_acc;
            ++num_good_pix;
          }
          
          disp_col_acc.next_col();
          weight_col_acc.next_col();
          template_col_acc.next_col();
        }
        disp_row_acc.next_row();
        weight_row_acc.next_row();
        template_row_acc.next_row();
      }
      
      // Normalize the weight image
      if (sum == 0) 
        vw_throw(LogicErr() << "subpixel_weight: Sum of weight image was zero.  This isn't supposed to happen!");
      else 
        weight /= sum;
      return num_good_pix;
    }

  public:
      typedef PixelDisparity<float> pixel_type;
      typedef pixel_type result_type;
      typedef ProceduralPixelAccessor<AffineSubpixelView> pixel_accessor;
        
    template <class DisparityViewT, class InputViewT>
    AffineSubpixelView(DisparityViewT const& disparity_map,
                       InputViewT const& left_image,
                       InputViewT const& right_image,
                       int kern_width, int kern_height,
                       bool do_horizontal_subpixel,
                       bool do_vertical_subpixel,
                       bool verbose) : m_disparity_map( edge_extend( disparity_map,ZeroEdgeExtension() )  ),
                                       m_left_image(left_image),
                                       m_right_image(right_image),
                                       m_kern_width(kern_width), m_kern_height(kern_height),
                                       m_do_h_subpixel(do_horizontal_subpixel),
                                       m_do_v_subpixel(do_vertical_subpixel),
                                       m_verbose(verbose) {
      // Basic assertions
      VW_ASSERT((left_image.impl().cols() == right_image.impl().cols()) &&
                (left_image.impl().rows() == right_image.impl().rows()) &&
                (disparity_map.impl().cols() == right_image.impl().cols()) &&
                (disparity_map.impl().cols() == right_image.impl().cols()),
                ArgumentErr() << "AffineSubpixelView::AffineSubpixelView(): input image dimensions and/or disparity_map dimensions do not agree.\n");
      
      VW_ASSERT((left_image.channels() == 1) && (left_image.impl().planes() == 1) &&
                (right_image.channels() == 1) && (right_image.impl().planes() == 1),
                ArgumentErr() << "AffineSubpixelView::AffineSubpixelView(): multi-channel, multi-plane images not supported.\n");

      m_left_log_image = laplacian_filter(gaussian_filter(channel_cast<float>(left_image.impl()),1.5));
      m_right_log_image = laplacian_filter(gaussian_filter(channel_cast<float>(right_image.impl()),1.5));

      m_left_bbox = BBox2i(0,0,m_left_image.cols(),m_left_image.rows());
      m_right_bbox = BBox2i(0,0,m_right_image.cols(),m_right_image.rows());

      // Robust cost function settings
      int kern_pixels = kern_height * kern_width;
      m_weight_threshold = kern_pixels / 2;

      // Weight template and workspace images are allocated and
      // computed up here out of the tight inner loop.  We rasterize
      // into these directly in the code below.
      m_weight_template = compute_gaussian_weight_image(kern_width, kern_height);
      m_weight.set_size(kern_width, kern_height);

    }

    // Standard ImageView interface methods
    inline int32 cols() const { return m_disparity_map.cols(); }
    inline int32 rows() const { return m_disparity_map.rows(); }
    inline int32 planes() const { return 1; }
    
    inline pixel_accessor origin() const { return pixel_accessor( *this, 0, 0 ); }
    
    inline pixel_type operator()(double x, double y, int32 p = 0) const {

      // Bail out if no subpixel computation has been requested 
      if (!m_do_h_subpixel && !m_do_v_subpixel) return m_disparity_map(x,y);

      if (m_left_cached_log_image) {
        return get_cached_pixel(x,y,p);
      } else {
        return get_uncached_pixel(x,y,p);
      }      
    }

    inline pixel_type get_uncached_pixel(double x, double y, int32 p = 0) const {
      vw_throw(NoImplErr() << "AffineSubpixelView: accessing pixels outside of a rasterization loop is not supported.");
      return PixelDisparity<float>();
    }
    
    pixel_type get_cached_pixel(double x, double y, int32 p = 0) const;

    /// \cond INTERNAL
    typedef AffineSubpixelView prerasterize_type;
    inline prerasterize_type prerasterize(BBox2i bbox) const { 
      prerasterize_type img(m_disparity_map,
                            m_left_image, m_right_image,
                            m_kern_width, m_kern_height,
                            m_do_h_subpixel, m_do_v_subpixel,
                            m_verbose); 
      img.cache(bbox);
      return img;
    }
    
    template <class DestT> inline void rasterize(DestT const& dest, BBox2i bbox) const {
      vw::rasterize(prerasterize(bbox), dest, bbox);
    }
    /// \endcond

  };

}} // namespace vw::stereo

#endif // __VW_STEREO_CORRELATOR_VIEW__         
