#ifndef __VW_STEREO_AFFINE_SUBPIXEL_VIEW__
#define __VW_STEREO_AFFINE_SUBPIXEL_VIEW__

#include <vw/Image/ImageView.h>
#include <vw/Stereo/DisparityMap.h>
#include <vw/Stereo/Correlate.h>


namespace vw {
namespace stereo {

  /// An image view for performing image correlation
  template<class InputViewT, class DisparityViewT>
  class AffineSubpixelView : public ImageViewBase<AffineSubpixelView<InputViewT, DisparityViewT> > {

    struct HuberError { 
      double m_b;
      HuberError(double b = 0.01) : m_b(b) {}
      
      double operator() (double delta_norm) const {
        if (delta_norm < m_b)
          return delta_norm*delta_norm;
        else
          return 2*m_b*delta_norm - m_b*m_b;
      }
    };

    DisparityViewT m_disparity_map;
    InputViewT m_left_image; 
    InputViewT m_right_image;
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
    HuberError m_robust_cost_fn;
    int m_weight_threshold;
    
    boost::shared_ptr< CropView<ImageView<float> > > m_x_deriv;
    boost::shared_ptr< CropView<ImageView<float> > > m_y_deriv;
    ImageView<float> m_weight_template;

    mutable ImageView<float> m_weight;

    inline ImageView<float> compute_gaussian_weight_image(int kern_width, int kern_height) const {
      
      int center_pix_x = kern_width/2;
      int center_pix_y = kern_height/2;
      int two_sigma_sqr = 2*pow(kern_width/3,2);
      
      ImageView<float> weight(kern_width, kern_height);
      for (int j = 0; j < kern_height; ++j) {
        for (int i = 0; i < kern_width; ++i ) {
          weight(i,j) = exp(-1 * (pow(i-center_pix_x,2) + pow(j-center_pix_y,2)) / two_sigma_sqr);
        }
      }
      return weight;
    }

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
        
    AffineSubpixelView(DisparityViewT const& disparity_map,
                           InputViewT const& left_image,
                           InputViewT const& right_image,
                           int kern_width, int kern_height,
                           bool do_horizontal_subpixel,
                           bool do_vertical_subpixel,
                           bool verbose) : m_disparity_map(disparity_map),
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

      m_left_log_image = laplacian_filter(gaussian_filter(channel_cast<float>(channels_to_planes(left_image.impl())),1.5));
      m_right_log_image = laplacian_filter(gaussian_filter(channel_cast<float>(channels_to_planes(right_image.impl())),1.5));

      m_left_bbox = BBox2i(0,0,m_left_image.cols(),m_left_image.rows());
      m_right_bbox = BBox2i(0,0,m_right_image.cols(),m_right_image.rows());

      // Robust cost function settings
      float thresh = 0.01;
      m_robust_cost_fn = HuberError(thresh);
      int kern_pixels = kern_height * kern_width;
      m_weight_threshold = kern_pixels / 4;

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
    
    inline pixel_type get_cached_pixel(double x, double y, int32 p = 0) const {

      // Skip pixels along the border
      if (x < m_left_bbox.min().x() || x >= m_left_bbox.max().x() || y < m_left_bbox.min().y() || y >= m_left_bbox.max().y()) 
        return PixelDisparity<float>();

      // Skip over pixels for which we have no initial disparity estimate
      if (m_disparity_map(x,y).missing())
        return PixelDisparity<float>();

      BBox2i current_window(x-m_kern_width/2, y-m_kern_height/2, m_kern_width, m_kern_height);
      Vector2 base_offset( -m_disparity_map(x,y).h() , -m_disparity_map(x,y).v() );          
        
      // Initialize our affine transform with the identity.  The
      // entries of d are laid out in row major order:
      // 
      //   | d(0) d(1) d(2) | 
      //   | d(3) d(4) d(5) |
      //   |  0    0    1   |
      //
      Vector<float,6> d;
      d(0) = 1.0;
      d(4) = 1.0;
      Vector2 offset;

      // Compute the derivative image patches
      CropView<CropView<ImageView<float> > > left_image_patch = crop(*m_left_cached_log_image, current_window);
      CropView<CropView<ImageView<float> > > I_x = crop(*m_x_deriv, current_window);
      CropView<CropView<ImageView<float> > > I_y = crop(*m_y_deriv, current_window);
        
      // Compute the base weight image
      int good_pixels = adjust_weight_image(m_weight, crop(m_disparity_map, current_window), m_weight_template);
        
      // Skip over pixels for which there are very few good matches
      // in the neighborhood.
      if (good_pixels < m_weight_threshold) 
        return PixelDisparity<float>();
                
      // Iterate until a solution is found or the max number of
      // iterations is reached.
      for (unsigned iter = 0; iter < 10; ++iter) {
        offset(0) = d[2];
        offset(1) = d[5];
        
        // First we check to see if our current subpixel translation
        // is less than one half of the window width.  If not, then
        // we are probably having trouble converging and we abort
        // this pixel!!
        if (norm_2(offset) > m_kern_width/2) 
          break;
          
        InterpolationView<EdgeExtensionView<CropView<ImageView<float> >, ZeroEdgeExtension>, BilinearInterpolation> right_interp_image =
          interpolate(*m_right_cached_log_image, BilinearInterpolation(), ZeroEdgeExtension());
        
        float x_base = x + m_disparity_map(x,y).h();
        float y_base = y + m_disparity_map(x,y).v();
        //          float error_total = 0;

        Matrix<float,6,6> rhs;
        Vector<float,6> lhs;
        for (int jj = -m_kern_height/2; jj <= m_kern_height/2; ++jj) {
          for (int ii = -m_kern_width/2; ii <= m_kern_width/2; ++ii) {
            int i = ii + m_kern_width/2;
            int j = jj + m_kern_height/2;

            // First we compute the pixel offset for the right image
            // and the error for the current pixel.
            float xx = x_base + d[0] * ii + d[1] * jj + offset(0);
            float yy = y_base + d[3] * ii + d[4] * jj + offset(1);
            float I_e_val = right_interp_image(xx,yy) - left_image_patch(i,j) + 1e-16; 

            // Apply the robust cost function.  We use a huber
            // function to gently remove outliers for small errors,
            // but we set a hard limit a 5 times the cost threshold
            // to remove major (salt&pepper) noise.
            //              if (fabs(I_e_val) > thresh*5 || I_e_val == 0.0)
            float robust_weight = sqrt(m_robust_cost_fn(fabs(I_e_val)))/fabs(I_e_val);
            //              error_total += pow(I_e_val,2);

            // We combine the error value with the derivative and
            // add this to the update equation.
            float I_x_val = robust_weight * m_weight(i,j) * I_x(i,j);
            float I_y_val = robust_weight * m_weight(i,j) * I_y(i,j);
            float I_x_sqr = robust_weight * m_weight(i,j) * I_x(i,j) * I_x(i,j);
            float I_y_sqr = robust_weight * m_weight(i,j) * I_y(i,j) * I_y(i,j);
            float I_x_I_y = robust_weight * m_weight(i,j) * I_x(i,j) * I_y(i,j);

            // Left hand side
            lhs(0) += ii * I_x_val * I_e_val;
            lhs(1) += jj * I_x_val * I_e_val;
            lhs(2) +=      I_x_val * I_e_val;
            lhs(3) += ii * I_y_val * I_e_val;
            lhs(4) += jj * I_y_val * I_e_val;
            lhs(5) +=      I_y_val * I_e_val;
              
            // Right Hand Side UL
            rhs(0,0) += ii*ii * I_x_sqr;
            rhs(0,1) += ii*jj * I_x_sqr;
            rhs(0,2) += ii    * I_x_sqr;
            rhs(1,1) += jj*jj * I_x_sqr;
            rhs(1,2) += jj    * I_x_sqr;
            rhs(2,2) +=         I_x_sqr;
            
            // Right Hand Side UR
            rhs(0,3) += ii*ii * I_x_I_y;
            rhs(0,4) += ii*jj * I_x_I_y;
            rhs(0,5) += ii    * I_x_I_y;
            rhs(1,4) += jj*jj * I_x_I_y;
            rhs(1,5) += jj    * I_x_I_y;
            rhs(2,5) +=         I_x_I_y;
            
            // Right Hand Side LR
            rhs(3,3) += ii*ii * I_y_sqr;
            rhs(3,4) += ii*jj * I_y_sqr;
            rhs(3,5) += ii    * I_y_sqr;
            rhs(4,4) += jj*jj * I_y_sqr;
            rhs(4,5) += jj    * I_y_sqr;
              rhs(5,5) +=         I_y_sqr;
          }
        }

        lhs *= -1;

        // Fill in symmetric entries
        rhs(1,0) = rhs(0,1);
        rhs(2,0) = rhs(0,2);
        rhs(2,1) = rhs(1,2);
        rhs(1,3) = rhs(0,4);
        rhs(2,3) = rhs(0,5);
        rhs(2,4) = rhs(1,5);
        rhs(3,0) = rhs(0,3);
        rhs(3,1) = rhs(1,3);
        rhs(3,2) = rhs(2,3);
        rhs(4,0) = rhs(0,4);
        rhs(4,1) = rhs(1,4);
        rhs(4,2) = rhs(2,4);
        rhs(4,3) = rhs(3,4);
        rhs(5,0) = rhs(0,5);
        rhs(5,1) = rhs(1,5);
        rhs(5,2) = rhs(2,5);
        rhs(5,3) = rhs(3,5);
        rhs(5,4) = rhs(4,5);

        //           {          
        //             ImageView<float> right_image_patch(kern_width, kern_height);
        //             for (int jj = -kern_height/2; jj <= kern_height/2; ++jj) {
        //               for (int ii = -kern_width/2; ii <= kern_width/2; ++ii) {
        //                 float xx = x_base + d[0] * ii + d[1] * jj + offset(0);
        //                 float yy = y_base + d[3] * ii + d[4] * jj + offset(1);
        //                 right_image_patch(ii+kern_width/2, jj+kern_width/2) = right_interp_image(xx,yy);
        //               }
        //             }
        //             std::ostringstream ostr;
        //             ostr << x << "_" << y << "-" << iter;
        //             write_image("small/left-"+ostr.str()+".tif", left_image_patch);
        //             write_image("small/right-"+ostr.str()+".tif", right_image_patch);
        //             write_image("small/weight-"+ostr.str()+".tif", w);
        //           }
        

        // Solves lhs = rhs * x, and stores the result in-place in lhs.
        Matrix<double,6,6> pre_rhs = rhs;
        Vector<double,6> pre_lhs = lhs;
        try { 
          solve_symmetric_nocopy(rhs,lhs);
        } catch (ArgumentErr &e) {
          std::cout << "Error @ " << x << " " << y << "\n";
          //             std::cout << "Exception caught: " << e.what() << "\n";
          //             std::cout << "PRERHS: " << pre_rhs << "\n";
          //             std::cout << "PRELHS: " << pre_lhs << "\n\n";
          //             std::cout << "RHS: " << rhs << "\n";
          //             std::cout << "LHS: " << lhs << "\n\n";
          //             std::cout << "DEBUG: " << rhs(0,1) << "   " << rhs(1,0) << "\n\n";
          //             exit(0);
        }
        d += lhs;

        //          std::cout << "Update: " << lhs << "     " << d << "     " << sqrt(error_total) << "    " << (sqrt(lhs[2]*lhs[2]+lhs[5]*lhs[5])) << "\n";

        // Termination condition
        if (norm_2(lhs) < 0.01) 
          break;
      }
      //        std::cout << "----> " << d << "\n\n";
      offset(0) = d[2];
      offset(1) = d[5];
      
      if (norm_2(offset) > 1.5 || 
          offset(0) != offset(0) ||  // Check to make sure the offset is not NaN...
          offset(1) != offset(1) ) { // ... ditto.
        return PixelDisparity<float>();
      } else {
        PixelDisparity<float> returnval = m_disparity_map(x,y);
        returnval.h() += offset(0);
        returnval.v() += offset(1);
        return returnval;
      }
    }

    void cache(BBox2i bbox) {
      //      std::cout << "bbox: " << bbox << "\n";
      m_left_bbox = bbox;
      m_left_bbox.min() -= Vector2i(m_kern_width/2+1,m_kern_height/2+1);
      m_left_bbox.max() += Vector2i(m_kern_width/2+1,m_kern_height/2+1);
      //      std::cout << "left bbox: " << m_left_bbox << "\n";

      ImageView<typename InputViewT::pixel_type> left_buf = crop( m_left_log_image, m_left_bbox );
      ImageView<typename InputViewT::pixel_type> x_deriv_buf = derivative_filter(crop( m_left_log_image, m_left_bbox ), 1, 0);
      ImageView<typename InputViewT::pixel_type> y_deriv_buf = derivative_filter(crop( m_left_log_image, m_left_bbox ), 0, 1);
      m_left_cached_log_image.reset(new CropView<ImageView<float> >( left_buf,
                                                                     BBox2i(-m_left_bbox.min().x(), -m_left_bbox.min().y(),
                                                                            m_left_image.cols(), m_left_image.rows()) ) );

      m_x_deriv.reset(new CropView<ImageView<float> >( x_deriv_buf,
                                                       BBox2i(-m_left_bbox.min().x(), -m_left_bbox.min().y(),
                                                              m_left_image.cols(), m_left_image.rows()) ) );

      m_y_deriv.reset(new CropView<ImageView<float> >( y_deriv_buf,
                                                       BBox2i(-m_left_bbox.min().x(), -m_left_bbox.min().y(),
                                                              m_left_image.cols(), m_left_image.rows()) ) );

      int num_good;
      BBox2 disp_range = disparity::get_disparity_range(crop(edge_extend(m_disparity_map,ZeroEdgeExtension()), m_left_bbox), num_good, false);
      //      std::cout << "disparity range: " << disp_range << "\n";

      m_right_bbox = bbox;
      m_right_bbox.min() -= Vector2i(m_kern_width/2+1,m_kern_height/2+1);
      m_right_bbox.min() += disp_range.min();
      m_right_bbox.max() += Vector2i(m_kern_width/2+1,m_kern_height/2+1);
      m_right_bbox.max() += disp_range.max();
      //      std::cout << "Right bbox: " << m_right_bbox << "\n";

      ImageView<typename InputViewT::pixel_type> right_buf = crop( m_right_log_image, m_right_bbox );
      m_right_cached_log_image.reset(new CropView<ImageView<float> >(right_buf, BBox2i(-m_right_bbox.min().x(), -m_right_bbox.min().y(),
                                                                                       m_right_image.cols(), m_right_image.rows()) ) );

    }

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
