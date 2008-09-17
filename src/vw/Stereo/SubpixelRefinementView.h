#ifndef __VW_STEREO_SUBPIXEL_REFINEMENT_VIEW__
#define __VW_STEREO_SUBPIXEL_REFINEMENT_VIEW__

#include <vw/Image/ImageView.h>
#include <vw/Stereo/DisparityMap.h>
#include <vw/Stereo/Correlate.h>


namespace vw {
namespace stereo {

  // We get a considerable speedup in our 2d subpixel correlation if
  // we go ahead and compute the pseudoinverse of the A matrix (where
  // each row in A is [ x^2 y^2 xy x y 1] (our 2d parabolic surface)
  // for the range of x = [-1:1] and y = [-1:1].
  static double pinvA_data[] = { 1.0/6,  1.0/6,  1.0/6, -1.0/3, -1.0/3, -1.0/3,  1.0/6,  1.0/6,  1.0/6,
                                 1.0/6, -1.0/3,  1.0/6,  1.0/6, -1.0/3,  1.0/6,  1.0/6, -1.0/3,  1.0/6,
                                 1.0/4,    0.0, -1.0/4,    0.0,    0.0,    0.0, -1.0/4,    0.0,  1.0/4,
                                 -1.0/6, -1.0/6, -1.0/6,    0.0,    0.0,   0.0,  1.0/6,  1.0/6,  1.0/6,
                                 -1.0/6,    0.0,  1.0/6, -1.0/6,    0.0, 1.0/6, -1.0/6,    0.0,  1.0/6,
                                 -1.0/9,  2.0/9, -1.0/9,  2.0/9,  5.0/9, 2.0/9, -1.0/9,  2.0/9, -1.0/9 }; 

  /// An image view for performing image correlation
  template<class InputViewT, class DisparityViewT>
  class SubpixelRefinementView : public ImageViewBase<SubpixelRefinementView<InputViewT, DisparityViewT> > {

    DisparityViewT m_disparity_map;
    InputViewT m_left_image; 
    InputViewT m_right_image;
    ImageViewRef<float> m_left_log_image; 
    ImageViewRef<float> m_right_log_image;

    boost::shared_ptr< CropView<ImageView<float> > > m_left_cached_log_image;
    boost::shared_ptr< CropView<ImageView<float> > > m_right_cached_log_image;

    // Settings
    int m_kern_width, m_kern_height;
    bool m_do_h_subpixel, m_do_v_subpixel;
    bool m_verbose;

    // Parabola fit
    vw::MatrixProxy<double,6,9> pinvA;


    double find_minimum(double lt, double mid, double rt) const {
      double a = (rt+lt)*0.5-mid;
      double b = (rt-lt)*0.5;
      return -b/(2.0*a);
    }

    /* 
     * Find the minimun of a 2d hyperbolic surface that is fit to the nine points 
     * around and including the peak in the disparity map.  This gives better 
     * subpixel resolution when both horizontal and vertical subpixel is requested.
     * 
     * The equation of the surface we are fitting is:
     *    z = ax^2 + by^2 + cxy + dx + ey + f
     */
    vw::Vector2 find_minimum_2d(Vector<double,9> const& points, Matrix<double,6,9> const& pinvA) const {
      
      vw::Vector2 offset;
      
      /* 
       * First, compute the parameters of the hyperbolic surface by fitting the nine points in 'points'
       * using a linear least squares fit.  This process is fairly fast, since we have already pre-computed
       * the inverse of the A matrix in Ax = b.
       */
      vw::Vector<double> x = pinvA * points;
      
      /* 
       * With these parameters, we have a closed form expression for the surface.  We compute the 
       * derivative, and find the point where the slope is zero.  This is our maximum.
       *
       * Max is at [x,y] where:
       *
       *   dz/dx = 2ax + cy + d = 0
       *   dz/dy = 2by + cx + e = 0
       * 
       * Of course, we optimize this computation a bit by unrolling it by hand beforehand.
       */
      double denom = 4 * x(0) * x(1) - (x(2) * x(2));
      
      offset(0) = ( x(2) * x(4) - 2 * x(1) * x(3) ) / denom;
      offset(1) = ( x(2) * x(3) - 2 * x(0) * x(4) ) / denom;
      
      return offset;
    }


  public:
      typedef PixelDisparity<float> pixel_type;
      typedef pixel_type result_type;
      typedef ProceduralPixelAccessor<SubpixelRefinementView> pixel_accessor;
      
    SubpixelRefinementView(DisparityViewT const& disparity_map,
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
                                           m_verbose(verbose), pinvA(pinvA_data) {
      // Basic assertions
      VW_ASSERT((left_image.impl().cols() == right_image.impl().cols()) &&
                (left_image.impl().rows() == right_image.impl().rows()) &&
                (disparity_map.impl().cols() == right_image.impl().cols()) &&
                (disparity_map.impl().cols() == right_image.impl().cols()),
                ArgumentErr() << "SubpixelRefinementView::SubpixelRefinementView(): input image dimensions and/or disparity_map dimensions do not agree.\n");
      
      VW_ASSERT((left_image.channels() == 1) && (left_image.impl().planes() == 1) &&
                (right_image.channels() == 1) && (right_image.impl().planes() == 1),
                ArgumentErr() << "SubpixelRefinementView::SubpixelRefinementView(): multi-channel, multi-plane images not supported.\n");

      m_left_log_image = laplacian_filter(gaussian_filter(channel_cast<float>(channels_to_planes(left_image.impl())),1.5));
      m_right_log_image = laplacian_filter(gaussian_filter(channel_cast<float>(channels_to_planes(right_image.impl())),1.5));
    }

    // Standard ImageView interface methods
    inline int32 cols() const { return m_disparity_map.cols(); }
    inline int32 rows() const { return m_disparity_map.rows(); }
    inline int32 planes() const { return 1; }
    
    inline pixel_accessor origin() const { return pixel_accessor( *this, 0, 0 ); }
    
    inline pixel_type operator()(double i, double j, int32 p = 0) const {
      
      // Bail out if no subpixel computation has been requested 
      if (!m_do_h_subpixel && !m_do_v_subpixel) return m_disparity_map(i,j);

      if (m_left_cached_log_image) {
        return get_cached_pixel(i,j,p);
      } else {
        return get_uncached_pixel(i,j,p);
      }      
    }

    
    inline pixel_type get_cached_pixel(double i, double j, int32 p = 0) const {

      if ( !m_disparity_map(i,j).missing() ) {
        int hdisp= (int)m_disparity_map(i,j).h();
        int vdisp= (int)m_disparity_map(i,j).v();
        
        double mid = compute_soad(*m_left_cached_log_image, *m_right_cached_log_image,
                                  j, i,
                                  hdisp,   vdisp,
                                  m_kern_width, m_kern_height);
        
        // If only horizontal subpixel resolution is requested 
        if (m_do_h_subpixel && !m_do_v_subpixel) {
          double lt= compute_soad(*m_left_cached_log_image, *m_right_cached_log_image,
                                  j, i,
                                  hdisp-1, vdisp,
                                  m_kern_width, m_kern_height);
          double rt= compute_soad(*m_left_cached_log_image, *m_right_cached_log_image,
                                  j, i,
                                  hdisp+1, vdisp,
                                  m_kern_width, m_kern_height);
          if ((mid <= lt && mid < rt) || (mid <= rt && mid < lt)) {
            PixelDisparity<float> returnval = m_disparity_map(i,j);
            returnval.h() += find_minimum(lt, mid, rt);
            return returnval;
          } else {
            return PixelDisparity<float>();
          }
        }
        
        // If only vertical subpixel resolution is requested 
        if (m_do_v_subpixel && !m_do_h_subpixel) {
            double up= compute_soad(*m_left_cached_log_image, *m_right_cached_log_image,
                                    j, i,
                                    hdisp, vdisp-1,
                                    m_kern_width, m_kern_height);
            double dn= compute_soad(*m_left_cached_log_image, *m_right_cached_log_image,
                                    j, i,
                                    hdisp, vdisp+1,
                                    m_kern_width, m_kern_height);
          if ((mid <= up && mid < dn) || (mid <= dn && mid < up)) {
            PixelDisparity<float> returnval = m_disparity_map(i,j);
            returnval.v() += find_minimum(up, mid, dn);
            return returnval;
          } else {
            return PixelDisparity<float>();
          }
        }
                
        // If both vertical and horizontal subpixel resolution is requested,
        // we try to fit a 2d hyperbolic surface using the 9 points surrounding the
        // peak SOAD value.  
        //
        // We place the soad values into a vector using the following indices
        // (i.e. index 4 is the max disparity value)
        // 
        //     0  3  6
        //     1  4  7
        //     2  5  8
        //
        if (m_do_v_subpixel && m_do_h_subpixel) {
          vw::Vector<double,9> points;
          points(0) = (double)compute_soad(*m_left_cached_log_image, *m_right_cached_log_image,
                                           j, i,
                                           hdisp-1, vdisp-1,
                                           m_kern_width, m_kern_height);
          points(1) = (double)compute_soad(*m_left_cached_log_image, *m_right_cached_log_image,
                                           j, i,
                                           hdisp-1, vdisp,
                                           m_kern_width, m_kern_height);
          points(2) = (double)compute_soad(*m_left_cached_log_image, *m_right_cached_log_image,
                                           j, i,
                                           hdisp-1, vdisp+1,
                                           m_kern_width, m_kern_height);
          points(3) = (double)compute_soad(*m_left_cached_log_image, *m_right_cached_log_image,
                                           j, i,
                                           hdisp, vdisp-1,
                                           m_kern_width, m_kern_height);
          points(4) = (double)mid;
          points(5) = (double)compute_soad(*m_left_cached_log_image, *m_right_cached_log_image,
                                           j, i,
                                           hdisp, vdisp+1,
                                           m_kern_width, m_kern_height);
          points(6) = (double)compute_soad(*m_left_cached_log_image, *m_right_cached_log_image,
                                           j, i,
                                           hdisp+1, vdisp-1,
                                           m_kern_width, m_kern_height);
          points(7) = (double)compute_soad(*m_left_cached_log_image, *m_right_cached_log_image,
                                           j, i,
                                           hdisp+1, vdisp,
                                           m_kern_width, m_kern_height);
          points(8) = (double)compute_soad(*m_left_cached_log_image, *m_right_cached_log_image,
                                           j, i,
                                           hdisp+1, vdisp+1,
                                           m_kern_width, m_kern_height);
          
          vw::Vector2 offset = find_minimum_2d(points, pinvA);
          
          // This prevents us from adding in large offsets for
          // poorly fit data.
          if (fabs(offset(0)) < 2.0 && fabs(offset(1)) < 2.0) {
            PixelDisparity<float> returnval = m_disparity_map(i,j);
            returnval.h() += offset(0);
            returnval.v() += offset(1);
            return returnval;
          } else {
            return PixelDisparity<float>();
          }
        }
      }
      // Missing pixel
      return m_disparity_map(i,j);
    }

    inline pixel_type get_uncached_pixel(double i, double j, int32 p = 0) const {

      if ( !m_disparity_map(i,j).missing() ) {
        int hdisp= (int)m_disparity_map(i,j).h();
        int vdisp= (int)m_disparity_map(i,j).v();
        
        double mid = compute_soad(m_left_log_image, m_right_log_image,
                                  j, i,
                                  hdisp,   vdisp,
                                  m_kern_width, m_kern_height);
        
        // If only horizontal subpixel resolution is requested 
        if (m_do_h_subpixel && !m_do_v_subpixel) {
          double lt= compute_soad(m_left_log_image, m_right_log_image,
                                  j, i,
                                  hdisp-1, vdisp,
                                  m_kern_width, m_kern_height);
          double rt= compute_soad(m_left_log_image, m_right_log_image,
                                  j, i,
                                  hdisp+1, vdisp,
                                  m_kern_width, m_kern_height);
          if ((mid <= lt && mid < rt) || (mid <= rt && mid < lt)) {
            PixelDisparity<float> returnval = m_disparity_map(i,j);
            returnval.h() += find_minimum(lt, mid, rt);
            return returnval;
          } else {
            return PixelDisparity<float>();
          }
        }
        
        // If only vertical subpixel resolution is requested 
        if (m_do_v_subpixel && !m_do_h_subpixel) {
            double up= compute_soad(m_left_log_image, m_right_log_image,
                                    j, i,
                                    hdisp, vdisp-1,
                                    m_kern_width, m_kern_height);
            double dn= compute_soad(m_left_log_image, m_right_log_image,
                                    j, i,
                                    hdisp, vdisp+1,
                                    m_kern_width, m_kern_height);
          if ((mid <= up && mid < dn) || (mid <= dn && mid < up)) {
            PixelDisparity<float> returnval = m_disparity_map(i,j);
            returnval.v() += find_minimum(up, mid, dn);
            return returnval;
          } else {
            return PixelDisparity<float>();
          }
        }
                
        // If both vertical and horizontal subpixel resolution is requested,
        // we try to fit a 2d hyperbolic surface using the 9 points surrounding the
        // peak SOAD value.  
        //
        // We place the soad values into a vector using the following indices
        // (i.e. index 4 is the max disparity value)
        // 
        //     0  3  6
        //     1  4  7
        //     2  5  8
        //
        if (m_do_v_subpixel && m_do_h_subpixel) {
          vw::Vector<double,9> points;
          points(0) = (double)compute_soad(m_left_log_image, m_right_log_image,
                                           j, i,
                                           hdisp-1, vdisp-1,
                                           m_kern_width, m_kern_height);
          points(1) = (double)compute_soad(m_left_log_image, m_right_log_image,
                                           j, i,
                                           hdisp-1, vdisp,
                                           m_kern_width, m_kern_height);
          points(2) = (double)compute_soad(m_left_log_image, m_right_log_image,
                                           j, i,
                                           hdisp-1, vdisp+1,
                                           m_kern_width, m_kern_height);
          points(3) = (double)compute_soad(m_left_log_image, m_right_log_image,
                                           j, i,
                                           hdisp, vdisp-1,
                                           m_kern_width, m_kern_height);
          points(4) = (double)mid;
          points(5) = (double)compute_soad(m_left_log_image, m_right_log_image,
                                           j, i,
                                           hdisp, vdisp+1,
                                           m_kern_width, m_kern_height);
          points(6) = (double)compute_soad(m_left_log_image, m_right_log_image,
                                           j, i,
                                           hdisp+1, vdisp-1,
                                           m_kern_width, m_kern_height);
          points(7) = (double)compute_soad(m_left_log_image, m_right_log_image,
                                           j, i,
                                           hdisp+1, vdisp,
                                           m_kern_width, m_kern_height);
          points(8) = (double)compute_soad(m_left_log_image, m_right_log_image,
                                           j, i,
                                           hdisp+1, vdisp+1,
                                           m_kern_width, m_kern_height);
          
          vw::Vector2 offset = find_minimum_2d(points, pinvA);
          
          // This prevents us from adding in large offsets for
          // poorly fit data.
          if (fabs(offset(0)) < 2.0 && fabs(offset(1)) < 2.0) {
            PixelDisparity<float> returnval = m_disparity_map(i,j);
            returnval.h() += offset(0);
            returnval.v() += offset(1);
            return returnval;
          } else {
            return PixelDisparity<float>();
          }
        }
      }
      // Missing pixel
      return m_disparity_map(i,j);
    }

    void cache(BBox2i bbox) {
      BBox2i left_bbox = bbox;
      left_bbox.min() -= Vector2i(1,1);
      left_bbox.max() += Vector2i(1,1);

      ImageView<typename InputViewT::pixel_type> left_buf = crop( m_left_log_image, left_bbox );
      m_left_cached_log_image.reset(new CropView<ImageView<float> >(left_buf, BBox2i(-left_bbox.min().x(), -left_bbox.min().y(),
                                                                                     m_left_image.cols(), m_left_image.rows()) ) );

      int num_good;
      BBox2 disp_range = disparity::get_disparity_range(crop(edge_extend(m_disparity_map,ZeroEdgeExtension()), left_bbox), num_good, false);


      BBox2i right_bbox = bbox;
      right_bbox.min() -= Vector2i(1,1);
      right_bbox.min() += disp_range.min();
      right_bbox.max() += Vector2i(1,1);
      right_bbox.max() += disp_range.max();

      ImageView<typename InputViewT::pixel_type> right_buf = crop( m_right_log_image, right_bbox );
      m_right_cached_log_image.reset(new CropView<ImageView<float> >(right_buf, BBox2i(-right_bbox.min().x(), -right_bbox.min().y(),
                                                                                       m_right_image.cols(), m_right_image.rows()) ) );

    }

    /// \cond INTERNAL
    typedef SubpixelRefinementView prerasterize_type;
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
