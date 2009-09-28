// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifndef __VW_STEREO_EM_SUBPIXEL_CORRELATOR_VIEW__
#define __VW_STEREO_EM_SUBPIXEL_CORRELATOR_VIEW__

#include <vw/Image.h>
#include <vw/Image/ImageViewRef.h>
#include <vw/Stereo/Correlate.h>
#include <vw/Stereo/PyramidCorrelator.h>
#include <vw/Stereo/DisparityMap.h>
#include <vw/Math.h>
#include <string>
#include <ostream>




// For the PixelDisparity math.
#include <boost/operators.hpp>


namespace vw {
  namespace stereo {
    
    class AffineTransformOrigin;
    
    template <class ImageT>
      ImageT subsample_img_by_two(ImageViewBase<ImageT> const& img);

    template <class PixelT>
      ImageView<PixelT > subsample_disp_map_by_two(ImageView<PixelT> const& input_disp);
    template <class PixelT>
      ImageView<PixelT > upsample_disp_map_by_two(ImageView<PixelT> const& input_disp, int up_width, int up_height);



   
    template <class ImagePixelT, class DisparityPixelT, class PreProcFuncT>
      class EMSubpixelCorrelatorView : public ImageViewBase<EMSubpixelCorrelatorView<ImagePixelT, DisparityPixelT, PreProcFuncT> > {
    public:
      //typedef PixelUncertainDisparity<typename PixelChannelType<DisparityPixelT>::type> pixel_type;      
      typedef DisparityPixelT pixel_type;      
      typedef pixel_type result_type;
      typedef ProceduralPixelAccessor<EMSubpixelCorrelatorView> pixel_accessor;
      
      template <class ImageT, class DisparityT>
	EMSubpixelCorrelatorView(ImageViewBase<ImageT> const& left_image, ImageViewBase<ImageT> const& right_image,
				 ImageViewBase<DisparityT> const& course_disparity, PreProcFuncT preproc_func);
      
      
      
      // Basic accessor functions
      void set_search_range(BBox2i range) { m_search_range = range; }
      BBox2i search_range() const { return m_search_range; }
      
      void set_kernel_size(Vector2i size) { m_kernel_size = size; }
      Vector2i kernel_size() const { return m_kernel_size; }
      
      void set_correlator_options(int cost_blur, stereo::CorrelatorType correlator_type) {
	m_cost_blur = cost_blur;
	m_correlator_type = correlator_type;
      }
      int cost_blur() const { return m_cost_blur; }
      stereo::CorrelatorType correlator_type() const { return m_correlator_type; }
      
      void set_cross_corr_threshold(float threshold) { m_cross_corr_threshold = threshold; }
      float cross_corr_threshold() const { return m_cross_corr_threshold; }
      
      void set_corr_score_threshold(float threshold) { m_corr_score_threshold = threshold; }
      float corr_score_threshold() const { return m_corr_score_threshold; }
      
      /// Turn on debugging output.  The debug_file_prefix string is
      /// used as a prefix for all debug image files.
      void set_debug_mode(std::string const& debug_file_prefix) { m_debug_prefix = debug_file_prefix; }
      
      // Standard ImageView interface methods
      inline int32 cols() const { return m_left_image.cols(); }
      inline int32 rows() const { return m_left_image.rows(); }
      inline int32 planes() const { return 1; }
      
      inline pixel_accessor origin() const { return pixel_accessor( *this, 0, 0 ); }
      
      inline pixel_type operator()(double i, double j, int32 p = 0) const {
	vw_throw(NoImplErr() << "CorrelatorView::operator()(double i, double j, int32 p) has not been implemented.");
	return pixel_type();
      }
      
      
      /// \cond INTERNAL
      typedef CropView<ImageView<pixel_type> > prerasterize_type; 
      inline prerasterize_type prerasterize(BBox2i bbox) const;

      template <class DestT> inline void rasterize(DestT const& dest, BBox2i bbox) const {
	vw::rasterize(prerasterize(bbox), dest, bbox);
      }
      

      
    private: 
      // Image references
      ImageViewRef<ImagePixelT> m_left_image, m_right_image;
      ImageViewRef<DisparityPixelT> m_course_disparity;
      PreProcFuncT m_preproc_func;
      
      // Settings
      BBox2i m_search_range;
      Vector2i m_kernel_size;
      float m_cross_corr_threshold;
      float m_corr_score_threshold;
      int m_cost_blur;
      CorrelatorType m_correlator_type;
      std::string m_debug_prefix;


      // private helper methods
      template <class ImageT, class DisparityT, class AffineT>
	inline void 
	m_subpixel_refine(ImageViewBase<ImageT> const& left_image, ImageViewBase<ImageT> const& right_image,
			  ImageViewBase<DisparityT> &disparity, ImageViewBase<AffineT> & affine_warps,
			  BBox2i const& ROI, bool final, bool debug = false) const;
      
      inline void 
	m_compute_gradient_hessian(Vector<double, 6> &gradient, Matrix<double, 6, 6> &hessian, 
				   int x, int y,
				   ImageView<double> const &weights, ImageView<ImagePixelT> const &errors,
				   ImageView<ImagePixelT> const &r_window_dx, ImageView<ImagePixelT> const &r_window_dy,
				   ImageView<ImagePixelT> const &r_window_dxdx, ImageView<ImagePixelT> const &r_window_dydy,
				   ImageView<ImagePixelT> const &r_window_dxdy) const;

      template <class ImageT1, class ImageT2, class ImageT3, class DImageT, typename WeightT, typename ErrorT>
	inline double
      	m_fit_affine(Matrix<double, 6, 6> &hessian, 
		     double *matrix_data_linear,
		     double *matrix_data_offset,
		     AffineTransformOrigin &T,
		     MatrixProxy<double, 2, 2> const& M_transform_linear, 
		     VectorProxy<double, 2> const& M_transform_offset, 
		     ImageViewBase<ImageT1> &r_window, 
		     ImageViewBase<ImageT2> const& l_window, 
		     ImageViewBase<ImageT3> const& right_image, 
		     ImageViewBase<WeightT> const& weights, 
		     ImageViewBase<ErrorT> &errors, 
		     ImageViewBase<DImageT> const& r_image_dx, 
		     ImageViewBase<DImageT> const& r_image_dy,
		     ImageViewBase<DImageT> const& r_image_dxdx,
		     ImageViewBase<DImageT> const& r_image_dydy,
		     ImageViewBase<DImageT> const& r_image_dxdy,
		     double x, double y,
		     BBox2i const& window_box,
		     bool debug = false
		     ) const;
      // private temporary variables used in computation
      ImageView<ImagePixelT> r_window_dx;
      ImageView<ImagePixelT> r_window_dy;
      
      ImageView<ImagePixelT> r_window_dxdx;
      ImageView<ImagePixelT> r_window_dydy;
      ImageView<ImagePixelT> r_window_dxdy;

      
    }; // end class


    /// Summarize a CorrelatorView object
    template <class ImagePixelT, class DisparityPixelT, class PreProcFuncT>
      inline std::ostream& operator<<( std::ostream& os, EMSubpixelCorrelatorView<ImagePixelT, DisparityPixelT, PreProcFuncT> const& view );
    
  }
}


#include "EMSubpixelCorrelatorView.hpp"

#endif
