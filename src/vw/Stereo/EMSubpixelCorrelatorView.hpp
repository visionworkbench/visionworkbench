#include "EMSubpixelCorrelatorView.h"
#include "AffineMixtureComponent.h"
#include "GaussianMixtureComponent.h"
#include "GammaMixtureComponent.h"
#include "UniformMixtureComponent.h"

#include <vw/Math.h>
#include <vw/Image.h>
#include <iostream>

//#define USE_GRAPHICS

#ifdef USE_GRAPHICS
#include "graphics.h"
#endif

using namespace std;

#ifndef __VW_STEREO_EM_SUBPIXEL_CORRELATOR_VIEW_CPP__
#define __VW_STEREO_EM_SUBPIXEL_CORRELATOR_VIEW_CPP__


namespace vw {
  template<> struct PixelFormatID<Vector2>   { static const PixelFormatEnum value = VW_PIXEL_GENERIC_2_CHANNEL; };
  
  namespace stereo {
    inline ImageView<float> compute_spatial_weight_image(int kern_width, int kern_height, float two_sigma_sqr) {
  
      int center_pix_x = kern_width/2;
      int center_pix_y = kern_height/2;
      float sum;
      
      sum = 0.0;
      ImageView<float> weight(kern_width, kern_height);
      for (int j = 0; j < kern_height; ++j) {
	for (int i = 0; i < kern_width; ++i ) {
	  weight(i,j) = exp(-1 * (pow(i-center_pix_x,2) + pow(j-center_pix_y,2)) / two_sigma_sqr);
	  sum = sum + weight(i,j);
	}
      }
   
      for (int j = 0; j < kern_height; ++j) {
	for (int i = 0; i < kern_width; ++i ) {
	  weight(i,j) = weight(i, j)/sum;
	}
      }
      
      return weight;
    }
    
    template<class WeightT, class DisparityT>
    inline int adjust_weight_image(ImageViewBase<WeightT> &weight,
				   ImageViewBase<DisparityT> const& disparity_map_patch,
				   ImageViewBase<WeightT> const& weight_template) {
    
      //    const float continuity_threshold_squared = 64;  // T = 8
      int center_pix_x = weight_template.impl().cols()/2;
      int center_pix_y = weight_template.impl().rows()/2;
      typename DisparityT::pixel_type center_pix = disparity_map_patch.impl()(center_pix_x, center_pix_y);

      float sum = 0;
      int num_good_pix = 0;
      typename WeightT::pixel_accessor weight_row_acc = weight.impl().origin();
      typename WeightT::pixel_accessor template_row_acc = weight_template.impl().origin();
      typename DisparityT::pixel_accessor disp_row_acc = disparity_map_patch.impl().origin();
      
      for (int j = 0; j < weight_template.impl().rows(); ++j) {
	typename WeightT::pixel_accessor weight_col_acc = weight_row_acc;
	typename WeightT::pixel_accessor template_col_acc = template_row_acc;
	typename DisparityT::pixel_accessor disp_col_acc = disp_row_acc;
	for (int i = 0; i < weight_template.impl().cols(); ++i ) {

	  // Mask is zero if the disparity map's pixel is missing...
	  if (!(*disp_col_acc).valid()) 
	    *weight_col_acc = 0;

	  // ... or if there is a large discontinuity ...
	  //         if (pow( (*disp_col_acc).child().x()-center_pix.child().x(),2) + pow( (*disp_col_acc).child().y()-center_pix.child().y(),2) >= continuity_threshold_squared)
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


    // EMSubpixelCorrelator::EMSubpixelCorrelator( ... )
    template<class ImagePixelT, class DisparityPixelT, class PreProcFuncT>
    template <class ImageT, class DisparityT>
    EMSubpixelCorrelatorView<ImagePixelT, DisparityPixelT, PreProcFuncT >::EMSubpixelCorrelatorView(ImageViewBase<ImageT> const& left_image, ImageViewBase<ImageT> const& right_image,
												    ImageViewBase<DisparityT> const& course_disparity, PreProcFuncT preproc_func) :
      m_left_image(left_image.impl()), m_right_image(right_image.impl()), 
      m_course_disparity(course_disparity.impl()), m_preproc_func(preproc_func)
    {
      // Basic assertions
      VW_ASSERT((left_image.impl().cols() == right_image.impl().cols()) &&
                (left_image.impl().rows() == right_image.impl().rows()),
                ArgumentErr() << "EMSubpixelCorrelatorView::EMSubpixelCorrelatorView(): input image dimensions do not agree.\n");
      
      VW_ASSERT((left_image.impl().cols() == course_disparity.impl().cols()) &&
		(left_image.impl().rows() == course_disparity.impl().rows()),
                ArgumentErr() << "EMSubpixelCorrelatorView::EMSubpixelCorrelatorView(): input image and course disparity map dimensions do not agree.\n");
      
      VW_ASSERT((left_image.impl().channels() == 1) && (left_image.impl().planes() == 1) &&
                (right_image.impl().channels() == 1) && (right_image.impl().planes() == 1),
                ArgumentErr() << "EMSubpixelCorrelatorView::EMSubpixelCorrelatorView(): multi-channel, multi-plane images not supported.\n");
      
      // Set some sensible default values
      m_search_range = BBox2i(-50,-50,100,100);
      m_kernel_size = Vector2i(25, 25);
      m_cross_corr_threshold = 2.0;
      m_corr_score_threshold = 1.3;
      m_cost_blur = 1;
      m_correlator_type = ABS_DIFF_CORRELATOR;  

      r_window_dx = ImageView<ImagePixelT>(m_kernel_size[0], m_kernel_size[1]);
      r_window_dy = ImageView<ImagePixelT>(m_kernel_size[0], m_kernel_size[1]);
      r_window_dxdx = ImageView<ImagePixelT>(m_kernel_size[0], m_kernel_size[1]);
      r_window_dydy = ImageView<ImagePixelT>(m_kernel_size[0], m_kernel_size[1]);
      r_window_dxdy = ImageView<ImagePixelT>(m_kernel_size[0], m_kernel_size[1]);
    }

    // prerasterize( ... ) const;
    template <class ImagePixelT, class DisparityPixelT, class PreProcFuncT>
    typename EMSubpixelCorrelatorView<ImagePixelT, DisparityPixelT, PreProcFuncT>::prerasterize_type
    EMSubpixelCorrelatorView<ImagePixelT, DisparityPixelT, PreProcFuncT>::prerasterize(BBox2i bbox) const {      
      vw_out(InfoMessage, "stereo") << "EMSubpixelCorrelatorView: rasterizing image block " << bbox << ".\n";
      bool subpixel_debug = false;
      
      // Find the range of disparity values for this patch.
      int num_good;
      BBox2i search_range;
      try {
        search_range = get_disparity_range(crop(m_course_disparity, bbox));
      } 
      catch ( std::exception &e ) {
        search_range = BBox2i();
      }      

      ImageWindow window;
#ifdef USE_GRAPHICS      
      if(subpixel_debug) {
	window = create_window("disparity");
      }
#endif
      
      // The area in the right image that we'll be searching is
      // determined by the bbox of the left image plus the search
      // range.
      BBox2i left_crop_bbox(bbox);
      BBox2i right_crop_bbox(bbox.min() + search_range.min(),
                             bbox.max() + search_range.max());
      
      // The correlator requires the images to be the same size. The
      // search bbox will always be larger than the given left image
      // bbox, so we just make the left bbox the same size as the
      // right bbox.
      left_crop_bbox.max() = left_crop_bbox.min() + Vector2i(right_crop_bbox.width(), right_crop_bbox.height());

      // Finally, we must adjust both bounding boxes to account for
      // the size of the kernel itself.
      right_crop_bbox.min() -= Vector2i(m_kernel_size[0], m_kernel_size[1]);
      right_crop_bbox.max() += Vector2i(m_kernel_size[0], m_kernel_size[1]);
      left_crop_bbox.min() -= Vector2i(m_kernel_size[0], m_kernel_size[1]);
      left_crop_bbox.max() += Vector2i(m_kernel_size[0], m_kernel_size[1]);
      
      // We crop the images to the expanded bounding box and edge
      // extend in case the new bbox extends past the image bounds.
      ImageView<ImagePixelT> left_image_patch, right_image_patch;      
      ImageView<pixel_type> disparity_map_patch;
      
      left_image_patch = crop(edge_extend(m_left_image, ZeroEdgeExtension()),
			      left_crop_bbox);
      right_image_patch = crop(edge_extend(m_right_image, ZeroEdgeExtension()), 
			       right_crop_bbox);
      disparity_map_patch = pixel_cast<pixel_type>(crop(edge_extend(m_course_disparity, ZeroEdgeExtension()), 
							left_crop_bbox));
      
      // Adjust the disparities to be relative to the cropped
      // image pixel locations
      for (int v = 0; v < disparity_map_patch.rows(); ++v) {
        for (int u = 0; u < disparity_map_patch.cols(); ++u) {
          if (disparity_map_patch(u,v).valid())  { 
            disparity_map_patch(u,v).child().x() -= search_range.min().x();
            disparity_map_patch(u,v).child().y() -= search_range.min().y();
          }
	}
      }
      
            
      int pyramid_levels = 3;
      double blur_sigma_progressive = .5; // 3*sigma = 1.5 pixels
      // create the pyramid first
      std::vector<ImageView<ImagePixelT> > left_pyramid(pyramid_levels), right_pyramid(pyramid_levels);
      std::vector<BBox2i> regions_of_interest(pyramid_levels);
      std::vector<ImageView<Matrix2x2> > warps(pyramid_levels);  
      
      
      left_pyramid[0] = channels_to_planes(left_image_patch);
      right_pyramid[0] = channels_to_planes(right_image_patch);
      regions_of_interest[0] = BBox2i(m_kernel_size[0], m_kernel_size[1],
				      bbox.width(),bbox.height());
      
      std::vector<ImageView<pixel_type> > disparity_map_pyramid(pyramid_levels);
      std::vector<ImageView<pixel_type> > disparity_map_upsampled(pyramid_levels);
      disparity_map_pyramid[0] = disparity_map_patch;
      
      // downsample the disparity map and the image pair     
      for (int i = 1; i < pyramid_levels; i++) {
	left_pyramid[i] = subsample(gaussian_filter(left_pyramid[i-1], blur_sigma_progressive), 2);
	right_pyramid[i] = subsample(gaussian_filter(right_pyramid[i-1], blur_sigma_progressive), 2);	
	
	disparity_map_pyramid[i] = subsample_disp_map_by_two(disparity_map_pyramid[i-1]);
	regions_of_interest[i] = BBox2i(regions_of_interest[i-1].min()/2, regions_of_interest[i-1].max()/2);
      }

#ifdef USE_GRAPHICS
      if(subpixel_debug) {
	for(int i = 0; i < pyramid_levels; i++) {
	  show_image(window, left_pyramid[i]);
	  usleep((int)(.2*1000*1000));
	}
      }
#endif
      
      ImageView<ImagePixelT> process_left_image = (left_pyramid[pyramid_levels-1]);
      ImageView<ImagePixelT> process_right_image = (right_pyramid[pyramid_levels-1]);

      warps[pyramid_levels-1].set_size(process_left_image.cols(), process_left_image.rows());
      for(int y = 0; y < warps[pyramid_levels-1].rows(); y++) {
	for(int x = 0; x < warps[pyramid_levels-1].cols(); x++) {
	  warps[pyramid_levels-1](x, y).set_identity();
	}
      }     
            
      vw_out() << "processing pyramid level " << pyramid_levels-1 << endl;
      m_subpixel_refine(edge_extend(process_left_image, ZeroEdgeExtension()), edge_extend(process_right_image, ZeroEdgeExtension()), 
			disparity_map_pyramid[pyramid_levels-1], warps[pyramid_levels-1],
			regions_of_interest[pyramid_levels-1], false , subpixel_debug && true);
      
      disparity_map_upsampled[pyramid_levels-1] = copy(disparity_map_pyramid[pyramid_levels-1]);
      
      
      for (int i = pyramid_levels-2; i>=0; i--){
	int up_width = left_pyramid[i].cols();
	int up_height = left_pyramid[i].rows();
	vw_out() << "processing pyramid level " << i << endl;
	warps[i] = copy(resize(warps[i+1], up_width , up_height, ConstantEdgeExtension(), NearestPixelInterpolation())); //linear interpolation here
	disparity_map_upsampled[i] = copy(upsample_disp_map_by_two(disparity_map_upsampled[i+1], up_width, up_height));
	
	if(i == 0) {
	  process_left_image = channels_to_planes(left_image_patch);
	  process_right_image = channels_to_planes(right_image_patch);
	}
	else {
	  process_left_image = (left_pyramid[i]);
	  process_right_image = (right_pyramid[i]);
	}
	
	m_subpixel_refine(edge_extend(process_left_image, ZeroEdgeExtension()), edge_extend(process_right_image, ZeroEdgeExtension()), 
			  disparity_map_upsampled[i], warps[i],
			  regions_of_interest[i], i == 0, subpixel_debug && i == 0);
#ifdef USE_GRAPHICS
	if(subpixel_debug) {
	  show_image(window, resize(.5 + select_plane(channels_to_planes(disparity_map_upsampled[i])/6., 0), 400, 400, ZeroEdgeExtension(), NearestPixelInterpolation()));
	  usleep(2*1000*1000);
	}
#endif
      }
      
      disparity_map_patch = disparity_map_upsampled[0];      
      
#ifdef USE_GRAPHICS      
      if(subpixel_debug) {
	show_image(window, .5 + select_plane(channels_to_planes(disparity_map_patch)/6., 0));
	usleep(10*1000*1000);
      }
#endif

      // Undo the above adjustment
      for (int v = 0; v < disparity_map_patch.rows(); ++v) {
        for (int u = 0; u < disparity_map_patch.cols(); ++u) {
          if (disparity_map_patch(u,v).valid())  {
            disparity_map_patch(u,v).child().x() += search_range.min().x();
            disparity_map_patch(u,v).child().y() += search_range.min().y();
          }
	}
      }

#ifdef USE_GRAPHICS      
      if(subpixel_debug) {
	destroy_window(window);
      }
#endif
      return crop(disparity_map_patch, BBox2i(m_kernel_size[0]-bbox.min().x(), 
					      m_kernel_size[1]-bbox.min().y(), 
					      m_left_image.cols(), 
					      m_left_image.rows() ));
    }


    
    
    // m_sub_pixel_refine
    /* this actually performs the subpixel refinement */
    template <class ImagePixelT, class DisparityPixelT, class PreProcFuncT>
    template <class ImageT, class DisparityT, class AffineT>
    inline void
    EMSubpixelCorrelatorView<ImagePixelT, DisparityPixelT, PreProcFuncT>::m_subpixel_refine(ImageViewBase<ImageT> const& left_image,
											    ImageViewBase<ImageT> const& right_image,
											    ImageViewBase<DisparityT> & disparity,
											    ImageViewBase<AffineT> & affine_warps,
											    BBox2i const& ROI, bool final, bool p_debug) const {
      // set up some constants that should probably be external to this method:
      int em_iter_max = 20; //20;
      double P_inlier_0 = .8;
      double P_inlier_min = 1e-5;
      double P_inlier_max = 1-1e-5;
      double sigma_p1_0 = sqrt(1e-8); // sqrt(1e-3) for apollo // sqrt(1e-8) worked best
      double sigma_p_min = 1e-5; // 1e-3 for apollo
      double sigma_n_0 = .25; // sqrt(1e-2) for apollo      
      double sigma_n_min = 1e-4; // 1e-3 for apollo      
      double mu_n_0 = 0.;
      
      // algorithm features to enable
      bool blur_posterior = false;
      double blur_posterior_sigma = 1.5;
      
      bool use_left_outliers = true;
      bool use_right_outliers = false;
      
      bool save_states = true;
            
      double epsilon_em = 1;
      double debug_view_mag = 1;


      ImageWindow r_window_p1_view = NULL;
      ImageWindow l_window_view = NULL;
      ImageWindow w_window_p1_view = NULL;
      ImageWindow w_window_o1_view = NULL;
      ImageWindow w_window_o2_view = NULL;
      ImageWindow errors_window_p1_view = NULL;
      // set up windows for debug views

#ifdef USE_GRAPHICS            
      if(p_debug) {
	r_window_p1_view = create_window("r_window_p1");
	l_window_view = create_window("l_window");  
	w_window_p1_view = create_window("weights_p1");  
	if(use_left_outliers) {
	  w_window_o1_view = create_window("weights_outlier1");  
	}
	if(use_right_outliers) {
	  w_window_o2_view = create_window("weights_outlier2");  
	}
	errors_window_p1_view = create_window("errors_p1");        
      }
#endif
      
      //TODO: this debug region should be passed in rather then set here
      //BBox2i debug_region(45, 45, 1000, 1000); // in left ref
      //BBox2i debug_region(40, 30, 1000, 1000); // in right ref
      //BBox2i debug_region(17, 12, 5, 5);
      //BBox2i debug_region(100, 30, 25, 25);
      //BBox2i debug_region(50, 25, 1000, 1000);
      BBox2i debug_region(15, 30, 1, 1);
      /*
      double region_scale = 20;      
      BBox2i debug_region(ROI.center() - Vector2i(ROI.width()/region_scale, ROI.height()/region_scale),
			  ROI.center() + Vector2i(ROI.width()/region_scale, ROI.height()/region_scale));
      */
      
#ifdef USE_GRAPHICS      
      if(p_debug) {
	show_image(l_window_view, debug_view_mag*left_image.impl());
	show_image(r_window_p1_view, debug_view_mag*right_image.impl());
	usleep(1*1000*1000);
      }
#endif
      
      // allocate memory once for all inner loop operations
      ImageView<double> weights_p1(m_kernel_size(0), m_kernel_size(1));
      ImageView<double> weights_outlier1(m_kernel_size(0), m_kernel_size(1));
      ImageView<double> weights_outlier2(m_kernel_size(0), m_kernel_size(1));
      ImageView<double> temp_sum(m_kernel_size(0), m_kernel_size(1));
      
      AffineMixtureComponent<typename ImageT::pixel_type, double> affine_comp(laplacian_filter(gaussian_filter(left_image.impl(), 1.)), // 1.0 worked best
									      laplacian_filter(gaussian_filter(right_image.impl(), 1.)),
									      m_kernel_size, sigma_p1_0, sigma_p_min);
      //GaussianMixtureComponent<typename ImageT::pixel_type, double> outlier_comp1(left_image.impl(), m_kernel_size, mu_n_0, sigma_n_0, sigma_n_min);
      //GaussianMixtureComponent<typename ImageT::pixel_type, double> outlier_comp2(left_image.impl(), m_kernel_size, mu_n_0, sigma_n_0, sigma_n_min);
      GammaMixtureComponent<typename ImageT::pixel_type, double> outlier_comp1(left_image.impl(), m_kernel_size, 1., .25, 1e-2, 1e-2); // k_0 = 1., theta_0 = .25 worked best
      ImageView<float> image_hack = 1. - left_image.impl(); //TODO: remove this hack
      GammaMixtureComponent<typename ImageT::pixel_type, double> outlier_comp2(image_hack.impl(), m_kernel_size, 1., .25, 1e-2, 1e-2);
      
      int kernel_width = m_kernel_size(0);
      int kernel_height = m_kernel_size(1);
      
      ImageView<ImagePixelT> l_window(kernel_width, kernel_height);
      ImageView<ImagePixelT> r_window_p1(kernel_width, kernel_height);
      
      Vector2 pos;
      Vector2 cor_pos;

      float  two_sigma_sqr = 2.0*pow(float(m_kernel_size(0))/5.0,2.0); //4.0//7.0 works well // using 5.
      ImageView<float> weight_template = compute_spatial_weight_image(m_kernel_size(0), m_kernel_size(1), two_sigma_sqr);
      ImageView<float> w_gaussian(m_kernel_size(0), m_kernel_size(1));
      
      int N = m_kernel_size(0)*m_kernel_size(1); // number of pixels in each window      

      // loop through all pixels
      int x, y;
      int num_pixels = 0;
      Stopwatch pixel_timer;
      bool debug = false;
      
      for(y = ROI.min()[1]; y < ROI.max()[1]; y++) {
	if(debug) {
	  cout << "Processing row " << y << endl;
	}
	if(y%10 == 0 && y != 0) {
	  vw_out() << "@ row " << y << ": average pixel took " 
		   << 1000*pixel_timer.elapsed_seconds()/(double)num_pixels 
		   << " over " <<  num_pixels << " pixels" << endl;
	}
	for(x = ROI.min()[0]; x < ROI.max()[0]; x++) {
	  if(p_debug) {
	    if(x - ROI.min()[0] >= debug_region.min()[0] && 
	       x - ROI.min()[0] < debug_region.max()[0] && 
	       y - ROI.min()[1] >= debug_region.min()[1] && 
	       y - ROI.min()[1] < debug_region.max()[1]) {
	      debug = true;
	    }
	    else {
	      debug = false;
	    }
	  }	  
	  pos = Vector2(x, y);	  
	  
	  if(!disparity.impl()(x, y).valid()) { // skip missing pixels in the course map
	    continue;
	  }
	  
	  // set up the current window and the left image cropped view
	  Vector2 window_center(x, y);
	  BBox2i window_box = BBox2i(window_center - elem_diff(m_kernel_size, 1)/2, window_center + elem_sum(m_kernel_size, 1)/2);
	  if(window_box.min()(0) < 0 || window_box.min()(1) < 0) {
	    continue;
	  }
	  
	  //adjust_weight_image(w_gaussian,  crop(edge_extend(disparity.impl()), window_box), weight_template);
	  w_gaussian = weight_template;
	  	  
	  pixel_timer.start();
	  l_window = crop(edge_extend(left_image.impl(), ZeroEdgeExtension()), window_box);
	  
	  if(debug) {
	    cout << "course estimate = " << pos << " + " << disparity.impl()(x, y) << endl;
	  }
	  if(debug) {
	    cout << "initial warp: " <<affine_warps.impl()(x, y) << endl;
	  }
	  
	  affine_comp.reset(window_box, edge_extend(disparity.impl())(x, y).child().x(), edge_extend(disparity.impl())(x, y).child().y(), 
			    affine_warps.impl()(x, y),
			    w_gaussian);
	  outlier_comp1.reset(window_box, w_gaussian);
	  outlier_comp2.reset(window_box, w_gaussian);
	  	  
	  double P_1;
	  double P_outlier1;
	  double P_outlier2;	  
	  
	  if(use_left_outliers || use_right_outliers) {
	    P_1 = P_inlier_0;
	  }
	  else {
	    P_1 = 1.;
	  }
	  
	  if(!use_left_outliers && !use_right_outliers) {  // don't use outliers
	    P_outlier1 = 0.;
	    P_outlier2 = 0.;
	  }
	  else if(use_left_outliers && !use_right_outliers) { // use left outliers only
	    P_outlier1 = (1 - P_inlier_0);
	    P_outlier2 = 0.;
	  }
	  else if(!use_left_outliers && use_right_outliers) { // use righ outliers only
	    P_outlier1 = 0.;
	    P_outlier2 = (1 - P_inlier_0);
	  }
	  else { // use both left and right outliers
	    P_outlier1 = .66666*(1 - P_inlier_0);
	    P_outlier2 = .33333*(1 - P_inlier_0);
	  }

	  double mu_n = mu_n_0;
	  double sigma_n = sigma_n_0;
	  double sigma_p1 = sigma_p1_0;
	  
	  // run EM
	  int em_iter;
	  double f_v1 = INFINITY;
	  double f_v1_last = INFINITY;
	  double f_v2 = INFINITY;
	  double f_v2_last = INFINITY;
	  double f_value = INFINITY;
	  double f_value_last = INFINITY;
	  Stopwatch em_timer;
	  
	  double sum_weights_p1 = 0;
	  double sum_weights_outlier1 = 0;
	  double sum_weights_outlier2 = 0;

	  double ll_w0, ll_w1, ll_w2;

	  if(debug) {
	    cout << "P_1 = " << P_1 << endl;
	    cout << "P_outlier1 = " << P_outlier1 << endl;
	    cout << "P_outlier2 = " << P_outlier2 << endl;
	    
	    affine_comp.print_status("affine.");
	    if(use_left_outliers) {
	      outlier_comp1.print_status("outlier1.");
	    }
	    if(use_right_outliers) {
		outlier_comp2.print_status("outlier2.");
	    }
	    cout << "iter = 0" << endl;
#ifdef USE_GRAPHICS
	    r_window_p1 = crop(transform(right_image.impl(), affine_comp.affine_transform()), window_box);	    
	    show_image(r_window_p1_view, debug_view_mag*resize(r_window_p1, 200, 200, ZeroEdgeExtension(), NearestPixelInterpolation())); //, NearestPixelInterpolation()));
	    show_image(l_window_view, debug_view_mag*resize(l_window, 200, 200, ZeroEdgeExtension(), NearestPixelInterpolation()));
	    show_image(w_window_p1_view, resize(affine_comp.weights(), 200, 200, ZeroEdgeExtension(), NearestPixelInterpolation()));
	    show_image(w_window_o1_view, resize(weights_outlier1, 200, 200, ZeroEdgeExtension(), NearestPixelInterpolation()));
	    show_image(w_window_o2_view, resize(weights_outlier2, 200, 200, ZeroEdgeExtension(), NearestPixelInterpolation()));
	    show_image(errors_window_p1_view, debug_view_mag*resize(normalize(affine_comp.errors()), 200, 200, ZeroEdgeExtension(), NearestPixelInterpolation()));
	    cout << "max_error = " << max_pixel_value(affine_comp.errors()) << endl;
#endif	    
	    usleep((int)(1*1000*1000));
	  }
	  
	  affine_comp.update_posterior();
	  outlier_comp1.update_posterior();
	  outlier_comp2.update_posterior();
	  em_timer.start();
	  for(em_iter = 0; em_iter < em_iter_max; em_iter++) {
	    Stopwatch inner_loop_init_timer;
	    inner_loop_init_timer.start();
	    // compute the weights with the posterior distribution (Q)
	    weights_p1 = copy(affine_comp.prob()*P_1);
	    if(blur_posterior) {
	      if(use_left_outliers) {
		weights_outlier1 = copy(outlier_comp1.prob()*P_outlier1);
		weights_outlier1 = gaussian_filter(weights_outlier1, blur_posterior_sigma);
	      }
	      if(use_right_outliers) {
		weights_outlier2 = copy(outlier_comp2.prob()*P_outlier2);
		weights_outlier2 = gaussian_filter(weights_outlier2, blur_posterior_sigma);
	      }
	    }
	    else {
	      if(use_left_outliers) {
		weights_outlier1 = copy(outlier_comp1.prob()*P_outlier1);
	      }
	      if(use_right_outliers) {
		weights_outlier2 = copy(outlier_comp2.prob()*P_outlier2);
	      }
	    }
	    
	    // normalize the weights
	    temp_sum = copy(weights_p1);
	    if(use_left_outliers) {
	      temp_sum += weights_outlier1;
	    }
	    if(use_right_outliers) {
	      temp_sum += weights_outlier2;
	    }
	    weights_p1 = weights_p1/temp_sum;
	    if(use_left_outliers) {
	      weights_outlier1 = weights_outlier1/temp_sum; 
	    }
	    if(use_right_outliers) {
	      weights_outlier2 = weights_outlier2/temp_sum; 
	    }
	    	    
	    // compute P(inlier) and P(outlier)
	    sum_weights_p1 = sum_of_pixel_values(weights_p1);
	    if(use_left_outliers) {
	      sum_weights_outlier1 = sum_of_pixel_values(weights_outlier1);
	    }
	    if(use_right_outliers) {
	      sum_weights_outlier2 = sum_of_pixel_values(weights_outlier2);
	    }
	    
	    P_1 = sum_weights_p1/(double)(N);
	    P_1 = std::min(P_1, P_inlier_max);
	    P_1 = std::max(P_1, P_inlier_min);   
	    
	    if(use_left_outliers) {
	      P_outlier1 = sum_weights_outlier1/(double)(N);
	      P_outlier1 = std::min(P_outlier1, P_inlier_max);
	      P_outlier1 = std::max(P_outlier1, P_inlier_min);
	    }
	    if(use_right_outliers) {
	      P_outlier2 = sum_weights_outlier2/(double)(N);
	      P_outlier2 = std::min(P_outlier2, P_inlier_max);
	      P_outlier2 = std::max(P_outlier2, P_inlier_min);
	    }
	    
	    // fit all the parameters (M-step)
	    affine_comp.weights() = copy(weights_p1);
	    if(sum_weights_p1 >= 1e-2) {
	      affine_comp.fit_parameters();
	    }
	    if(use_left_outliers) {
	      outlier_comp1.weights() = copy(weights_outlier1);
	      outlier_comp1.fit_parameters();
	    }
	    if(use_right_outliers) {
	      outlier_comp2.weights() = copy(weights_outlier2);
	      outlier_comp2.fit_parameters();
	    }
	    
	    // update the posterior distributions to the parameter values
	    // note that for outlier_comp2 the posterior distribution is computed after transforming the right image by the updated affine transform
	    affine_comp.update_posterior();
	    if(use_left_outliers) {
	      outlier_comp1.update_posterior();
	    }
	    if(use_right_outliers) {
	      //outlier_comp2.set_affine_transform(affine_comp.affine_transform());
	      outlier_comp2.update_posterior();
	    }
	    // update the negative expected-log-likelihood (log-loss)
	    f_value_last = f_value;
	    
	    f_value = affine_comp.log_likelihood() - sum_weights_p1*log(P_1);
	    if(use_left_outliers) {
	      f_value += outlier_comp1.log_likelihood() - sum_weights_outlier1*log(P_outlier1);
	    }
	    if(use_right_outliers) {
	      f_value += outlier_comp2.log_likelihood() - sum_weights_outlier2*log(P_outlier2);
	    }

	    if(debug && save_states) {
	      if(em_iter == 0) {
		std::stringstream fname;
		fname << "state_"<< x << "_" << y << "_left.tif";
		write_image(fname.str(), l_window);
	      }
	      std::stringstream fname;
	      fname << "state_" << x << "_" << y << "_" <<  em_iter << "_weights.tif";
	      write_image(fname.str(), affine_comp.weights());
	    }
	    
	    if(debug) {
	      cout << "P(center) = " << weights_p1(x - window_box.min().x(),  y - window_box.min().y()) << endl;
	      cout << "P_1 = " << P_1 << endl;
	      cout << "P_outlier1 = " << P_outlier1 << endl;
	      cout << "P_outlier2 = " << P_outlier2 << endl;
	      
	      affine_comp.print_status("affine.");
	      if(use_left_outliers) {
		outlier_comp1.print_status("outlier1.");
	      }
	      if(use_right_outliers) {
		outlier_comp2.print_status("outlier2.");
	      }
	      
	      cout << "f_value = " << f_value << endl;	      
	    }	    
	    
	    if(debug) {
	      cout << "iter = " << em_iter+1 << endl;

#ifdef USE_GRAPHICS
	      r_window_p1 = crop(transform(right_image.impl(), affine_comp.affine_transform()), window_box);	    

	      double norm = 1.; //std::max<double>(max_pixel_value(outlier_comp1.prob()), max_pixel_value(affine_comp.prob()));
	      cout << "normalization constant = " << norm << endl;
	      show_image(r_window_p1_view, debug_view_mag*resize(r_window_p1, 200, 200, ZeroEdgeExtension(), NearestPixelInterpolation())); //, NearestPixelInterpolation()));
	      show_image(l_window_view, debug_view_mag*resize(l_window, 200, 200, ZeroEdgeExtension(), NearestPixelInterpolation()));
	      show_image(w_window_p1_view, resize((affine_comp.weights()/norm), 200, 200, ZeroEdgeExtension(), NearestPixelInterpolation()));	      
	      if(use_left_outliers) {
		show_image(w_window_o1_view, resize((outlier_comp1.weights()/norm), 200, 200, ZeroEdgeExtension(), NearestPixelInterpolation()));
	      }
	      if(use_right_outliers) {
		show_image(w_window_o2_view, resize(weights_outlier2, 200, 200, ZeroEdgeExtension(), NearestPixelInterpolation()));
	      }
	      show_image(errors_window_p1_view, debug_view_mag*resize(20*affine_comp.errors(), 200, 200, ZeroEdgeExtension(), NearestPixelInterpolation()));
#endif	      
	      usleep((int)(.5*1000*1000));

	      
	      cout << endl;
	    }	    
	    
	    // check termination condition for em loop
	    if(fabs(f_value_last - f_value) < epsilon_em) {
	      break;
	    }
	  } // end em loop
	  em_timer.stop();
	  	    
	  //cout << "EM converged in " << em_iter << " iterations (" << 1000*em_timer.elapsed_seconds() << ")" <<  endl;
	  
	  pos = Vector2(x, y);	  
	  cor_pos = affine_comp.affine_transform().reverse(pos);
	  Vector2f delta = (cor_pos - pos);
	  
	  disparity.impl()(x, y).child().x() = delta.x();
       	  disparity.impl()(x, y).child().y() = delta.y();
	  //disparity.impl()(x, y).cov_h() = affine_comp.hessian()(4, 4);
	  //disparity.impl()(x, y).cov_hv() = affine_comp.hessian()(4, 5);
	  //disparity.impl()(x, y).cov_v() = affine_comp.hessian()(5, 5);	

	  affine_warps.impl()(x, y) = affine_comp.affine_transform_mat();

	  	  
	  if(debug) {	  
	    cout << "refined estimate = " << disparity.impl()(x, y) << endl;
	    cout << "refined warp: " <<affine_warps.impl()(x, y) << endl;
	  }
	  
	  if(final) {
	    double P_center = weights_p1(x - window_box.min().x(),  y - window_box.min().y());
	    if(P_1 <= .25) { // if most pixels are outliers, set this one as missing
	      disparity.impl()(x, y).invalidate();
	    }
	    else if(P_1 <= .5){ // if at least half are inliers, decide based on the actual center pixel weight
	      if(P_center <= .95) {
		disparity.impl()(x, y).invalidate();
	      }
	    }
	    else if(P_1 <= .95) {
	      if(P_center <= .75) {
		disparity.impl()(x, y).invalidate();
	      }
	    }
	  }
	  pixel_timer.stop();
	  num_pixels++;
	} // end x loop
      } // end y loop
      
#ifdef USE_GRAPHICS      
      if(p_debug) {
	destroy_window(r_window_p1_view);
	destroy_window(l_window_view);
	destroy_window(w_window_p1_view);
	if(use_left_outliers) {
	  destroy_window(w_window_o1_view);
	}
	if(use_right_outliers) {
	  destroy_window(w_window_o2_view);
	}
	destroy_window(errors_window_p1_view);
      }      
#endif
    }   
    
    
    
    // operator<<( ... )
    template <class ImagePixelT, class DisparityPixelT, class PreProcFuncT>
    inline std::ostream& operator<<( std::ostream& os, EMSubpixelCorrelatorView<ImagePixelT, DisparityPixelT, PreProcFuncT> const& view) {
      os << "------------------------- EMSubpixelCorrelatorView ----------------------\n";
      os << "\tsearch range: " << view.search_range() << "\n";
      os << "\tkernel size : " << view.kernel_size() << "\n";
      os << "\txcorr thresh: " << view.cross_corr_threshold() << "\n";
      os << "\tcost blur: " << view.cost_blur() << "\n";
      os << "\tcorrelator type: " << view.correlator_type() << "\n";
      os << "\tcorrscore rejection thresh: " << view.corr_score_threshold() << "\n";
      os << "---------------------------------------------------------------\n";
      return os;
    }



    
    //image subsampling by two
    template <class ImageT>
    ImageT subsample_img_by_two(ImageViewBase<ImageT> const& img) {
      
      //determine the new size of the image
      //      printf("img: orig_w = %d, orig_h = %d\n", img.cols(), img.rows()); 
      int new_width;
      int new_height;
    
      new_width = img.impl().cols()/2;
      new_height = img.impl().rows()/2;
	
      //      printf("interp img: w = %d, h = %d, new_w = %d, new_h = %d\n", img.cols(), img.rows(), new_width, new_height); 
	
      ImageT outImg(new_width, new_height, img.impl().planes());		
      ImageViewRef<typename ImageT::pixel_type>  interpImg = interpolate(img.impl());
      int32 i, j, p;
      
      for (p = 0; p < outImg.planes() ; p++) {
	for (i = 0; i < outImg.cols(); i++) {
	  for (j = 0; j < outImg.rows(); j++) {
              
	    outImg(i,j,p) = interpImg(2*i, 2*j, p);
	      
	  }
	}
      }
      return outImg;
    }

    // disparity map down-sampling by two
    template <class PixelT>
    ImageView<PixelT > subsample_disp_map_by_two(ImageView<PixelT> const& input_disp) {
      // determine the new size of the image 
      int new_width = input_disp.cols()/2;
      int new_height = input_disp.rows()/2;
	
      ImageView<PixelT> outDisp(new_width, new_height);	
      ImageViewRef<PixelT> disp = edge_extend(input_disp, ConstantEdgeExtension() );
	
      for (int j = 0; j < outDisp.rows(); j++) {  
	for (int i = 0; i < outDisp.cols(); i++) {
	    
	  int old_i = i*2;
	  int old_j = j*2;
	  typename PixelChannelType<PixelT>::type h = disp(old_i  , old_j  ).child().x() + 
	    disp(old_i  , old_j+1).child().x() + 
	    disp(old_i+1, old_j  ).child().x() + 
	    disp(old_i+1, old_j+1).child().x();
	  typename PixelChannelType<PixelT>::type v = disp(old_i  , old_j  ).child().y() + 
	    disp(old_i  , old_j+1).child().y() + 
	    disp(old_i+1, old_j  ).child().y() + 
	    disp(old_i+1, old_j+1).child().y();
	    
	  int num_valid = 0;
	  if (disp(old_i  , old_j  ).valid()) ++num_valid;
	  if (disp(old_i+1, old_j  ).valid()) ++num_valid;
	  if (disp(old_i  , old_j+1).valid()) ++num_valid;
	  if (disp(old_i+1, old_j+1).valid()) ++num_valid;
	    
	  if (num_valid == 0)
	    outDisp(i,j) = PixelT();
	  else 
	    outDisp(i,j) = PixelT(h/typename PixelChannelType<PixelT>::type(num_valid) / 2.0,
						 v/typename PixelChannelType<PixelT>::type(num_valid) / 2.0);
	}
      }
	
      return outDisp;
    }
    
    
    // disparity map up-sampling by two
    template <class PixelT>
    ImageView<PixelT > upsample_disp_map_by_two(ImageView<PixelT> const& input_disp, int up_width, int up_height) {      
      ImageView<PixelT> outDisp(up_width, up_height);	
      ImageViewRef<PixelT> disp = edge_extend(input_disp, ConstantEdgeExtension());	
      
      for (uint j = 0; j < (uint)outDisp.rows(); ++j) {
        for (uint i = 0; i < (uint)outDisp.cols(); ++i) {
          int x = math::impl::_floor(float(i)/2.0), y = math::impl::_floor(float(j)/2.0);

          if ( i%2 == 0 && j%2 == 0) {
            if(!disp(x, y).valid())
              outDisp(i,j) = PixelT();
            else 
              outDisp(i,j) = PixelT( 2 * disp(x,y).child().x(),
                                                    2 * disp(x,y).child().y() );
          }
          
          else if (j%2 == 0) {
            if (!disp(x,y).valid() && !disp(x+1,y).valid())
              outDisp(i,j) = PixelT();
            else {
              if ( disp(x,y).valid() && disp(x+1,y).valid() ) {
                outDisp(i,j) = PixelT( 2 * (0.5 * disp(x,y).child().x() + 0.5 * disp(x+1,y).child().x() ),
                                                      2 * (0.5 * disp(x,y).child().y() + 0.5 * disp(x+1,y).child().y() ));
              } else if (disp(x,y).valid() && !disp(x+1,y).valid() ) {
                outDisp(i,j) = PixelT( 2 * disp(x,y).child().x(),
                                                      2 * disp(x,y).child().y() );
              } else if (!disp(x,y).valid() && disp(x+1,y).valid() ) {
                outDisp(i,j) = PixelT( 2 * disp(x+1,y).child().x(),
                                                      2 * disp(x+1,y).child().y() );
              } 
            }
          }	 

          else if (i%2 == 0) {
            if (!disp(x,y).valid() && !disp(x,y+1).valid())
              outDisp(i,j) = PixelT();
            else {
              if ( disp(x,y).valid() && disp(x,y+1).valid() ) {
                outDisp(i,j) = PixelT( 2 * (0.5 * disp(x,y).child().x() + 0.5 * disp(x,y+1).child().x() ),
                                                      2 * (0.5 * disp(x,y).child().y() + 0.5 * disp(x,y+1).child().y() ));
              } else if (disp(x,y).valid() && !disp(x,y+1).valid() ) {
                outDisp(i,j) = PixelT( 2 * disp(x,y).child().x(),
                                                      2 * disp(x,y).child().y() );
              } else if (!disp(x,y).valid() && disp(x,y+1).valid() ) {
                outDisp(i,j) = PixelT( 2 * disp(x,y+1).child().x(),
                                                      2 * disp(x,y+1).child().y() );
              } 
            }
          }	 
          
          else {
            if ( disp(x,y).valid() && disp(x,y+1).valid() && disp(x+1,y).valid() && disp(x+1,y+1).valid() ) {

              // All good pixels
              float normx = float(i)/2.0-x, normy = float(j)/2.0-y, norm1mx = 1.0-normx, norm1my = 1.0-normy;
              outDisp(i,j) = PixelT( 2 * (disp(x  ,y  ).child().x() * norm1mx*norm1my + 
                                                         disp(x+1,y  ).child().x() * normx*norm1my + 
                                                         disp(x  ,y+1).child().x() * norm1mx*normy + 
                                                         disp(x+1,y+1).child().x() * normx*normy),
                                                    2 * (disp(x  ,y  ).child().y() * norm1mx*norm1my + 
                                                         disp(x+1,y  ).child().y() * normx*norm1my + 
                                                         disp(x  ,y+1).child().y() * norm1mx*normy + 
                                                         disp(x+1,y+1).child().y() * normx*normy) ); 
	      
            } else if ( !disp(x,y).valid() && !disp(x,y+1).valid() && !disp(x+1,y).valid() && !disp(x+1,y+1).valid() ) {

              // no good pixels
              outDisp(i, j) = PixelT();

            } else {
              // some good & some bad pixels
              //
              // We fudge things a little bit here by picking the
              // first good pixel.  This isn't exactly correct, but
              // it's close enough in the handful of cases where we
              // are near some missing pixels, and we just need an
              // approximately valid value.  The subpixel refinement
              // will correct any minor mistakes we introduce here.
              if ( disp(x,y).valid() )
                outDisp(i,j) = PixelT(2 * disp(x,y).child().x(),
						 2 * disp(x,y).child().y());
              else if ( disp(x+1,y).valid() )
                outDisp(i,j) = PixelT(2 * disp(x+1,y).child().x(),
						 2 * disp(x+1,y).child().y());
              else if ( disp(x,y+1).valid() )
                outDisp(i,j) = PixelT(2 * disp(x,y+1).child().x(),
						 2 * disp(x,y+1).child().y());
              else if ( disp(x+1,y+1).valid() )
                outDisp(i,j) = PixelT(2 * disp(x+1,y+1).child().x(),
						 2 * disp(x+1,y+1).child().y());
            }


          }


        }
      }
      return outDisp;
    }
  } // end namepsace stereo
} // end namespace vw


#endif
