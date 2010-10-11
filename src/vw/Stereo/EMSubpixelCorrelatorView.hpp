// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifndef __VW_STEREO_EM_SUBPIXEL_CORRELATOR_VIEW_HPP__
#define __VW_STEREO_EM_SUBPIXEL_CORRELATOR_VIEW_HPP__

#include <vw/Stereo/AffineMixtureComponent.h>
#include <vw/Stereo/GammaMixtureComponent.h>
#include <vw/Stereo/DisparityMap.h>

#include <vw/Math.h>
#include <vw/FileIO/DiskImageResource.h>
#include <iostream>

#ifdef USE_GRAPHICS
#include <vw/gui/Graphics.h>
#endif

namespace vw {
  namespace stereo {

    namespace detail {
      ImageView<PixelMask<Vector2f> >
      subsample_disp_map_by_two(ImageView<PixelMask<Vector2f> > const& input_disp);

      ImageView<PixelMask<Vector2f> >
      upsample_disp_map_by_two(ImageView<PixelMask<Vector2f> > const& input_disp,
                               int up_width, int up_height);
    }


    inline ImageView<float> compute_spatial_weight_image(int kern_width, int kern_height, float two_sigma_sqr) {

      int center_pix_x = kern_width/2;
      int center_pix_y = kern_height/2;
      float sum;

      sum = 0.0;
      ImageView<float> weight(kern_width, kern_height);
      for (int j = 0; j < kern_height; ++j) {
        for (int i = 0; i < kern_width; ++i ) {
          weight(i,j) = exp(-1 * (pow(i-center_pix_x,2.) + pow(j-center_pix_y,2.)) / two_sigma_sqr);
          sum += weight(i,j);
        }
      }

      return weight / sum;
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
    template<class ImagePixelT>    template <class ImageT, class DisparityT>
    EMSubpixelCorrelatorView<ImagePixelT>::EMSubpixelCorrelatorView(ImageViewBase<ImageT> const& left_image, ImageViewBase<ImageT> const& right_image,
                                                                    ImageViewBase<DisparityT> const& course_disparity, int debug) :
      m_left_image(left_image.impl()), m_right_image(right_image.impl()),
      m_course_disparity(course_disparity.impl()), debug_level(debug)
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
      m_kernel_size = Vector2i(25, 25);
      pyramid_levels = 3;
      em_iter_max = 20; //20;
      P_inlier_0 = .8;
      P_inlier_min = 1e-5;
      P_inlier_max = 1-1e-5;
      sigma_p1_0 = sqrt(1e-8); // sqrt(1e-3) for apollo // sqrt(1e-8) worked best
      sigma_p1_min = 1e-5; // 1e-3 for apollo
      sigma_n_0 = .25; // sqrt(1e-2) for apollo
      sigma_n_min = 1e-4; // 1e-3 for apollo
      mu_n_0 = 0.;
      epsilon_em = 1;
      // affine model defaults
      inner_iter_max = 8; //8
      epsilon_inner = 1e-8; // 1e-8
      affine_min_det = .1; //.1
      affine_max_det = 1.9; //1.9
      debug_region = BBox2i(-1,-1, 0,0);
    }



    // prerasterize( ... ) const;
    template <class ImagePixelT>
    typename EMSubpixelCorrelatorView<ImagePixelT>::prerasterize_type
    EMSubpixelCorrelatorView<ImagePixelT>::prerasterize(BBox2i bbox) const {
      vw_out(InfoMessage, "stereo") << "EMSubpixelCorrelatorView: rasterizing image block " << bbox << ".\n";

      // Find the range of disparity values for this patch.
      // int num_good; // not used
      BBox2i search_range;
      try {
        search_range = get_disparity_range(crop(m_course_disparity, bbox));
      }
      catch ( std::exception &e ) {
        search_range = BBox2i();
      }


#ifdef USE_GRAPHICS
      ImageWindow window;
      if(debug_level >= 0) {
        window = vw_create_window("disparity");
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
      ImageView<disparity_pixel> disparity_map_patch_in;
      ImageView<result_type> disparity_map_patch_out;


      left_image_patch = crop(edge_extend(m_left_image, ZeroEdgeExtension()),
                              left_crop_bbox);
      right_image_patch = crop(edge_extend(m_right_image, ZeroEdgeExtension()),
                               right_crop_bbox);
      disparity_map_patch_in = crop(edge_extend(m_course_disparity, ZeroEdgeExtension()),
                                    left_crop_bbox);
      disparity_map_patch_out.set_size(disparity_map_patch_in.cols(), disparity_map_patch_in.rows());


      // Adjust the disparities to be relative to the cropped
      // image pixel locations
      for (int v = 0; v < disparity_map_patch_in.rows(); ++v) {
        for (int u = 0; u < disparity_map_patch_in.cols(); ++u) {
          if (disparity_map_patch_in(u,v).valid())  {
            disparity_map_patch_in(u,v).child().x() -= search_range.min().x();
            disparity_map_patch_in(u,v).child().y() -= search_range.min().y();
          }
        }
      }


      double blur_sigma_progressive = .5; // 3*sigma = 1.5 pixels

      // create the pyramid first
      std::vector<ImageView<ImagePixelT> > left_pyramid(pyramid_levels), right_pyramid(pyramid_levels);
      std::vector<BBox2i> regions_of_interest(pyramid_levels);
      std::vector<ImageView<Matrix2x2> > warps(pyramid_levels);
      std::vector<ImageView<disparity_pixel> > disparity_map_pyramid(pyramid_levels);


      // initialize the pyramid at level 0
      left_pyramid[0] = channels_to_planes(left_image_patch);
      right_pyramid[0] = channels_to_planes(right_image_patch);
      disparity_map_pyramid[0] = disparity_map_patch_in;
      regions_of_interest[0] = BBox2i(m_kernel_size[0], m_kernel_size[1],
                                      bbox.width(),bbox.height());


      // downsample the disparity map and the image pair to initialize the intermediate levels
      for(int i = 1; i < pyramid_levels; i++) {
        left_pyramid[i] = subsample(gaussian_filter(left_pyramid[i-1], blur_sigma_progressive), 2);
        right_pyramid[i] = subsample(gaussian_filter(right_pyramid[i-1], blur_sigma_progressive), 2);

        disparity_map_pyramid[i] = detail::subsample_disp_map_by_two(disparity_map_pyramid[i-1]);
        regions_of_interest[i] = BBox2i(regions_of_interest[i-1].min()/2, regions_of_interest[i-1].max()/2);
      }

      // initialize warps at the lowest resolution level
      warps[pyramid_levels-1].set_size(left_pyramid[pyramid_levels-1].cols(),
                                       left_pyramid[pyramid_levels-1].rows());
      for(int y = 0; y < warps[pyramid_levels-1].rows(); y++) {
        for(int x = 0; x < warps[pyramid_levels-1].cols(); x++) {
          warps[pyramid_levels-1](x, y).set_identity();
        }
      }

#ifdef USE_GRAPHICS
      vw_initialize_graphics(0, NULL);
      if(debug_level >= 0) {
        for(int i = 0; i < pyramid_levels; i++) {
          vw_show_image(window, left_pyramid[i]);
          usleep((int)(.2*1000*1000));
        }
      }
#endif

      // go up the pyramid; first run refinement, then upsample result for the next level
      for(int i = pyramid_levels-1; i >=0; i--) {
        vw_out() << "processing pyramid level "
                  << i << " of " << pyramid_levels-1 << std::endl;

        if(debug_level >= 0) {
          std::stringstream stream;
          stream << "pyramid_level_" << i << ".tif";
          write_image(stream.str(), disparity_map_pyramid[i]);
        }

        ImageView<ImagePixelT> process_left_image = left_pyramid[i];
        ImageView<ImagePixelT> process_right_image = right_pyramid[i];

        if(i > 0) { // in this case take refine the upsampled disparity map from the previous level,
          // and upsample for the next level
          m_subpixel_refine(edge_extend(process_left_image, ZeroEdgeExtension()), edge_extend(process_right_image, ZeroEdgeExtension()),
                            disparity_map_pyramid[i], disparity_map_pyramid[i], warps[i],
                            regions_of_interest[i], false, debug_level == i);

          // upsample the warps and the refined map for the next level of processing
          int up_width = left_pyramid[i-1].cols();
          int up_height = left_pyramid[i-1].rows();
          warps[i-1] = copy(resize(warps[i], up_width , up_height, ConstantEdgeExtension(), NearestPixelInterpolation())); //upsample affine transforms
          disparity_map_pyramid[i-1] = copy(detail::upsample_disp_map_by_two(disparity_map_pyramid[i], up_width, up_height));
        }
        else { // here there is no next level so we refine directly to the output patch
          m_subpixel_refine(edge_extend(process_left_image, ZeroEdgeExtension()), edge_extend(process_right_image, ZeroEdgeExtension()),
                            disparity_map_pyramid[i], disparity_map_patch_out, warps[i],
                            regions_of_interest[i], true, debug_level == i);
        }
      }

#ifdef USE_GRAPHICS
      if(debug_level >= 0) {
        vw_show_image(window, .5 + select_plane(channels_to_planes(disparity_map_patch_out)/6., 0));
        usleep(10*1000*1000);
      }
#endif

      // Undo the above adjustment
      for (int v = 0; v < disparity_map_patch_out.rows(); ++v) {
        for (int u = 0; u < disparity_map_patch_out.cols(); ++u) {
          if (disparity_map_patch_out(u,v).valid())  {
            disparity_map_patch_out(u,v).child().x() += search_range.min().x();
            disparity_map_patch_out(u,v).child().y() += search_range.min().y();
          }
        }
      }

#ifdef USE_GRAPHICS
      if(debug_level >= 0 ) {
        vw_destroy_window(window);
      }
#endif

      return crop(disparity_map_patch_out, BBox2i(m_kernel_size[0]-bbox.min().x(),
                                                  m_kernel_size[1]-bbox.min().y(),
                                                  m_left_image.cols(),
                                                  m_left_image.rows()));
    }




    // m_sub_pixel_refine
    /* this actually performs the subpixel refinement */
    template <class ImagePixelT>
    template <class ImageT, class DisparityT1, class DisparityT2, class AffineT>
    inline void
    EMSubpixelCorrelatorView<ImagePixelT>::m_subpixel_refine(ImageViewBase<ImageT> const& left_image,
                                                             ImageViewBase<ImageT> const& right_image,
                                                             ImageViewBase<DisparityT1> & disparity_in,
                                                             ImageViewBase<DisparityT2> & disparity_out,
                                                             ImageViewBase<AffineT> & affine_warps,
                                                             BBox2i const& ROI, bool final, bool p_debug) const {

      // algorithm features to enable
      bool blur_posterior = false;
      double blur_posterior_sigma = 1.5;

      bool use_left_outliers = true;
      bool use_right_outliers = false;

      bool save_states = false;


      // set up windows for debug views
#ifdef USE_GRAPHICS
      ImageWindow r_window_p1_view = NULL;
      ImageWindow l_window_view = NULL;
      ImageWindow w_window_p1_view = NULL;
      ImageWindow w_window_o1_view = NULL;
      ImageWindow w_window_o2_view = NULL;
      ImageWindow errors_window_p1_view = NULL;
      if(p_debug) {
        r_window_p1_view = vw_create_window("r_window_p1");
        l_window_view = vw_create_window("l_window");
        w_window_p1_view = vw_create_window("weights_p1");
        if(use_left_outliers) {
          w_window_o1_view = vw_create_window("weights_outlier1");
        }
        if(use_right_outliers) {
          w_window_o2_view = vw_create_window("weights_outlier2");
        }
        errors_window_p1_view = vw_create_window("errors_p1");
      }
#endif

#ifdef USE_GRAPHICS
      if(p_debug) {
        vw_show_image(l_window_view, debug_view_mag*left_image.impl());
        vw_show_image(r_window_p1_view, debug_view_mag*right_image.impl());
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
                                                                              m_kernel_size, sigma_p1_0, sigma_p1_min);
      affine_comp.set_max_iter(inner_iter_max);
      affine_comp.set_epsilon_convergence(epsilon_inner);
      affine_comp.set_min_determinant(affine_min_det);
      affine_comp.set_max_determinant(affine_max_det);
      //affine_comp.set_debug(true);

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
        if(y%10 == 0 && y != 0) {
          vw_out() << "@ row " << y << ": average pixel took "
                    << 1000*pixel_timer.elapsed_seconds()/(double)num_pixels
                    << " over " <<  num_pixels << " pixels" << std::endl;
        }
        for(x = ROI.min()[0]; x < ROI.max()[0]; x++) {
          pos = Vector2(x, y);

          if(p_debug) {
            if(debug_region.contains(pos - ROI.min())) {
              debug = true;
            }
            else {
              debug = false;
            }
          }

          if(!disparity_in.impl()(x, y).valid()) { // skip missing pixels in the course map
            disparity_out.impl()(x, y).invalidate();
            continue;
          }

          // set up the current window and the left image cropped view
          Vector2 window_center(x, y);
          BBox2i window_box = BBox2i(window_center - elem_diff(m_kernel_size, 1)/2, window_center + elem_sum(m_kernel_size, 1)/2);
          if(window_box.min()(0) < 0 || window_box.min()(1) < 0) {
            continue;
          }

          //adjust_weight_image(w_gaussian,  crop(edge_extend(disparity_in.impl()), window_box), weight_template);
          w_gaussian = weight_template;

          pixel_timer.start();
          l_window = crop(edge_extend(left_image.impl(), ZeroEdgeExtension()), window_box);

          if(debug) {
            vw_out() << "course estimate = " << pos << " + "
                      << disparity_in.impl()(x, y) << std::endl;
            vw_out() << "initial warp: "
                      << affine_warps.impl()(x, y) << std::endl;
          }

          affine_comp.reset(window_box, edge_extend(disparity_in.impl())(x, y).child().x(), edge_extend(disparity_in.impl())(x, y).child().y(),
                            affine_warps.impl()(x, y),
                            w_gaussian);
          outlier_comp1.reset(window_box, w_gaussian);
          outlier_comp2.reset(window_box, w_gaussian);

          double P_1;
          double P_outlier1;
          double P_outlier2;

          if (use_left_outliers || use_right_outliers)
            P_1 = P_inlier_0;

          else
            P_1 = 1.;

          if (!use_left_outliers && !use_right_outliers) {
            // don't use outliers
            P_outlier1 = 0.;
            P_outlier2 = 0.;
          } else if (use_left_outliers && !use_right_outliers) {
            // use left outliers only
            P_outlier1 = (1 - P_inlier_0);
            P_outlier2 = 0.;
          } else if (!use_left_outliers && use_right_outliers) {
            // use righ outliers only
            P_outlier1 = 0.;
            P_outlier2 = (1 - P_inlier_0);
          } else {
            // use both left and right outliers
            P_outlier1 = .66666*(1 - P_inlier_0);
            P_outlier2 = .33333*(1 - P_inlier_0);
          }

          // double mu_n = mu_n_0; // not used
          // double sigma_n = sigma_n_0; // not used
          // double sigma_p1 = sigma_p1_0; // not used

          // run EM
          int em_iter;
          // double f_v1 = INFINITY; // not used
          // double f_v1_last = INFINITY; // not used
          // double f_v2 = INFINITY; // not used
          // double f_v2_last = INFINITY; // not used
          double f_value = INFINITY;
          double f_value_last = INFINITY;
          Stopwatch em_timer;

          double sum_weights_p1 = 0;
          double sum_weights_outlier1 = 0;
          double sum_weights_outlier2 = 0;

          // double ll_w0, ll_w1, ll_w2; // not used

          if (debug) {
            vw_out() << "P_1 = " << P_1 << std::endl;
            vw_out() << "P_outlier1 = " << P_outlier1 << std::endl;
            vw_out() << "P_outlier2 = " << P_outlier2 << std::endl;

            affine_comp.print_status("affine.");
            if(use_left_outliers)
              outlier_comp1.print_status("outlier1.");
            if(use_right_outliers)
              outlier_comp2.print_status("outlier2.");
            vw_out() << "iter = 0" << std::endl;

            vw_out() << "RMS error = " << sqrt(sum_of_pixel_values(pow(affine_comp.errors(),2))/affine_comp.errors().cols()/affine_comp.errors().rows()) << std::endl;
#ifdef USE_GRAPHICS
            r_window_p1 = crop(transform(right_image.impl(), affine_comp.affine_transform()), window_box);
            vw_show_image(r_window_p1_view, debug_view_mag*resize(r_window_p1, 200, 200, ZeroEdgeExtension(), NearestPixelInterpolation())); //, NearestPixelInterpolation()));
            vw_show_image(l_window_view, debug_view_mag*resize(l_window, 200, 200, ZeroEdgeExtension(), NearestPixelInterpolation()));
            vw_show_image(w_window_p1_view, resize(affine_comp.weights(), 200, 200, ZeroEdgeExtension(), NearestPixelInterpolation()));
            vw_show_image(w_window_o1_view, resize(weights_outlier1, 200, 200, ZeroEdgeExtension(), NearestPixelInterpolation()));
            vw_show_image(w_window_o2_view, resize(weights_outlier2, 200, 200, ZeroEdgeExtension(), NearestPixelInterpolation()));
            vw_show_image(errors_window_p1_view, debug_view_mag*resize(normalize(affine_comp.errors()), 200, 200, ZeroEdgeExtension(), NearestPixelInterpolation()));
            vw_out() << "max_error = " << max_pixel_value(affine_comp.errors()) << std::endl;
            usleep((int)(1*1000*1000));
#endif

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
            if(use_left_outliers)
              temp_sum += weights_outlier1;
            if(use_right_outliers)
              temp_sum += weights_outlier2;
            weights_p1 = weights_p1/temp_sum;
            if(use_left_outliers)
              weights_outlier1 = weights_outlier1/temp_sum;
            if(use_right_outliers)
              weights_outlier2 = weights_outlier2/temp_sum;

            // compute P(inlier) and P(outlier)
            sum_weights_p1 = sum_of_pixel_values(weights_p1);
            if(use_left_outliers)
              sum_weights_outlier1 = sum_of_pixel_values(weights_outlier1);
            if(use_right_outliers)
              sum_weights_outlier2 = sum_of_pixel_values(weights_outlier2);

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
              vw_out() << "P(center) = " << weights_p1(x - window_box.min().x(),  y - window_box.min().y()) << std::endl;
              vw_out() << "P_1 = " << P_1 << std::endl;
              vw_out() << "P_outlier1 = " << P_outlier1 << std::endl;

              affine_comp.print_status("affine.");
              if(use_left_outliers)
                outlier_comp1.print_status("outlier1.");
              if(use_right_outliers)
                outlier_comp2.print_status("outlier2.");

              vw_out() << "f_value = " << f_value << std::endl;
            }

            if (debug) {
              vw_out() << "iter = " << em_iter+1 << std::endl;
              vw_out() << "RMS error = " << sqrt(sum_of_pixel_values(pow(affine_comp.errors(),2))/affine_comp.errors().cols()/affine_comp.errors().rows()) << std::endl;
#ifdef USE_GRAPHICS
              r_window_p1 = crop(transform(right_image.impl(), affine_comp.affine_transform()), window_box);

              double norm = 1.; //std::max<double>(max_pixel_value(outlier_comp1.prob()), max_pixel_value(affine_comp.prob()));
              vw_out() << "normalization constant = " << norm << std::endl;
              vw_show_image(r_window_p1_view, debug_view_mag*resize(r_window_p1, 200, 200, ZeroEdgeExtension(), NearestPixelInterpolation())); //, NearestPixelInterpolation()));
              vw_show_image(l_window_view, debug_view_mag*resize(l_window, 200, 200, ZeroEdgeExtension(), NearestPixelInterpolation()));
              vw_show_image(w_window_p1_view, resize((affine_comp.weights()/norm), 200, 200, ZeroEdgeExtension(), NearestPixelInterpolation()));
              if(use_left_outliers) {
                vw_show_image(w_window_o1_view, resize((outlier_comp1.weights()/norm), 200, 200, ZeroEdgeExtension(), NearestPixelInterpolation()));
              }
              if(use_right_outliers) {
                vw_show_image(w_window_o2_view, resize(weights_outlier2, 200, 200, ZeroEdgeExtension(), NearestPixelInterpolation()));
              }
              vw_show_image(errors_window_p1_view, debug_view_mag*resize(20*affine_comp.errors(), 200, 200, ZeroEdgeExtension(), NearestPixelInterpolation()));
              usleep((int)(.5*1000*1000));
#endif
              vw_out() << std::endl;
            }

            // check termination condition for em loop
            if(fabs(f_value_last - f_value) < epsilon_em) {
              break;
            }
          } // end em loop
          em_timer.stop();

          //vw_out() << "EM converged in " << em_iter << " iterations (" << 1000*em_timer.elapsed_seconds() << ")" <<  std::endl;

          pos = Vector2(x, y);
          cor_pos = affine_comp.affine_transform().reverse(pos);
          Vector2f delta = (cor_pos - pos);

          disparity_out.impl()(x, y).validate();
          disparity_out.impl()(x, y).child().x() = delta.x();
          disparity_out.impl()(x, y).child().y() = delta.y();
          if(disparity_out.impl()(x, y).child().size() == 5) { // if the disparity map given has space for uncertainty, output that as well
            disparity_out.impl()(x, y).child()(2) = affine_comp.hessian()(4, 4);
            disparity_out.impl()(x, y).child()(3) = affine_comp.hessian()(4, 5);
            disparity_out.impl()(x, y).child()(4) = affine_comp.hessian()(5, 5);
          }
          affine_warps.impl()(x, y) = affine_comp.affine_transform_mat();


          if(debug) {
            vw_out() << "refined estimate = "
                      << disparity_out.impl()(x, y) << std::endl;
            vw_out() << "refined warp: "
                      << affine_warps.impl()(x, y) << std::endl;
          }

          if(final) {
            double P_center = weights_p1(x - window_box.min().x(),  y - window_box.min().y());
            if(P_1 <= .25) { // if most pixels are outliers, set this one as missing
              disparity_out.impl()(x, y).invalidate();
            }
            else if(P_1 <= .5){ // if at least half are inliers, decide based on the actual center pixel weight
              if(P_center <= .95) {
                disparity_out.impl()(x, y).invalidate();
              }
            }
            else if(P_1 <= .95) {
              if(P_center <= .75) {
                disparity_out.impl()(x, y).invalidate();
              }
            }
          }
          pixel_timer.stop();
          num_pixels++;
        } // end x loop
      } // end y loop

#ifdef USE_GRAPHICS
      if(p_debug) {
        vw_destroy_window(r_window_p1_view);
        vw_destroy_window(l_window_view);
        vw_destroy_window(w_window_p1_view);
        if(use_left_outliers) {
          vw_destroy_window(w_window_o1_view);
        }
        if(use_right_outliers) {
          vw_destroy_window(w_window_o2_view);
        }
        vw_destroy_window(errors_window_p1_view);
      }
#endif
    }



    // operator<<( ... )
    template <class ImagePixelT>
    inline std::ostream& operator<<( std::ostream& os, EMSubpixelCorrelatorView<ImagePixelT> const& view) {
      os << "------------------------- EMSubpixelCorrelatorView ----------------------\n";
      os << "\tkernel size : " << view.kernel_size() << "\n";
      os << "---------------------------------------------------------------\n";
      return os;
    }

}} // namespace vw::stereo


#endif
