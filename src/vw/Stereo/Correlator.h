#ifndef __VW_STEREO_CORRELATOR_H__
#define __VW_STEREO_CORRELATOR_H__

#include <vw/Image/ImageView.h>
#include <vw/Stereo/DisparityMap.h>
#include <vw/Image/Transform.h>

namespace vw {
namespace stereo {

  class Correlator {

    BBox2i m_initial_search_range;
    Vector2i m_kernel_size;
    float m_slog_width;
    float m_cross_correlation_threshold;
    int m_pyramid_min_image_dimension;

    std::string m_debug_prefix;

    // Reduce the image size by a factor of two by averaging the pixels
    template <class PixelT>
    ImageView<PixelT> subsample_by_two(ImageView<PixelT> &img) {
      
      ImageView<PixelT> outImg(img.cols()/2, img.rows()/2,img.planes());		
      int32 i, j, p;
    
      for (p = 0; p < outImg.planes() ; p++) {
        for (i = 0; i < outImg.cols(); i++) {
          for (j = 0; j < outImg.rows(); j++) {  
            outImg(i,j,p) = 0.0f;
            outImg(i,j,p) += img(2*i     , 2*j    ,p);
            outImg(i,j,p) += img(2*i + 1 , 2*j    ,p);
            outImg(i,j,p) += img(2*i     , 2*j + 1,p);
            outImg(i,j,p) += img(2*i + 1 , 2*j + 1,p);
            outImg(i,j,p) /= 4;
          }
        }
      }
      return outImg;
    }

    // Iterate over the nominal blocks, creating output blocks for correlation
    BBox2i compute_matching_blocks(BBox2i const& nominal_block, BBox2i search_range, 
                                   BBox2i &left_block, BBox2i &right_block);

    std::vector<BBox2i> compute_search_ranges(ImageView<PixelDisparity<float> > const& prev_disparity_map, 
                                              std::vector<BBox2i> nominal_blocks);

    ImageView<PixelDisparity<float> > do_correlation(std::vector<ImageView<uint8> > left_slog_pyramid, 
                                                     std::vector<ImageView<uint8> > right_slog_pyramid,
                                                     std::vector<ImageView<bool> > left_masks, 
                                                     std::vector<ImageView<bool> > right_masks);

    ImageView<PixelDisparity<float> > correlate(ImageView<uint8> left_image, ImageView<uint8> right_image, 
                                                BBox2i search_range, Vector2i offset);
                                                

  public:

    /// Correlator Constructor
    ///
    /// Set pyramid_levels to 0 to force the use of a single pyramid level (essentially disabling pyramid correlation).
    Correlator(BBox2i initial_search_range, Vector2i kernel_size, int slog_width = 1.5, float cross_correlation_threshold, int pyramid_min_image_dimension = 256) :
      m_initial_search_range(initial_search_range), m_kernel_size(kernel_size), m_slog_width(slog_width), m_cross_correlation_threshold(cross_correlation_threshold), m_pyramid_min_image_dimension(pyramid_min_image_dimension) {
      m_debug_prefix = "";
    }

    /// Turn on debugging output.  The debug_file_prefix string is
    /// used as a prefix for all debug image files.
    void set_debug_mode(std::string const& debug_file_prefix) { m_debug_prefix = debug_file_prefix; }

    
    template <class ViewT>
    ImageView<PixelDisparity<float> > operator() (ImageViewBase<ViewT> const& left_image,
                                                  ImageViewBase<ViewT> const& right_image) {
      typedef typename ViewT::pixel_type pixel_type;
      typedef typename PixelChannelType<pixel_type>::type channel_type;
      
      VW_ASSERT(left_image.impl().cols() == right_image.impl().cols() && left_image.impl().rows() == right_image.impl().rows(),
                ArgumentErr() << "Correlator(): input image dimensions do not match.");


      // Compute the number of pyramid levels needed to reach the
      // minimum image resolution set by the user.  
      // 
      // Level 0 : finest (high resolution) image
      // Level n : coarsest (low resolution) image ( n == pyramid_levels - 1 )
      int pyramid_levels = 1;
      if (m_pyramid_min_image_dimension != 0) {
        int dimension = std::min(left_image.impl().cols(), left_image.impl().rows());
        while (dimension > m_pyramid_min_image_dimension) { pyramid_levels++; dimension /= 2;}
      } 

      vw_out(InfoMessage) << "Initializing pyramid correlator with " << pyramid_levels << " levels.\n";

      // Build the image pyramid
      std::vector<ImageView<channel_type> > left_pyramid(pyramid_levels), right_pyramid(pyramid_levels);
      std::vector<ImageView<uint8> > left_slog_pyramid(pyramid_levels), right_slog_pyramid(pyramid_levels);
      std::vector<ImageView<bool> > left_masks(pyramid_levels), right_masks(pyramid_levels);

      left_pyramid[0] = channels_to_planes(left_image);
      right_pyramid[0] = channels_to_planes(right_image);
      left_masks[0] = disparity::generate_mask(left_image);
      right_masks[0] = disparity::generate_mask(right_image);
      left_slog_pyramid[0] = vw::channel_cast<uint8>(vw::threshold(vw::laplacian_filter(vw::gaussian_filter(channel_cast<float>(left_image),m_slog_width)), 0.0));
      right_slog_pyramid[0] = vw::channel_cast<uint8>(vw::threshold(vw::laplacian_filter(vw::gaussian_filter(channel_cast<float>(right_image),m_slog_width)), 0.0));
      
      // Apply the SLOG filter to each level.  
      // 
      // TODO: Refactor so that the user can provide other
      // pre-processing routines via a functor, such as LoG.
      for (int n = 1; n < pyramid_levels; ++n) {
        left_pyramid[n] = subsample_by_two(left_pyramid[n-1]);
        right_pyramid[n] = subsample_by_two(right_pyramid[n-1]);
        left_masks[n] = disparity::generate_mask(left_pyramid[n]);
        right_masks[n] = disparity::generate_mask(right_pyramid[n]);
        left_slog_pyramid[n] = vw::channel_cast<uint8>(vw::threshold(vw::laplacian_filter(vw::gaussian_filter(left_pyramid[n],m_slog_width)), 0.0));
        right_slog_pyramid[n] = vw::channel_cast<uint8>(vw::threshold(vw::laplacian_filter(vw::gaussian_filter(right_pyramid[n],m_slog_width)), 0.0));
      }

      // Free up some memory.
      left_pyramid.clear();
      right_pyramid.clear();

      return do_correlation(left_slog_pyramid, right_slog_pyramid, left_masks, right_masks);
    }
      
  };
  
}} // namespace vw::stereo

#endif // __VW_STEREO_CORRELATOR_H__
