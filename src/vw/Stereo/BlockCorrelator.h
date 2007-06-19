#ifndef __VW_STEREO_BLOCK_CORRELATOR__
#define __VW_STEREO_BLOCK_CORRELATOR__

#include <vw/Stereo/Correlate.h>
#include <vw/Stereo/OptimizedCorrelator.h>
#include <vw/Stereo/DisparityMap.h>
#include <vw/Image/Manipulation.h>
#include <vw/Math/BBox.h>
#include <vw/Core/Debugging.h>

namespace vw {
namespace stereo {

  class BlockCorrelator {
    
    int m_kernel_width, m_kernel_height;
    int m_min_h, m_max_h, m_min_v, m_max_v;
    int m_block_size;
    int m_verbose;
    double m_cross_corr_threshold;
    float m_corrscore_rejection_threshold;
    int m_use_horiz_subpixel;
    int m_use_vert_subpixel;
    
  public:
    BlockCorrelator(int minH,	/* left bound disparity search window*/
                    int maxH,	/* right bound disparity search window*/
                    int minV,	/* bottom bound disparity search window */ 
                    int maxV,	/* top bound disparity search window */
                    int kernWidth,	/* size of the kernel */
                    int kernHeight,       
                    int verbose,
                    double cross_corr_threshold,
                    float corrscore_rejection_threshold,
                    int block_size,
                    int useSubpixelH, int useSubpixelV) {

      m_kernel_width = kernWidth;
      m_kernel_height = kernHeight;
      m_min_h = minH;
      m_max_h = maxH;
      m_min_v = minV;
      m_max_v = maxV;  
      m_verbose = verbose;
      m_block_size = block_size;
      
      m_cross_corr_threshold = cross_corr_threshold;
      m_corrscore_rejection_threshold = corrscore_rejection_threshold;
      m_use_horiz_subpixel = useSubpixelH;
      m_use_vert_subpixel = useSubpixelV;
    }

    // Constrain a set of nominal image blocks to the specified
    // disparity search range.
    void constrain_to_search_range(std::vector<BBox2i> &nominal_blocks, int width, int height) {      
      int min_x = std::max(0, -m_min_h);
      int max_x = std::min(width, width-m_max_h);
      int min_y = std::max(0, -m_min_v);
      int max_y = std::min(height, height-m_max_v);

      BBox2i workspace( Vector2i(min_x, min_y), Vector2i(max_x, max_y) );

      // Check for bounding boxes that con be removed entirely
      // because there is no overlap with the workspace.
      std::vector<BBox2i>::iterator iter = nominal_blocks.begin();
      while ( iter != nominal_blocks.end() ) {
        if ( !workspace.intersects(*iter) ) {
          nominal_blocks.erase(iter);
          iter = nominal_blocks.begin();
        } else {
          ++iter;
        }
      }

      // Check for bounding boxes with partial overlap.  Here we
      // just resize the bounding boxes so that the area outside the
      // workspace is eliminated.
      iter = nominal_blocks.begin();
      while ( iter != nominal_blocks.end() ) {
        iter->min().x() = (iter->min().x() < min_x) ? min_x : iter->min().x();
        iter->max().x() = (iter->max().x() > max_x) ? max_x : iter->max().x();
        iter->min().y() = (iter->min().y() < min_y) ? min_y : iter->min().y();
        iter->max().y() = (iter->max().y() > max_y) ? max_y : iter->max().y();
        ++iter;
      }
    }

    // Iterate over the nominal blocks, creating output blocks for correlation
    //
    // To compute the block size, we must 1) Adjust the and size of
    // the blocks for the right image to account for the disparity
    // search range, 2) adjust the size of the left image blocks to
    // be equal to the size of the right image blocks, and 3) buffer
    // both blocks by the kernel size.
    void compute_blocks(std::vector<BBox2i> const& nominal_blocks, 
                        std::vector<BBox2i> &left_blocks, std::vector<BBox2i> &right_blocks) {
      
      left_blocks.resize(nominal_blocks.size());
      right_blocks.resize(nominal_blocks.size());

      for (int i = 0; i < nominal_blocks.size(); ++i) {
        left_blocks[i] = nominal_blocks[i];
        right_blocks[i] = BBox2i(Vector2i(nominal_blocks[i].min().x()+m_min_h,
                                          nominal_blocks[i].min().y()+m_min_v),
                                 Vector2i(nominal_blocks[i].max().x()+m_max_h,
                                          nominal_blocks[i].max().y()+m_max_v));
        left_blocks[i].max() = Vector2i(nominal_blocks[i].min().x() + right_blocks[i].width(),
                                        nominal_blocks[i].min().y() + right_blocks[i].height());
        right_blocks[i].min() -= Vector2i(m_kernel_width, m_kernel_height);
        right_blocks[i].max() += Vector2i(m_kernel_width, m_kernel_height);
        left_blocks[i].min() -= Vector2i(m_kernel_width, m_kernel_height);
        left_blocks[i].max() += Vector2i(m_kernel_width, m_kernel_height);
      }
    }

    template <class ViewT>
    ImageView<PixelDisparity<float> > operator()(vw::ImageViewBase<ViewT> const& left_image, 
                                                 vw::ImageViewBase<ViewT> const& right_image,
                                                 bool bit_image) {
      
      VW_ASSERT(left_image.impl().cols() == right_image.impl().cols() &&
                left_image.impl().rows() == right_image.impl().rows(),
                ArgumentErr() << "block_correlator: input image dimensions do not agree.\n");
      
      VW_ASSERT(left_image.impl().channels() == 1 && left_image.impl().planes()==1 &&
                right_image.impl().channels() == 1 && right_image.impl().planes()==1,
                ArgumentErr() << "block_correlator: does not support multi-channel, multi-plane images.\n");
      
      //Run the correlator and record how long it takes to run.
      double begin__ = Time();
      
      // Allocate memory for the resulting disparity image.
      ImageView<PixelDisparity<float> > disparity_map(left_image.impl().cols(), left_image.impl().rows());

      // Divide the left image into blocks and compute the bounding
      // boxes in the left and right image given the disparity range
      // that this object was initialized with.
      std::vector<BBox2i> nominal_left_blocks = image_blocks(left_image.impl(), m_block_size, m_block_size);
      std::vector<BBox2i> left_blocks, right_blocks;
      compute_blocks(nominal_left_blocks, left_blocks, right_blocks);

      // If we only plan to process a single block, we do it here.
      if (left_blocks.size() == 1) {
        ImageView<typename PixelChannelType<typename ViewT::pixel_type>::type> l_image = channels_to_planes(left_image);
        ImageView<typename PixelChannelType<typename ViewT::pixel_type>::type> r_image = channels_to_planes(right_image);        

        // Create an optimized correlator to run on a single block.
        vw::stereo::OptimizedCorrelator correlator( m_min_h, m_max_h, 
                                                    m_min_v, m_max_v,
                                                    m_kernel_width, m_kernel_height,
                                                    true, m_cross_corr_threshold,
                                                    m_corrscore_rejection_threshold,
                                                    m_use_horiz_subpixel,
                                                    m_use_vert_subpixel );
        
        disparity_map = correlator( l_image, r_image );

      // Otherwise, we process the blocks individually
      } else {

        // Create an optimized correlator to run on each block.
        vw::stereo::OptimizedCorrelator correlator( 0, m_max_h-m_min_h, 
                                                    0, m_max_v-m_min_v,
                                                    m_kernel_width, m_kernel_height,
                                                    true, m_cross_corr_threshold,
                                                    m_corrscore_rejection_threshold,
                                                    m_use_horiz_subpixel,
                                                    m_use_vert_subpixel );


        // Run the image blocks through the correlator.  Fill the disparity map as you go.
        std::cout << "Processing Blocks: \n";
        for (int i = 0; i < left_blocks.size();++i) {
          std::cout << i+1 << " / " << left_blocks.size() << ".        \r";
          ImageView<PixelDisparity<float> > disparity_subregion;
          ImageView<typename ViewT::pixel_type> left_subimage = crop(edge_extend(channels_to_planes(left_image), ReflectEdgeExtension()), left_blocks[i]);
          ImageView<typename ViewT::pixel_type> right_subimage = crop(edge_extend(channels_to_planes(right_image), ReflectEdgeExtension()), right_blocks[i]);
          
          disparity_subregion = correlator( left_subimage, right_subimage );
          
          // Debugging:
          //         std::ostringstream ostream;
          //         ostream << i;
          //         double min_h_disp, min_v_disp, max_h_disp, max_v_disp;
          //         disparity::get_disparity_range(disparity_subregion, min_h_disp, max_h_disp, min_v_disp, max_v_disp,true);
          //         write_image( "disparityblock-" + ostream.str() + "-H.jpg", normalize(clamp(select_channel(disparity_subregion,0), min_h_disp, max_h_disp)));
          //         write_image( "disparityblock-" + ostream.str() + "-V.jpg", normalize(clamp(select_channel(disparity_subregion,1), min_v_disp, max_v_disp)));
          
          // Adjust the disparity range in the sub-region and place it in the overall disparity map.
          for (int v = 0; v < disparity_subregion.rows(); ++v) {
            for (int u = 0; u < disparity_subregion.cols(); ++u) {
              if (!disparity_subregion(u,v).missing()) {
                disparity_subregion(u,v).h() += m_min_h;
                disparity_subregion(u,v).v() += m_min_v;
              }
            }
          }
          crop(disparity_map, nominal_left_blocks[i]) = crop(disparity_subregion, 
                                                             m_kernel_width,m_kernel_height, 
                                                             nominal_left_blocks[i].width(),
                                                             nominal_left_blocks[i].height());
        }
      }
     
      double lapse__ = Time() - begin__;
      std::cout << "\tTotal block correlation took " << lapse__ << " seconds.\n";
      
      return disparity_map;
  }
};



}} // namespace vw::stereo

#endif // __VW_STEREO_BLOCK_CORRELATOR__
