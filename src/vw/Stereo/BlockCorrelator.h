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
    
    int m_lKernWidth, m_lKernHeight;
    int m_lMinH, m_lMaxH, m_lMinV, m_lMaxV;
    int m_block_size;
    int m_verbose;
    double m_crossCorrThreshold;
    int m_useHorizSubpixel;
    int m_useVertSubpixel;
    
  public:
    BlockCorrelator(int minH,	/* left bound disparity search window*/
                    int maxH,	/* right bound disparity search window*/
                    int minV,	/* bottom bound disparity search window */ 
                    int maxV,	/* top bound disparity search window */
                    int kernWidth,	/* size of the kernel */
                    int kernHeight,       
                    int verbose,
                    double crosscorrThreshold,
                    int block_size,
                    int useSubpixelH, int useSubpixelV) {

      m_lKernWidth = kernWidth;
      m_lKernHeight = kernHeight;
      m_lMinH = minH;
      m_lMaxH = maxH;
      m_lMinV = minV;
      m_lMaxV = maxV;  
      m_verbose = verbose;
      m_block_size = block_size;
      
      m_crossCorrThreshold = crosscorrThreshold;
      m_useHorizSubpixel = useSubpixelH;
      m_useVertSubpixel = useSubpixelV;
    }

    // Constrain a set of nominal image blocks to the specified
    // disparity search range.
    void constrain_to_search_range(std::vector<BBox2i> &nominal_blocks, int width, int height) {      
      int min_x = std::max(0, -m_lMinH);
      int max_x = std::min(width, width-m_lMaxH);
      int min_y = std::max(0, -m_lMinV);
      int max_y = std::min(height, height-m_lMaxV);
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
        right_blocks[i] = BBox2i(Vector2i(nominal_blocks[i].min().x()+m_lMinH,
                                          nominal_blocks[i].min().y()+m_lMinV),
                                 Vector2i(nominal_blocks[i].max().x()+m_lMaxH,
                                          nominal_blocks[i].max().y()+m_lMaxV));
        left_blocks[i].max() = Vector2i(nominal_blocks[i].min().x() + right_blocks[i].width(),
                                        nominal_blocks[i].min().y() + right_blocks[i].height());
        right_blocks[i].min() -= Vector2i(m_lKernWidth, m_lKernHeight);
        right_blocks[i].max() += Vector2i(m_lKernWidth, m_lKernHeight);
        left_blocks[i].min() -= Vector2i(m_lKernWidth, m_lKernHeight);
        left_blocks[i].max() += Vector2i(m_lKernWidth, m_lKernHeight);
      }
    }

    template <class PixelT>
    ImageView<PixelDisparity<float> > operator()(vw::ImageView<PixelT>& left_image, 
                                                 vw::ImageView<PixelT>& right_image,
                                                 bool bit_image) {
      
      VW_ASSERT(left_image.cols() == right_image.cols() &&
                left_image.rows() == right_image.rows(),
                ArgumentErr() << "block_correlator: input image dimensions do not agree.\n");
      
      VW_ASSERT(left_image.channels() == 1 && left_image.planes()==1 &&
                right_image.channels() == 1 && right_image.planes()==1,
                ArgumentErr() << "block_correlator: does not support multi-channel, multi-plane images.\n");
      
      //Run the correlator and record how long it takes to run.
      double begin__ = Time();
      
      // Ask the worker threads for the actual results of the disparity correlation
      ImageView<typename PixelChannelType<PixelT>::type> l_image = channels_to_planes(left_image);
      ImageView<typename PixelChannelType<PixelT>::type> r_image = channels_to_planes(right_image);        

      // Allocate memory for the resulting disparity image.
      ImageView<PixelDisparity<float> > disparity_map(left_image.cols(), left_image.rows());

      // Divide the left image into blocks and compute the bounding
      // boxes in the left and right image given the disparity range
      // that this object was initialized with.
      std::vector<BBox2i> nominal_left_blocks = image_blocks(l_image, m_block_size, m_block_size);
      constrain_to_search_range(nominal_left_blocks, l_image.cols(), l_image.rows());
      std::vector<BBox2i> left_blocks, right_blocks;
      compute_blocks(nominal_left_blocks, left_blocks, right_blocks);

      // Create an optimized correlator to run on each block.
      vw::stereo::OptimizedCorrelator correlator( 0, m_lMaxH-m_lMinH, 
                                                  0, m_lMaxV-m_lMinV,
                                                  m_lKernWidth, m_lKernHeight,
                                                  true, m_crossCorrThreshold,
                                                  m_useHorizSubpixel,
                                                  m_useVertSubpixel );

      // Run the image blocks through the correlator.  Fill the disparity map as you go.
      std::cout << "Processing Blocks: \n";
      for (int i = 0; i < left_blocks.size();++i) {
        std::cout << i+1 << " / " << left_blocks.size() << ".        \r";
        ImageView<PixelDisparity<float> > disparity_subregion;
        ImageView<PixelT> left_subimage = crop(edge_extend(l_image, ReflectEdgeExtend()), left_blocks[i]);
        ImageView<PixelT> right_subimage = crop(edge_extend(r_image, ReflectEdgeExtend()), right_blocks[i]);
        
        disparity_subregion = correlator( left_subimage, right_subimage, bit_image );
        
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
              disparity_subregion(u,v).h() += m_lMinH;
              disparity_subregion(u,v).v() += m_lMinV;
            }
          }
        }
        crop(disparity_map, nominal_left_blocks[i]) = crop(disparity_subregion, m_lKernWidth,m_lKernHeight, 
                                                           nominal_left_blocks[i].width(),
                                                           nominal_left_blocks[i].height());
      }
      //      std::cout << "done.                 \n";
     
      double lapse__ = Time() - begin__;
      std::cout << "\tTotal block correlation took " << lapse__ << " seconds.\n";
      
      return disparity_map;
  }
};



}} // namespace vw::stereo

#endif // __VW_STEREO_BLOCK_CORRELATOR__
