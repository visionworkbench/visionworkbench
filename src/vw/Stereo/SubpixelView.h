// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifndef __VW_STEREO_AFFINE_SUBPIXEL_VIEW__
#define __VW_STEREO_AFFINE_SUBPIXEL_VIEW__

#include <vw/Image/ImageView.h>
#include <vw/Image/ImageViewRef.h>
#include <vw/Stereo/DisparityMap.h>
#include <vw/Stereo/Correlate.h>
#include <vw/Image/Filter.h>
#include <vw/Image/Interpolation.h>

namespace vw {
namespace stereo {

  /// An image view for performing image correlation
  template <class PreprocFilterT>
  class SubpixelView : public ImageViewBase<SubpixelView<PreprocFilterT> > {

    ImageViewRef<PixelDisparity<float> > m_disparity_map;
    ImageViewRef<float> m_left_image; 
    ImageViewRef<float> m_right_image;

    // General Settings
    int m_kern_width, m_kern_height;
    bool m_do_h_subpixel, m_do_v_subpixel;
    //bool m_do_affine_subpixel;
    int m_do_affine_subpixel;
    PreprocFilterT m_preproc_filter;
    bool m_verbose;

  public:
      typedef PixelDisparity<float> pixel_type;
      typedef pixel_type result_type;
      typedef ProceduralPixelAccessor<SubpixelView> pixel_accessor;
      
  
      //image subsampling by two
      ImageView<float> subsample_img_by_two(ImageView<float> &img) const {
      
      //determine the new size of the image
      printf("img: orig_w = %d, orig_h = %d\n", img.cols(), img.rows()); 
      int new_width;
      int new_height;
    
      new_width = img.cols()/2;
      new_height = img.rows()/2;

      printf("img: w = %d, h = %d, new_w = %d, new_h = %d\n", img.cols(), img.rows(), new_width, new_height); 
     
      ImageView<float> outImg(new_width, new_height, img.planes());		
      int32 i, j, p;
      int32 nw_i, new_j;
 
      ImageView<float> g_img;
      g_img = gaussian_filter(img, 1.5);
      printf("Gaussian blurring\n");

      for (p = 0; p < outImg.planes() ; p++) {
        for (i = 0; i < outImg.cols(); i++) {
          for (j = 0; j < outImg.rows(); j++) {
              
            outImg(i,j,p) = g_img(2*i, 2*j, p);
            /*
            outImg(i,j,p) = 0.0f;
            outImg(i,j,p) += img(2*i     , 2*j    ,p);
            outImg(i,j,p) += img(2*i + 1 , 2*j    ,p);
            outImg(i,j,p) += img(2*i     , 2*j + 1,p);
            outImg(i,j,p) += img(2*i + 1 , 2*j + 1,p);
            outImg(i,j,p) /= 4;
            */
          }
        }
      }
      
      return outImg;
    }

    //disparity map down-sampling by two
    ImageView<PixelDisparity<float> > subsample_disp_map_by_two(ImageView<PixelDisparity<float> >&disp) const {
      
      //determine the new size of the image 
      printf("disp: orig_w = %d, orig_h = %d\n", disp.cols(), disp.rows()); 
      int new_width;
      int new_height;
    
      new_width = disp.cols()/2;
      new_height = disp.rows()/2;
      printf("disp: new_w = %d, new_h = %d\n", new_width, new_height); 
      
      ImageView<PixelDisparity<float> > outDisp(new_width, new_height);		
      int32 i, j;
      int32 new_i, new_j;
      
      for (i = 0; i < new_width; i++) {
	for (j = 0; j < new_height; j++) {  
            
            new_i = i*2;
            new_j = j*2;
            
            outDisp(i,j).h() = 0.0f;
            outDisp(i,j).h() += disp(new_i     , new_j    ).h();
            outDisp(i,j).h() += disp(new_i + 1 , new_j    ).h();
            outDisp(i,j).h() += disp(new_i     , new_j + 1).h();
	    outDisp(i,j).h() += disp(new_i + 1 , new_j + 1).h();
      
            outDisp(i,j).h() /= 8;
            
            outDisp(i,j).v() = 0.0f;
            outDisp(i,j).v() += disp(new_i     , new_j    ).v();
            outDisp(i,j).v() += disp(new_i + 1 , new_j    ).v();
            outDisp(i,j).v() += disp(new_i     , new_j + 1).v();
	    outDisp(i,j).v() += disp(new_i + 1 , new_j + 1).v();

            outDisp(i,j).v() /= 8;

            outDisp(i,j)[2] = disp(new_i, new_j)[2];
        }
     }
    
     return outDisp;
    }

     //disparity map up-sampling by two
    ImageView<PixelDisparity<float> > upsample_disp_map_by_two(ImageView<PixelDisparity<float> > input_disp, int up_width, int up_height) const {

      ImageView<PixelDisparity<float> > outDisp(up_width, up_height);	
      ImageViewRef<PixelDisparity<float> > disp = edge_extend(input_disp, 
                                                              ConstantEdgeExtension());	

      for (unsigned j = 0; j < outDisp.rows(); ++j) {
        for (unsigned i = 0; i < outDisp.cols(); ++i) {

          float x = math::impl::_floor(float(i)/2.0), y = math::impl::_floor(float(j)/2.0);
          float normx = float(i)/2-x, normy = float(j)/2-y, norm1mx = 1.0-normx, norm1my = 1.0-normy;
          float im_00 = disp(x,y).h();
          float im_01 = disp(x,y+1).h();
          float im_10 = disp(x+1,y).h();
          float im_11 = disp(x+1,y+1).h();
          
          outDisp(i,j) = PixelDisparity<float>(2 * (disp(x,y).h()   * norm1mx*norm1my + 
                                                    disp(x+1,y).h() * normx*norm1my + 
                                                    disp(x,y+1).h() * norm1mx*normy + 
                                                    disp(x+1,y+1).h() * normx*normy),

                                               2 * (disp(x,y).v()   * norm1mx*norm1my + 
                                                    disp(x+1,y).v() * normx*norm1my + 
                                                    disp(x,y+1).v() * norm1mx*normy + 
                                                    disp(x+1,y+1).v() * normx*normy), 
                                               
                                               ceil(
                                                    disp(x,y).missing()   * norm1mx*norm1my + 
                                                    disp(x+1,y).missing() * normx*norm1my + 
                                                    disp(x,y+1).missing() * norm1mx*normy + 
                                                    disp(x+1,y+1).missing() * normx*normy
                                                    )
                                               );
        }
      }

      #if 0
      ImageView<PixelDisparity<float> > outDisp(up_width, up_height);			
      int32 i, j;
      int32 old_i, old_j;
      int valid;
      float new_disp_v, new_disp_h;

      for (i = 0; i < up_width; i++) {
	 
         old_i = i/2;
         if (old_i > disp.cols()-1){
	     old_i = disp.cols()-1;
         }   

         for (j = 0; j < up_height; j++) {
           
            old_j = j/2;
            if (old_j > disp.rows()-1){
	        old_j = disp.rows()-1;
            }          
            
            new_disp_h = 2*disp(old_i,old_j).h();
            new_disp_v = 2*disp(old_i,old_j).v();
            valid =  disp(old_i, old_j)[2];
	    
            outDisp(i,j).v() = new_disp_v;
            outDisp(i,j).h() = new_disp_h; 
            outDisp(i,j)[2] = valid;
            
            /*          
            if ( ((outDisp(i,j).h() == 0) || (outDisp(i,j).h() == 0)) && (valid!=0) ) { 
               printf("i=%d, j=%d, new_v=%f, new_h=%f, old_v = %f, old_h = %f, valid=%d\n",
		      i,j, outDisp(i,j).v(), outDisp(i,j).h(), disp(old_i, old_j).v(), disp(old_i, old_j).h(), valid);
            }
            */
            
         }
      }     
      #endif

      return outDisp;
    }

    template <class DisparityViewT, class InputViewT>
    SubpixelView(DisparityViewT const& disparity_map,
                 InputViewT const& left_image,
                 InputViewT const& right_image,
                 int kern_width, int kern_height,
                 bool do_horizontal_subpixel,
                 bool do_vertical_subpixel,
                 //bool do_affine_subpixel,
                 int do_affine_subpixel,
                 PreprocFilterT preproc_filter,
                 bool verbose) : m_disparity_map(disparity_map),
                                 m_left_image(left_image),
                                 m_right_image(right_image),
                                 m_kern_width(kern_width), m_kern_height(kern_height),
                                 m_do_h_subpixel(do_horizontal_subpixel),
                                 m_do_v_subpixel(do_vertical_subpixel),
                                 m_do_affine_subpixel(do_affine_subpixel),
                                 m_preproc_filter(preproc_filter),
                                 m_verbose(verbose) {
      // Basic assertions
      VW_ASSERT((left_image.impl().cols() == right_image.impl().cols()) &&
                (left_image.impl().rows() == right_image.impl().rows()) &&
                (disparity_map.impl().cols() == right_image.impl().cols()) &&
                (disparity_map.impl().cols() == right_image.impl().cols()),
                ArgumentErr() << "SubpixelView::SubpixelView(): input image dimensions and/or disparity_map dimensions do not agree.\n");
      
      VW_ASSERT((left_image.channels() == 1) && (left_image.impl().planes() == 1) &&
                (right_image.channels() == 1) && (right_image.impl().planes() == 1),
                ArgumentErr() << "SubpixelView::SubpixelView(): multi-channel, multi-plane images not supported.\n");
    }

    // Standard ImageView interface methods
    inline int32 cols() const { return m_left_image.cols(); }
    inline int32 rows() const { return m_left_image.rows(); }
    inline int32 planes() const { return 1; }
    
    inline pixel_accessor origin() const { return pixel_accessor( *this, 0, 0 ); }
    
    inline pixel_type operator()(float x, float y, int32 p = 0) const {
      vw_throw(NoImplErr() << "SubpixelView::operator() is not yet implemented.");
      return PixelDisparity<float>(); // Never reached
    }


    /// \cond INTERNAL
    typedef CropView<ImageView<pixel_type> > prerasterize_type;
    inline prerasterize_type prerasterize(BBox2i bbox) const { 

      // Find the range of disparity values for this patch.
      int num_good;
      BBox2i search_range = disparity::get_disparity_range(crop(m_disparity_map, bbox),num_good,false);

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
      right_crop_bbox.min() -= Vector2i(m_kern_width, m_kern_height);
      right_crop_bbox.max() += Vector2i(m_kern_width, m_kern_height);
      left_crop_bbox.min() -= Vector2i(m_kern_width, m_kern_height);
      left_crop_bbox.max() += Vector2i(m_kern_width, m_kern_height);

      // We crop the images to the expanded bounding box and edge
      // extend in case the new bbox extends past the image bounds.
      ImageView<float> left_image_patch, right_image_patch;
      if (m_do_affine_subpixel > 0) { // affine subpixel does its own pre-processing
        left_image_patch = crop(edge_extend(m_left_image,ZeroEdgeExtension()), left_crop_bbox);
        right_image_patch = crop(edge_extend(m_right_image,ZeroEdgeExtension()), right_crop_bbox);
      } else { // parabola subpixel does the same preprocessing as the pyramid correlator
        left_image_patch = m_preproc_filter(crop(edge_extend(m_left_image,ZeroEdgeExtension()), left_crop_bbox));
        right_image_patch = m_preproc_filter(crop(edge_extend(m_right_image,ZeroEdgeExtension()), right_crop_bbox));
      }
      ImageView<PixelDisparity<float> > disparity_map_patch = crop(edge_extend(m_disparity_map, ZeroEdgeExtension()), left_crop_bbox);

      // Adjust the disparities to be relative to the cropped
      // image pixel locations
      for (int v = 0; v < disparity_map_patch.rows(); ++v)
        for (int u = 0; u < disparity_map_patch.cols(); ++u)
          if (!disparity_map_patch(u,v).missing())  {
            disparity_map_patch(u,v).h() -= search_range.min().x();
            disparity_map_patch(u,v).v() -= search_range.min().y();
          }

      //       std::ostringstream ostr;
      //       ostr << "__" << bbox.min().x() << "_" << bbox.min().y() << ".tif";
      //       write_image("left"+ostr.str(), left_image_patch);
      //       write_image("right"+ostr.str(), right_image_patch);
      
     
    
      switch (m_do_affine_subpixel){

       
      case 0 : // Parabola Subpixel
        subpixel_correlation_parabola(disparity_map_patch,
                                      left_image_patch,
                                      right_image_patch,
                                      m_kern_width, m_kern_height,
                                      m_do_h_subpixel, m_do_v_subpixel,
                                      m_verbose);
        break;
      case 1: // Robust Subpixel
        subpixel_correlation_affine_2d(disparity_map_patch,
                                       left_image_patch,
                                       right_image_patch,
                                       m_kern_width, m_kern_height,
                                       m_do_h_subpixel, m_do_v_subpixel,
                                       m_verbose);
        break;
      case 2: // Bayes Subpixel
        subpixel_correlation_affine_2d_bayesian(disparity_map_patch,
                                                left_image_patch,
                                                right_image_patch,
                                                m_kern_width, m_kern_height,
                                                m_do_h_subpixel, m_do_v_subpixel,
                                                m_verbose);
        break;
      case 3: // Bayes EM  Subpixel
        {
	  
        int pyramid_levels = 2;
        printf("test0\n");
        
        //create the pyramid first
        std::vector<ImageView<float> > left_pyramid(pyramid_levels), right_pyramid(pyramid_levels);
        left_pyramid[0] = channels_to_planes(left_image_patch);
        right_pyramid[0] = channels_to_planes(right_image_patch);
	printf("test1\n");

        std::vector<ImageView<PixelDisparity<float> > > disparity_map_pyramid(pyramid_levels);
        std::vector<ImageView<PixelDisparity<float> > > disparity_map_upsampled(pyramid_levels);
        disparity_map_pyramid[0] = disparity_map_patch;

        printf("test2\n");
        //downsample the disparity map and the image pair
        for (int i = 1; i < pyramid_levels; i++) {
          left_pyramid[i] = subsample_img_by_two(left_pyramid[i-1]);
          right_pyramid[i] = subsample_img_by_two(right_pyramid[i-1]);
          disparity_map_pyramid[i] = subsample_disp_map_by_two(disparity_map_pyramid[i-1]);
        }
        
        subpixel_correlation_affine_2d_EM(disparity_map_pyramid[pyramid_levels-1],
                                          left_pyramid[pyramid_levels-1],
                                          right_pyramid[pyramid_levels-1],
                                          m_kern_width, m_kern_height,
                                          m_do_h_subpixel, m_do_v_subpixel,
                                          m_verbose);	
        disparity_map_upsampled[pyramid_levels-1] = disparity_map_pyramid[pyramid_levels-1];

        for (int i = pyramid_levels-2; i>=0; i--){

          // For Debugging
          // std::ostringstream ostr2;
          // ostr2 << "subsamp-" << i << ".exr";
          // write_image(ostr2.str(), disparity_map_upsampled[i+1]);

          int up_width = left_pyramid[i].cols();
          int up_height = left_pyramid[i].rows();
          std::cout << "\n\n--> Upsampling to " << up_width << " " << up_height << "\n";
          disparity_map_upsampled[i] = upsample_disp_map_by_two(disparity_map_upsampled[i+1], 
                                                                up_width, up_height);
          
          // For Debugging
          // std::ostringstream ostr;
          // ostr << "upsamp-" << i << ".exr";
          // write_image(ostr.str(), disparity_map_upsampled[i]);

          printf("disp_map_width = %d\n", disparity_map_upsampled[i].cols());
          printf("disp_map_height = %d\n", disparity_map_upsampled[i].rows());
          printf("left_pyramid_width = %d\n", left_pyramid[i].cols());
          printf("left_pyramid_height = %d\n", left_pyramid[i].rows());
          
          subpixel_correlation_affine_2d_EM(disparity_map_upsampled[i],
                                            left_pyramid[i],
                                            right_pyramid[i],
                                            m_kern_width, m_kern_height,
                                            m_do_h_subpixel, m_do_v_subpixel,
                                            m_verbose);
	}

        printf("disp_map_width = %d, height = %d\n", disparity_map_patch.cols(), disparity_map_patch.rows());
        printf("disp_map_up_width = %d, height = %d\n", disparity_map_upsampled[0].cols(), disparity_map_upsampled[0].rows());
	
        //disparity_map_patch = disparity_map_upsampled[0];
        for (int i = 0; i < disparity_map_patch.cols(); i++) {
          for (int j = 0; j < disparity_map_patch.rows(); j++) {  
            
	    disparity_map_patch(i,j).h() = disparity_map_upsampled[0](i,j).h();
	    disparity_map_patch(i,j).v() = disparity_map_upsampled[0](i,j).v();
	    disparity_map_patch(i,j)[2]  = disparity_map_upsampled[0](i,j)[2];
          
          }
       }
       
	/*
        subpixel_correlation_affine_2d_EM(disparity_map_patch,
                                          left_image_patch,
                                          right_image_patch,
                                          m_kern_width, m_kern_height,
                                          m_do_h_subpixel, m_do_v_subpixel,
                                          m_verbose);
	*/
	  
	}
        break;
      default:
        vw_throw(ArgumentErr() << "Unknown subpixel correlation type: " << m_do_affine_subpixel << ".");
        break;
      }
      
      // Undo the above adjustment
      for (int v = 0; v < disparity_map_patch.rows(); ++v)
        for (int u = 0; u < disparity_map_patch.cols(); ++u)
          if (!disparity_map_patch(u,v).missing())  {
            disparity_map_patch(u,v).h() += search_range.min().x();
            disparity_map_patch(u,v).v() += search_range.min().y();
          }

      // This may seem confusing, but we must crop here so that the
      // good pixel data is placed into the coordinates specified by
      // the bbox.  This allows rasterize to touch those pixels
      // using the coordinates inside the bbox.  The pixels outside
      // those coordinates are invalid, and they never get accessed.
      return crop(disparity_map_patch, BBox2i(m_kern_width-bbox.min().x(), 
                                              m_kern_height-bbox.min().y(), 
                                              m_left_image.cols(), 
                                              m_left_image.rows() ));
    }
    
    template <class DestT> inline void rasterize(DestT const& dest, BBox2i bbox) const {
      vw::rasterize(prerasterize(bbox), dest, bbox);
    }
    /// \endcond

  };

}} // namespace vw::stereo

#endif // __VW_STEREO_CORRELATOR_VIEW__         
