#ifndef __VW_STEREO_CORRELATOR_VIEW__
#define __VW_STEREO_CORRELATOR_VIEW__
#include <vw/Stereo/Correlate.h>
#include <vw/Stereo/Correlator.h>
#include <vw/Stereo/DisparityMap.h>
#include <vw/Image/ImageViewRef.h>
#include <vw/Image/Manipulation.h>
#include <vw/Math/BBox.h>
#include <vw/Core/Debugging.h>

namespace vw {
  namespace stereo {

    struct CorrelationSettings
    {
      CorrelationSettings() {}
      CorrelationSettings(int minH,	// left bound disparity search window
			  int maxH,	// right bound disparity search window
			  int minV,	// bottom bound disparity search window
			  int maxV,	// top bound disparity search window
			  int kernWidth,	// size of the kernel
			  int kernHeight,       
			  int verbose,
			  double crosscorrThreshold,
        float corrscore_rejection_threshold,
        double slog_width,
			  int useSubpixelH, int useSubpixelV,
			  bool bit_image)
      {
        m_lKernWidth = kernWidth;
        m_lKernHeight = kernHeight;
        m_lMinH = minH;
        m_lMaxH = maxH;
        m_lMinV = minV;
        m_lMaxV = maxV;  
        m_verbose = verbose;
        m_crossCorrThreshold = crosscorrThreshold;
        m_corrscore_rejection_threshold = corrscore_rejection_threshold;
        m_slog_width = slog_width;
        m_useHorizSubpixel = useSubpixelH;
        m_useVertSubpixel = useSubpixelV;
        m_bit_image = bit_image;
      }
      int m_lKernWidth, m_lKernHeight;
      int m_lMinH, m_lMaxH, m_lMinV, m_lMaxV;
      int m_verbose;
      double m_crossCorrThreshold;
      float m_corrscore_rejection_threshold;
      double m_slog_width;
      int m_useHorizSubpixel;
      int m_useVertSubpixel;
      bool m_bit_image;
    };

    /// An image view for performing image correlation
    template <class PixelT>
    class CorrelatorView : public ImageViewBase<CorrelatorView<PixelT> >
    {    
    public:
      typedef PixelDisparity<float> pixel_type;
      typedef pixel_type result_type;
      typedef ProceduralPixelAccessor<CorrelatorView> pixel_accessor;
      
      template <class ImageT>
      CorrelatorView(ImageViewBase<ImageT> const& left_image, ImageViewBase<ImageT> const& right_image, const CorrelationSettings &settings) :
        m_left_image(left_image.impl()), m_right_image(right_image.impl()), m_settings(settings)
      {
        VW_ASSERT((left_image.impl().cols() == right_image.impl().cols()) &&
                  (left_image.impl().rows() == right_image.impl().rows()),
                  ArgumentErr() << "CorrelatorView::CorrelatorView(): input image dimensions do not agree.\n");
        VW_ASSERT((left_image.channels() == 1) && (left_image.impl().planes() == 1) &&
                  (right_image.channels() == 1) && (right_image.impl().planes() == 1),
                  ArgumentErr() << "CorrelatorView::CorrelatorView(): multi-channel, multi-plane images not supported.\n");
      }

      inline int32 cols() const { return m_left_image.cols(); }
      inline int32 rows() const { return m_left_image.rows(); }
      inline int32 planes() const { return 1; }
      
      inline pixel_accessor origin() const { return pixel_accessor( *this, 0, 0 ); }

      inline pixel_type operator()(double i, double j, int32 p = 0) const
      {
        vw_throw(NoImplErr() << "CorrelatorView::operator()(double i, double j, int32 p) has not been implemented.");
        return pixel_type();
      }

      /// \cond INTERNAL
      typedef CropView<ImageView<pixel_type> > prerasterize_type;
      inline prerasterize_type prerasterize(BBox2i bbox) const
      {
        vw_out(InfoMessage) << "\n\tBlock: " << bbox << "\n";
        
        BBox2i search_bbox(Vector2i(m_settings.m_lMinH, m_settings.m_lMinV),
                           Vector2i(m_settings.m_lMaxH, m_settings.m_lMaxV));

        // The area in the right image that we'll be searching is
        // determined by the bbox of the left image plus the search
        // range.
        BBox2i left_crop_bbox(bbox);
        BBox2i right_crop_bbox(bbox.min() + search_bbox.min(),
                               bbox.max() + search_bbox.max());
        
        // The correlator requires the images to be the same size. The
        // search bbox will always be larger than the given left image
        // bbox, so we just make the left bbox the same size as the
        // right bbox.
        left_crop_bbox.max() = left_crop_bbox.min() +
          Vector2i(right_crop_bbox.width(), right_crop_bbox.height());

        // Finally, we must adjust both bounding boxes to account for
        // the size of the kernel itself.
        right_crop_bbox.min() -= Vector2i(m_settings.m_lKernWidth, m_settings.m_lKernHeight);
        right_crop_bbox.max() += Vector2i(m_settings.m_lKernWidth, m_settings.m_lKernHeight);
        left_crop_bbox.min() -= Vector2i(m_settings.m_lKernWidth, m_settings.m_lKernHeight);
        left_crop_bbox.max() += Vector2i(m_settings.m_lKernWidth, m_settings.m_lKernHeight);

//         std::cout << "\nCorrelatorView::prerasterize(): search_bbox: " << search_bbox << std::endl;
//         std::cout << "\n                             left_crop_bbox: " << left_crop_bbox << std::endl;
//         std::cout << "\n                            right_crop_bbox: " << right_crop_bbox << std::endl;

        // We crop the images to the expanded bounding box and edge
        // extend in case the new bbox extends past the image bounds.
        ImageView<PixelT> cropped_left_image = crop(edge_extend(m_left_image, ReflectEdgeExtension()), left_crop_bbox);
        ImageView<PixelT> cropped_right_image = crop(edge_extend(m_right_image, ReflectEdgeExtension()), right_crop_bbox);

        // We have all of the settings adjusted.  Now we just have to
        // run the correlator.
        ImageView<pixel_type> disparity_map;

        vw::stereo::Correlator correlator(BBox2(0,0,search_bbox.width(),search_bbox.height()),
                                          Vector2i(m_settings.m_lKernWidth, m_settings.m_lKernHeight),
                                          m_settings.m_slog_width,
                                          m_settings.m_crossCorrThreshold,
                                          m_settings.m_corrscore_rejection_threshold,
                                          m_settings.m_useHorizSubpixel,
                                          m_settings.m_useVertSubpixel);
        // For debugging: this saves the disparity map at various pyramid levels to disk.
        //        correlator.set_debug_mode("debug");
        disparity_map = correlator(cropped_left_image, cropped_right_image);

        // Adjust the disparities to be relative to the uncropped
        // image pixel locations
        for (int v = 0; v < disparity_map.rows(); ++v)
          for (int u = 0; u < disparity_map.cols(); ++u)
            if (!disparity_map(u,v).missing())  {
              disparity_map(u,v).h() += m_settings.m_lMinH;
              disparity_map(u,v).v() += m_settings.m_lMinV;
            }

        // This may seem confusing, but we must crop here so that the
        // good pixel data is placed into the coordinates specified by
        // the bbox.  This allows rasterize to touch those pixels
        // using the coordinates inside the bbox.  The pixels outside
        // those coordinates are invalid, but they never get accessed.
        return CropView<ImageView<pixel_type> > (disparity_map, BBox2i(m_settings.m_lKernWidth-bbox.min().x(), 
                                                                       m_settings.m_lKernHeight-bbox.min().y(), 
                                                                       bbox.width(), bbox.height()));
      }

      template <class DestT> inline void rasterize(DestT const& dest, BBox2i bbox) const {
        vw::rasterize(prerasterize(bbox), dest, bbox);
      }
      /// \endcond

    private:
      ImageViewRef<PixelT> m_left_image, m_right_image;
      CorrelationSettings m_settings;
    };
  }
} // namespace vw::stereo

#endif // __VW_STEREO_CORRELATOR_VIEW__         
