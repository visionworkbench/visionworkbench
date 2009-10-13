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

    static ImageView<float> subsample_img_by_two(ImageView<float> &img);
    static ImageView<PixelMask<Vector2f> >
      subsample_disp_map_by_two(ImageView<PixelMask<Vector2f> > const& input_disp);
    static ImageView<PixelMask<Vector2f> >
      upsample_disp_map_by_two(ImageView<PixelMask<Vector2f> > const& input_disp,
                               int up_width, int up_height);

    template <class ImagePixelT>
      class EMSubpixelCorrelatorView : public ImageViewBase<EMSubpixelCorrelatorView<ImagePixelT> > {
    public:
      typedef PixelMask<Vector2f> pixel_type;
      typedef pixel_type result_type;
      typedef ProceduralPixelAccessor<EMSubpixelCorrelatorView> pixel_accessor;

      template <class ImageT, class DisparityT>
        EMSubpixelCorrelatorView(ImageViewBase<ImageT> const& left_image, ImageViewBase<ImageT> const& right_image,
                                 ImageViewBase<DisparityT> const& course_disparity, bool debug = false);

      // EM parameter setters
      void set_kernel_size(Vector2i size) { m_kernel_size = size; }
      void set_em_iter_max(int iter) { em_iter_max = iter; }
      void set_em_epsilon(double epsilon) { epsilon_em = epsilon; }
      void set_prob_inlier_0(double P) { P_inlier_0 = P; }
      void set_prob_inlier_min(double P) { P_inlier_min = P; }
      void set_prob_inlier_max(double P) { P_inlier_max = P; }
      void set_sigma_affine_0(double s) { sigma_p1_0 = s; }
      void set_sigma_affine_min(double s) { sigma_p1_min = s; }
      // affine model parameter setters
      void set_inner_iter_max(double iter) { inner_iter_max = iter; }
      void set_inner_epsilon(double epsilon) { epsilon_inner = epsilon; }
      void set_min_determinant(double min) { affine_min_det = min; }
      void set_max_determinant(double max) { affine_max_det = max; }


      Vector2i kernel_size() const { return m_kernel_size; }

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
      ImageViewRef<pixel_type> m_course_disparity;

      // EM settings
      Vector2i m_kernel_size;
      int em_iter_max;
      double epsilon_em;
      double P_inlier_0;
      double P_inlier_min;
      double P_inlier_max;
      double sigma_p1_0;
      double sigma_p1_min;
      double sigma_n_0;
      double sigma_n_min;
      double mu_n_0;
      // affine model settings
      int inner_iter_max;
      double epsilon_inner;
      double affine_min_det;
      double affine_max_det;

      bool m_debug;

      // private helper methods
      template <class ImageT, class DisparityT, class AffineT>
        inline void
        m_subpixel_refine(ImageViewBase<ImageT> const& left_image, ImageViewBase<ImageT> const& right_image,
                          ImageViewBase<DisparityT> &disparity, ImageViewBase<AffineT> & affine_warps,
                          BBox2i const& ROI, bool final, bool debug = false) const;
    }; // end class


    /// Summarize a CorrelatorView object
    template <class ImagePixelT>
      inline std::ostream& operator<<( std::ostream& os, EMSubpixelCorrelatorView<ImagePixelT> const& view );

  }
}

#include "EMSubpixelCorrelatorView.hpp"

#endif
