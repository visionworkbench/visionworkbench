// __BEGIN_LICENSE__
//  Copyright (c) 2006-2013, United States Government as represented by the
//  Administrator of the National Aeronautics and Space Administration. All
//  rights reserved.
//
//  The NASA Vision Workbench is licensed under the Apache License,
//  Version 2.0 (the "License"); you may not use this file except in
//  compliance with the License. You may obtain a copy of the License at
//  http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
// __END_LICENSE__


#ifndef __VW_STEREO_EM_SUBPIXEL_CORRELATOR_VIEW__
#define __VW_STEREO_EM_SUBPIXEL_CORRELATOR_VIEW__

#include <vw/Math/Matrix.h>
#include <vw/Math/LinearAlgebra.h>
#include <vw/Image/Filter.h>
#include <vw/Image/ImageView.h>
#include <vw/Image/ImageViewBase.h>
#include <vw/Image/ImageViewRef.h>
#include <vw/Image/ImageMath.h>
#include <vw/Image/Manipulation.h>
#include <vw/Image/PixelMask.h>
#include <vw/Image/Transform.h>
#include <vw/Stereo/AffineMixtureComponent.h>
#include <vw/Stereo/GaussianMixtureComponent.h>

#include <ostream>

// For the PixelDisparity math.
#include <boost/operators.hpp>

namespace vw {
  namespace stereo {

    template <class ImagePixelT>
      class EMSubpixelCorrelatorView : public ImageViewBase<EMSubpixelCorrelatorView<ImagePixelT> > {
    public:
      typedef PixelMask<Vector2f> disparity_pixel;
      typedef PixelMask<Vector<float, 5> > pixel_type;
      typedef pixel_type result_type;
      typedef ProceduralPixelAccessor<EMSubpixelCorrelatorView> pixel_accessor;

      template <class Image1T, class Image2T, class DisparityT>
      EMSubpixelCorrelatorView(ImageViewBase<Image1T> const& left_image, ImageViewBase<Image2T> const& right_image,
                               ImageViewBase<DisparityT> const& course_disparity, int debug = -1);


      struct TfmCovarianceDisparityTo3DFunctor : public vw::ReturnFixedType<float>  {
        float operator()(Vector<float, 3> const & u_in, Matrix<float, 3, 2> const & j_in) const {
          Matrix2x2 temp_1;
          temp_1(0,0) = u_in(0);
          temp_1(0,1) = u_in(1);
          temp_1(1,0) = u_in(1);
          temp_1(1,1) = u_in(2);
          Matrix3x3 temp_2 = j_in*temp_1*transpose(j_in);
          Vector3f s;
          svd(temp_2, s);
          //std::cout << s(0) << " " << s(1) << " " << s(2) << std::endl;
          return pow(s(0),.25); // svd's are eigen-values squares for symmetric matrices; std deviation is sqrt of eigen values
        }
      };


      struct ExtractDisparityFunctor : public vw::ReturnFixedType<PixelMask<Vector2f> >  {
        PixelMask<Vector2f> operator()(PixelMask<Vector<float, 5> > const & in) const {
          PixelMask<Vector2f> temp(subvector(in.child(), 0, 2));
          if(!in.valid()) {
            temp.invalidate();
          }
          return temp;
        }
      };

      struct ExtractUncertaintyFunctor : public vw::ReturnFixedType<Vector3f>  {
        Vector3f operator()(PixelMask<Vector<float, 5> > const & in) const {
          Vector3f temp(subvector(in.child(), 2, 3));
          return temp;
        }
      };

      struct SpectralRadiusUncertaintyFunctor : public vw::ReturnFixedType<float> {
        float operator()(Vector<float, 3> const & in) const {
          // assume a covariance matrix of the form
          // | a  b |
          // | b  c |

          double a = in[0];
          double b = in[1];
          double c = in[2];

          // compute the eigen values of the matrix
          double desc = sqrt(a*a + 4*b*b - 2*a*c + c*c);
          double e_1 = .5*(a+c - desc);
          double e_2 = .5*(a+c + desc);

          return sqrt(std::max<float>(e_1, e_2));
        }
      };

      void set_pyramid_levels(int levels) { pyramid_levels = levels; }
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
      void set_inner_iter_max(int iter) { inner_iter_max = iter; }
      void set_inner_epsilon(double epsilon) { epsilon_inner = epsilon; }
      void set_min_determinant(double min) { affine_min_det = min; }
      void set_max_determinant(double max) { affine_max_det = max; }
      void set_debug_region(BBox2i r) { debug_region = r; }


      Vector2i kernel_size() const { return m_kernel_size; }

      // Standard ImageView interface methods
      inline int32 cols() const { return m_left_image.cols(); }
      inline int32 rows() const { return m_left_image.rows(); }
      inline int32 planes() const { return 1; }

      inline pixel_accessor origin() const { return pixel_accessor( *this, 0, 0 ); }

      inline pixel_type operator()(double /*i*/, double /*j*/, int32 /*p*/ = 0) const {
        vw_throw(NoImplErr() << "CorrelatorView::operator()(double i, double j, int32 p) has not been implemented.");
        return pixel_type();
      }

      /// \cond INTERNAL
      typedef CropView<ImageView<result_type> > prerasterize_type;
      inline prerasterize_type prerasterize(BBox2i const& bbox) const;

      template <class DestT> inline void rasterize(DestT const& dest, BBox2i const& bbox) const {
        vw::rasterize(prerasterize(bbox), dest, bbox);
      }

    private:
      // Image references
      ImageViewRef<ImagePixelT> m_left_image, m_right_image;
      ImageViewRef<disparity_pixel> m_course_disparity;

      // global settings
      int pyramid_levels;

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

      int debug_level;
      BBox2i debug_region;

      // private helper methods
      template <class ImageT, class DisparityT1, class DisparityT2, class AffineT>
        inline void
        m_subpixel_refine(ImageViewBase<ImageT> const& left_image, ImageViewBase<ImageT> const& right_image,
                          ImageViewBase<DisparityT1> &disparity_in, ImageViewBase<DisparityT2> &disparity_out,
                          ImageViewBase<AffineT> & affine_warps,
                          BBox2i const& ROI, bool final, bool debug = false) const;
    }; // end class


    /// Summarize a CorrelatorView object
    template <class ImagePixelT>
      inline std::ostream& operator<<( std::ostream& os, EMSubpixelCorrelatorView<ImagePixelT> const& view );

  }
}

#include "EMSubpixelCorrelatorView.hpp"

#endif
