// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// \file LDRtoHDR.h
///
/// Functions for stitching multiple LDR images of the same scene into
/// one HDR image.
///
/// These functions will either take or return camera response curve
/// coefficients.  The camera response curve is approximated by an
/// N-th order polynomial:
///
/// psi(x) = sum_i ( a_i * x^i )  for i = [0..N]
///
/// The camera response coefficients are stored in a std::vector of
/// vw::Vector<double> data types.  The number of components in the
/// std::vector must match the number of channels in the image.
///

#ifndef __VW_HDR_LDRTOHDR_H__
#define __VW_HDR_LDRTOHDR_H__

#include <vw/Image/ImageViewRef.h>
#include <vw/HDR/CameraCurve.h>

#include <vector>

namespace vw {
namespace hdr {

  /// Converts each pixel value in the image to scaled illuminance
  /// values based on a set of polynomial response curves, one
  /// curve for each image channel.
  template <class SrcPixelT>
  class HighDynamicRangeView : public ImageViewBase<HighDynamicRangeView<SrcPixelT> > {

  public:

    typedef typename CompoundChannelCast<SrcPixelT, double>::type pixel_type;
    typedef pixel_type result_type;
    typedef ProceduralPixelAccessor<HighDynamicRangeView> pixel_accessor;

  private:
    std::vector<ImageViewRef<SrcPixelT> > m_views;
    CameraCurveFn m_curves;
    std::vector<double> m_brightness_vals;

  public:

    HighDynamicRangeView( std::vector<ImageViewRef<SrcPixelT> > const& views,
                          CameraCurveFn const& curves,
                          std::vector<double> brightness_vals) :
      m_views(views), m_curves(curves), m_brightness_vals(brightness_vals) {
    }

    inline int32 cols() const { return m_views[0].cols(); }
    inline int32 rows() const { return m_views[0].rows(); }
    inline int32 planes() const { return m_views[0].planes(); }

    inline pixel_accessor origin() const { return pixel_accessor(*this); }

    inline pixel_type operator()( int i, int j, int p=0 ) const {
      pixel_type hdr_pix;
      double weight_sum = 0;

      // Bring all images into same domain and average pixels across images using
      // a weighting function that favors pixels in middle of dynamic range.
      for ( unsigned c = 0; c < m_views.size(); ++c ) {

        // We will use a gaussian weighting scheme that peak at 0.5 and
        // falls off to very close to zero at 0.0 and 1.0.
        //
        // Although it was constructed by "eyeballing it" in MATLAB,
        // we find that this scheme works very well in practice.  It
        // is certainly better than our old "linear" weighting scheme:
        //
        //        double weight = 2.0 * (-abs(0.5 - gray) + 0.5);
        //
        // Gaussian Weighting:
        PixelGray<double> gray(m_views[c](i,j,p));  // Convert to grayscale
        double weight = exp(-pow((gray-0.5),2)/(0.07));

        // The camera response function returns a relative luminance
        // value between 0.0 and 2.0.
        pixel_type src_val = m_curves( m_views[c](i,j,p) );
        hdr_pix += weight * m_brightness_vals[c] * src_val;
        weight_sum += weight;
      }

      // Divide by sum of weights
      return hdr_pix / weight_sum;


      // NOTE: This is code was used to handle the saturation cases,
      // but it doesn't seem to be necessary with the new weighting
      // function.  It should be deleted at some point in the future.
//       if (weight_sum > 0.01)
//         // Divide by sum of weights
//         return hdr_pix / weight_sum;
//        else {
//         int nth_view = m_views.size()-1;
//         pixel_type bracket_val1 = m_curves(m_views[0](i,j,p) );
//         pixel_type bracket_val2 = m_curves(m_views[nth_view](i,j,p) );
//         double gray1 = 2.0*m_brightness_vals[0]*PixelGray<double>(bracket_val1);  // Convert to grayscale
//         double gray2 = 2.0*m_brightness_vals[nth_view]*PixelGray<double>(bracket_val2);  // Convert to grayscale

//         // For pixels that lie at the very fringes of the dynamic
//         // range, the weight will be very close to zero, so dividing
//         // by the weight will lead to a poorly conditioned average
//         // value as a result.  For these values, we instead take the
//         // pixel values directly from the brightest or darkest image.
//         //
//         // Case 1: Fully saturated (bright pixel)
//         if (gray1 > m_brightness_vals[0] && gray2 > m_brightness_vals[nth_view]) {
//           if (gray1 > gray2)
//             return 2.0*m_brightness_vals[0]*bracket_val1;
//           else
//             return 2.0*m_brightness_vals[nth_view]*bracket_val2;

//         // Case 2: Full unsaturated (black pixel)
//         } else {
//           if (gray1 < gray2)
//             return 2.0*m_brightness_vals[0]*bracket_val1;
//           else
//             return 2.0*m_brightness_vals[nth_view]*bracket_val2;
//         }

//       }
    }

    /// \COND INTERNAL
    typedef HighDynamicRangeView<SrcPixelT> prerasterize_type;
    inline prerasterize_type prerasterize( BBox2i const& bbox ) const {

      // Rasterize each of the individual LDR images if necessary.
      std::vector<ImageViewRef<SrcPixelT> > m_prerast_views(m_views.size());
      for (unsigned i = 0; i < m_views.size(); ++i)
        m_prerast_views[i] = m_views[i].prerasterize(bbox);

      return prerasterize_type( m_prerast_views, m_curves, m_brightness_vals );
    }
    template <class DestT> inline void rasterize( DestT const& dest, BBox2i const& bbox ) const { vw::rasterize( prerasterize(bbox), dest, bbox ); }
    /// \endcond
  };

//     // First pick a reasonable middle threshold value for picking
//     // whether a pixel is saturated as white or black.
//     double min, max;
//     ImageView<PixelGray<double> > bright_or_dark = LuminanceView<ViewT>(images[images.size() / 2], curves, brightness_values[images.size() / 2]);
//     min_max_channel_values(bright_or_dark, min, max);
//     double threshold = (min + max) / 2;

//     // Now, determine the absolute minimum and maximum values present
//     // in the image.
//     ImageView<PixelGray<double> > gray = hdr_image;
//     min_max_channel_values(gray, min, max);

//     // Assign very bright/dark pixels
//     for ( unsigned i = 0; i < sat.size(); ++i ) {
//       int col = sat[i].x();
//       int row = sat[i].y();
//       double value;
//       if (bright_or_dark(col, row)[0] > threshold) {
//         value = max;
//       } else {
//         value = min;
//       }
//       for ( int32 channel = 0; channel < n_channels; ++channel ) {
//         hdr_image(col, row)[channel] = value;
//       }
//     }

  /// Stitches the set of LDR images into one HDR image.  The camera
  /// response curve is estimated from the provided images, however
  /// you must supply a vector of brightness values corresponding to
  /// the sensitivity of the exposure of each image.  (only the ratio
  /// of the brightness values is important here, although you are
  /// best off working in physical units of lumens, etc.)  Brightness
  /// is a function of camera ISO (sensitivity or sensor gain),
  /// shutter speed, and aperture.
  ///
  /// This code roughly implements Mitsunaga nod Nayar technique for
  /// reconstructing HDR images from a bracketed set of LDR
  /// exposures. See:
  ///
  /// T. Mitsunaga and S.K.Nayar. "Radiometric Self Calibration,"
  /// in Proceedings of the IEEE Conference on Computer Vision and
  /// Pattern Recoognition, Fort Collins, CO. June 1999.
  ///
//   template <class ViewT>
//   ImageView<double> process_ldr_images(std::vector<ViewT> const &images,
//                                        std::vector<Vector<double> > &ret_curves,
//                                        std::vector<double> brightness_values) {

//     typedef typename PixelChannelType<PixelT>::type channel_type;

//     int32 n_channels = PixelNumChannels<PixelT>::value;
//     std::vector<Matrix<double> > pairs(n_channels);
//     std::vector<Vector<double> > curves(n_channels);

//     // Sample each image channel
//     for ( int32 i = 0; i < n_channels; ++i ) {
//       pairs[i] = generate_ldr_intensity_pairs(images, brightness_values, VW_HDR_NUM_PAIRS, i);
//     }
//     vw_out() << "Brightness values: \n";
//     for (int i = 0; i < brightness_values.size(); ++i) {
//       vw_out() << brightness_values[i] << "\n";
//     }

//      std::cout << pairs[0].cols() << " x " << pairs[0].rows() << "\n";

//     for (int i = 0; i < pairs[0].rows(); ++i) {
//       vw_out() << pairs[0](i,0) << " " << pairs[0](i,1) << " " <<  pairs[0](i,2) << "\n";
//     }

//     // Compute camera response curve for each channel. See CameraCurve.h.
//     for ( int32 i = 0; i < n_channels; ++i ) {
//       estimate_inverse_camera_curve(pairs[i], curves[i], VW_HDR_RESPONSE_POLYNOMIAL_ORDER);
//     }

//     ret_curves = curves;

//     // Stitch LDR image set into HDR image
//     return ldr_to_hdr(images, curves, brightness_values);
//   }

}} // vw::hdr

#endif  // __LDRtoHDR_H__
