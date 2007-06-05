#ifndef __VW_STEREO_CORRELATE_H__
#define __VW_STEREO_CORRELATE_H__

#include <vw/Core/FundamentalTypes.h>
#include <vw/Math/Matrix.h>
#include <vw/Stereo/DisparityMap.h>

#define VW_STEREO_MISSING_PIXEL -32000

namespace vw {
namespace stereo {

VW_DEFINE_EXCEPTION(CorrelatorErr, vw::Exception);

// Timing (profiling) utilites
static inline double Time(void) {
  struct timeval tv;
  gettimeofday(&tv, 0);
  return ((double) tv.tv_sec + (double) tv.tv_usec / 1.0e6);
}

// Return absolute difference of two bit images
inline unsigned char correlator_absolute_difference(unsigned char val1, unsigned char val2) {
  return val1 ^ val2;
}

// Return absolute difference of two bit images
inline float correlator_absolute_difference(float val1, float val2) {
  float diff = val1 - val2;
  if (diff > 0) return diff; else return -diff;
}


/// Compute the sum of the absolute difference between a template
/// region taken from img1 and the window centered at (c,r) in img0.
template <class ChannelT>
inline double compute_soad(ChannelT *img0, ChannelT *img1,
                           int r, int c,                   // row and column in img0
                           int hdisp, int vdisp,           // Current disparity offset from (c,r) for img1
                           int kern_width, int kern_height,// Kernel dimensions
                           int width, int height) {        // Image dimensions
  
  r -= kern_height/2;
  c -= kern_width/2;
  if (r<0         || c<0       || r+kern_height>=height       || c+kern_width>=width ||
      r+vdisp < 0 || c+hdisp<0 || r+vdisp+kern_height>=height || c+hdisp+kern_width>=width) {
    return VW_STEREO_MISSING_PIXEL;
  }

  ChannelT *new_img0 = img0;
  ChannelT *new_img1 = img1;

  new_img0 += c + r*width;
  new_img1 += (c+hdisp) + (r+vdisp)*width;
  
  typename vw::AccumulatorType<ChannelT>::type ret = 0;
  for (int rr= 0; rr< kern_height; rr++) {
    for (int cc= 0; cc< kern_width; cc++) {
      ret += correlator_absolute_difference(new_img0[cc], new_img1[cc]);
    }
    new_img0 += width;
    new_img1 += width;
  }
  return double(ret);
}

/// For a given set of images, compute the optimal disparity (minimum
/// SOAD) at position left_image(i,j) for the given correlation window
/// settings.
/// 
/// The left_image and right_image must have the same dimensions, but
/// this is only checked here if debugging is enabled.
template <class ChannelT>
inline PixelDisparity<float> compute_disparity(ImageView<ChannelT> &left_image,
                                               ImageView<ChannelT> &right_image,
                                               int i, int j,
                                               int kern_width, int kern_height,
                                               int min_h_disp, int max_h_disp,
                                               int min_v_disp, int max_v_disp) {

  const double default_soad = 1.0e10;     // Impossibly large value
  double min_soad = default_soad;
  PixelDisparity<float> best_disparity; // Starts as a missing pixel
  for (int ii = min_h_disp; ii <= max_h_disp; ++ii) {
    for (int jj = min_v_disp; jj <= max_v_disp; ++jj) {
      double soad = compute_soad(&(left_image(0,0)), &(right_image(0,0)),
                                 j, i, ii, jj,kern_width, kern_height, 
                                 left_image.cols(), left_image.rows());
      if (soad != VW_STEREO_MISSING_PIXEL && soad < min_soad) {
        min_soad = soad;
        best_disparity = PixelDisparity<float>(ii, jj);
      }
    }
  }
  return best_disparity;
}

static double find_minimum(double lt, double mid, double rt) {
  double a = (rt+lt)*0.5-mid;
  double b = (rt-lt)*0.5;
  return -b/(2.0*a);
}

/* 
 * Find the minimun of a 2d hyperbolic surface that is fit to the nine points 
 * around and including the peak in the disparity map.  This gives better 
 * subpixel resolution when both horizontal and vertical subpixel is requested.
 * 
 * The equation of the surface we are fitting is:
 *    z = ax^2 + by^2 + cxy + dx + ey + f
 */
template <class VectorT, class MatrixT>
static vw::Vector2 find_minimum_2d(vw::VectorBase<VectorT> &points, vw::MatrixBase<MatrixT> &pinvA) {

  vw::Vector2 offset;

  /* 
   * First, compute the parameters of the hyperbolic surface by fitting the nine points in 'points'
   * using a linear least squares fit.  This process is fairly fast, since we have already pre-computed
   * the inverse of the A matrix in Ax = b.
   */
  vw::Vector<double> x = pinvA * points;
  
  /* 
   * With these parameters, we have a closed form expression for the surface.  We compute the 
   * derivative, and find the point where the slope is zero.  This is our maximum.
   *
   * Max is at [x,y] where:
   *
   *   dz/dx = 2ax + cy + d = 0
   *   dz/dy = 2by + cx + e = 0
   * 
   * Of course, we optimize this computation a bit by unrolling it by hand beforehand.
   */
  double denom = 4 * x(0) * x(1) - (x(2) * x(2));
  
  offset(0) = ( x(2) * x(4) - 2 * x(1) * x(3) ) / denom;
  offset(1) = ( x(2) * x(3) - 2 * x(1) * x(4) ) / denom;

  return offset;
}

template <class ChannelT> 
void subpixel_correlation(ImageView<PixelDisparity<float> > &disparity_map,
                          ImageView<ChannelT> const& left_image,
                          ImageView<ChannelT> const& right_image,
                          int kern_width, int kern_height,
                          bool do_horizontal_subpixel = true,
                          bool do_vertical_subpixel = true) {

  VW_ASSERT(left_image.cols() == right_image.cols() && left_image.cols() == disparity_map.cols() &&
            left_image.rows() == right_image.rows() && left_image.rows() == disparity_map.rows(),
            ArgumentErr() << "subpixel_correlation: input image dimensions do not agree.\n");

  int height = disparity_map.rows();
  int width = disparity_map.cols();

  ChannelT *new_img0 = &(left_image(0,0));
  ChannelT *new_img1 = &(right_image(0,0));
  
  // Bail out if no subpixel computation has been requested 
  if (!do_horizontal_subpixel && !do_vertical_subpixel) return;
  
  // We get a considerable speedup in our 2d subpixel correlation if
  // we go ahead and compute the pseudoinverse of the A matrix (where
  // each row in A is [ x^2 y^2 xy x y 1] (our 2d hyperbolic surface)
  // for the range of x = [-1:1] and y = [-1:1].
  static double pinvA_data[] = { 1.0/6,  1.0/6,  1.0/6, -1.0/3, -1.0/3, -1.0/3,  1.0/6,  1.0/6,  1.0/6,
                                 1.0/6, -1.0/3,  1.0/6,  1.0/6, -1.0/3,  1.0/6,  1.0/6, -1.0/3,  1.0/6,
                                 1.0/4,    0.0, -1.0/4,    0.0,    0.0,    0.0, -1.0/4,    0.0,  1.0/4,
                                 -1.0/6, -1.0/6, -1.0/6,    0.0,    0.0,   0.0,  1.0/6,  1.0/6,  1.0/6,
                                 -1.0/6,    0.0,  1.0/6, -1.0/6,    0.0, 1.0/6, -1.0/6,    0.0,  1.0/6,
                                 -1.0/9,  2.0/9, -1.0/9,  2.0/9,  5.0/9, 2.0/9, -1.0/9,  2.0/9, -1.0/9 }; 
  vw::MatrixProxy<double,6,9> pinvA(pinvA_data);
  for (int r = 0; r < height; r++) {
    if (r%100 == 0) {
      printf("\tPerforming sub-pixel correlation... %0.2f%%\r", double(r)/height * 100);
      fflush(stdout);
    }
    
    for (int c = 0; c < width; c++) {
      
      if ( !disparity_map(c,r).missing() ) {
        int hdisp= (int)disparity_map(c,r).h();
        int vdisp= (int)disparity_map(c,r).v();
        
        double mid = compute_soad(new_img0, new_img1,
                                  r, c,
                                  hdisp,   vdisp,
                                  kern_width, kern_height,
                                  width, height);
        
        // If only horizontal subpixel resolution is requested 
        if (do_horizontal_subpixel && !do_vertical_subpixel) {
        double lt= compute_soad(new_img0, new_img1,
                                  r, c,
                                  hdisp-1, vdisp,
                                  kern_width, kern_height,
                                  width, height);
          double rt= compute_soad(new_img0, new_img1,
                                  r, c,
                                  hdisp+1, vdisp,
                                  kern_width, kern_height,
                                  width, height);
          
          if ((mid <= lt && mid < rt) || (mid <= rt && mid < lt)) {
            disparity_map(c,r).h() += find_minimum(lt, mid, rt);
          } else {
            disparity_map(c,r) = PixelDisparity<float>();
          }
        }
      
        // If only vertical subpixel resolution is requested 
        if (do_vertical_subpixel && !do_horizontal_subpixel) {
          double up= compute_soad(new_img0, new_img1,
                                  r, c,
                                  hdisp, vdisp-1,
                                  kern_width, kern_height,
                                  width, height);
          double dn= compute_soad(new_img0, new_img1,
                                  r, c,
                                  hdisp, vdisp+1,
                                  kern_width, kern_height,
                                  width, height);
          
          if ((mid <= up && mid < dn) || (mid <= dn && mid < up)) {
            disparity_map(c,r).v() += find_minimum(up, mid, dn);
          } else {
            disparity_map(c,r) = PixelDisparity<float>();
          }
        }
        
        
        // If both vertical and horizontal subpixel resolution is requested,
        // we try to fit a 2d hyperbolic surface using the 9 points surrounding the
        // peak SOAD value.  
        //
        // We place the soad values into a vector using the following indices
        // (i.e. index 4 is the max disparity value)
        // 
        //     0  3  6
        //     1  4  7
        //     2  5  8
        //
        if (do_vertical_subpixel && do_horizontal_subpixel) {
          vw::Vector<double,9> points;
          
          points(0) = (double)compute_soad(new_img0, new_img1,
                                           r, c,
                                           hdisp-1, vdisp-1,
                                           kern_width, kern_height,
                                           width, height);
          points(1) = (double)compute_soad(new_img0, new_img1,
                                           r, c,
                                           hdisp-1, vdisp,
                                           kern_width, kern_height,
                                           width, height);
          points(2) = (double)compute_soad(new_img0, new_img1,
                                           r, c,
                                           hdisp-1, vdisp+1,
                                           kern_width, kern_height,
                                           width, height);
          points(3) = (double)compute_soad(new_img0, new_img1,
                                           r, c,
                                           hdisp, vdisp-1,
                                           kern_width, kern_height,
                                           width, height);
          points(4) = (double)mid;
          points(5) = (double)compute_soad(new_img0, new_img1,
                                           r, c,
                                           hdisp, vdisp+1,
                                           kern_width, kern_height,
                                           width, height);
          points(6) = (double)compute_soad(new_img0, new_img1,
                                           r, c,
                                           hdisp+1, vdisp-1,
                                           kern_width, kern_height,
                                           width, height);
          points(7) = (double)compute_soad(new_img0, new_img1,
                                           r, c,
                                           hdisp+1, vdisp,
                                           kern_width, kern_height,
                                           width, height);
          points(8) = (double)compute_soad(new_img0, new_img1,
                                           r, c,
                                           hdisp+1, vdisp+1,
                                           kern_width, kern_height,
                                           width, height);
            
          vw::Vector2 offset = find_minimum_2d(points, pinvA);
          
          // This prevents us from adding in large offsets for
          // poorly fit data.
          if (fabs(offset(0)) < 3.0 && fabs(offset(1)) < 3.0) {
            disparity_map(c,r).h() += offset(0);
            disparity_map(c,r).v() += offset(1);
          } else {
            disparity_map(c,r) = PixelDisparity<float>();
            //            std::cout << "Bad offset: " << offset(0) << " " << offset(1) << "\n";
          }
        }
      } 
    }
  }
  printf("\tPerforming sub-pixel correlation... Done.           \n");

}

/// This routine cross checks L2R and R2L, placing the final version
/// of the disparity map in L2R.
static void cross_corr_consistency_check(ImageView<PixelDisparity<float> > &L2R, 
                                         ImageView<PixelDisparity<float> > &R2L,
                                         double cross_corr_threshold, bool verbose = false) {
  int32 xx,yy;
  int xOffset, yOffset;
  int count = 0, match_count = 0;
  
  if (verbose)
    vw_out(InfoMessage) << "\tCrosscorr threshold: " << cross_corr_threshold << "\n";
  if (cross_corr_threshold < 0) 
    vw_throw( vw::ArgumentErr() << "CrossCorrConsistencyCheck2D: the crosscorr threshold was less than 0." );
  
  for(xx = 0; xx < L2R.cols(); xx++) {     
    for(yy = 0; yy < L2R.rows(); yy++) {
      
      int xOffset = (int)L2R(xx,yy).h();
      int yOffset = (int)L2R(xx,yy).v();
      
      // Check to make sure we are within the image bounds
      if(xx+xOffset < 0 || yy+yOffset < 0 ||
         xx+xOffset >= R2L.cols() || yy+yOffset >= R2L.rows()) {
        L2R(xx,yy) = PixelDisparity<float>();  // Default constructor is missing pixel.
      }

      // Check for missing pixels
      else if ( L2R(xx,yy).missing() ||
                R2L(xx+xOffset, yy+yOffset).missing() ) {
        L2R(xx,yy) = PixelDisparity<float>();  // Default constructor is missing pixel.
      }
      
      // Check for correlation consistency
      //
      // Since the hdisp for the R2L and L2R buffers will be opposite 
      // in sign, we determine their similarity by *summing* them, rather
      // than differencing them as you might expect.
      else if (cross_corr_threshold >= fabs(L2R(xx,yy).h() + R2L(xx+xOffset,yy+yOffset).h()) &&
               cross_corr_threshold >= fabs(L2R(xx,yy).v() + R2L(xx+xOffset,yy+yOffset).v())) {
        count++;
        match_count++;
      }
      
      // Otherwise, the pixel is bad.
      else {
        match_count++;
        L2R(xx,yy) = PixelDisparity<float>();  // Default constructor is missing pixel.
      } 
    }
  } 
  if (verbose) 
    vw_out(InfoMessage) << "\tCross-correlation retained " << count << " / " << match_count << " matches (" << ((float)count/match_count*100) <<" percent).\n";
}

}} // namespace vw::stereo

#endif // __VW_STEREO_CORRELATOR_H__
