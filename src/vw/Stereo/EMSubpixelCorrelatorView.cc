#include <vw/Stereo/EMSubpixelCorrelatorView.h>

namespace vw {
namespace stereo {
namespace detail {

// disparity map down-sampling by two
ImageView<PixelMask<Vector2f> >
subsample_disp_map_by_two(ImageView<PixelMask<Vector2f> > const& input_disp)  {

  // determine the new size of the image
  int new_width = input_disp.cols()/2;
  int new_height = input_disp.rows()/2;

  ImageView<PixelMask<Vector2f> > outDisp(new_width, new_height);
  ImageViewRef<PixelMask<Vector2f> > disp = edge_extend(input_disp, ConstantEdgeExtension() );

  for (vw::int32 j = 0; j < outDisp.rows(); j++) {
    for (vw::int32 i = 0; i < outDisp.cols(); i++) {

      vw::int32 old_i = i*2;
      vw::int32 old_j = j*2;

      PixelMask<Vector2f> sum(0,0);
      int num_valid = 0;

      if ( is_valid(disp(old_i,old_j)) ) {
        sum += disp(old_i,old_j);
        num_valid++;
      }
      if ( is_valid(disp(old_i,old_j+1)) ) {
        sum += disp(old_i,old_j+1);
        num_valid++;
      }
      if ( is_valid(disp(old_i+1,old_j)) ) {
        sum += disp(old_i+1,old_j);
        num_valid++;
      }
      if ( is_valid(disp(old_i+1,old_j+1)) ) {
        sum += disp(old_i+1,old_j+1);
        num_valid++;
      }

      if (num_valid == 0)
        invalidate( outDisp(i,j) );
      else
        outDisp(i,j) = sum/float(2*num_valid);
    }
  }

  return outDisp;
}

// disparity map up-sampling by two
ImageView<PixelMask<Vector2f> >
upsample_disp_map_by_two(ImageView<PixelMask<Vector2f> > const& input_disp,
                         int up_width, int up_height) {

  ImageView<PixelMask<Vector2f> > outDisp(up_width, up_height);
  ImageViewRef<PixelMask<Vector2f> > disp = edge_extend(input_disp, ConstantEdgeExtension());

  for (vw::int32 j = 0; j < outDisp.rows(); ++j) {
    for (vw::int32 i = 0; i < outDisp.cols(); ++i) {
      float x = math::impl::_floor(float(i)/2.0), y = math::impl::_floor(float(j)/2.0);

      if ( i%2 == 0 && j%2 == 0) {
        if ( !is_valid(disp(x,y)) )
          invalidate( outDisp(i,j) );
        else
          outDisp(i,j) = 2 * disp(x,y);
      }

      else if (j%2 == 0) {
        if ( !is_valid(disp(x,y)) && !is_valid(disp(x+1,y)) )
          invalidate(outDisp(i,j));
        else {
          if ( is_valid(disp(x,y)) && is_valid(disp(x+1,y)) ) {
            outDisp(i,j) = disp(x,y) + disp(x+1,y);
          } else if ( is_valid(disp(x,y)) && !is_valid(disp(x+1,y)) ) {
            outDisp(i,j) = 2 * disp(x,y);
          } else if ( !is_valid(disp(x,y)) && is_valid(disp(x+1,y)) ) {
            outDisp(i,j) = 2 * disp(x+1,y);
          }
        }
      }

      else if (i%2 == 0) {
        if ( !is_valid(disp(x,y)) && !is_valid(disp(x,y+1)) )
          invalidate( outDisp(i,j) );
        else {
          if ( is_valid(disp(x,y)) && is_valid(disp(x,y+1)) ) {
            outDisp(i,j) = disp(x,y) + disp(x,y+1);
          } else if ( is_valid(disp(x,y)) && !is_valid(disp(x,y+1)) ) {
            outDisp(i,j) = 2*disp(x,y);
          } else if ( !is_valid(disp(x,y)) && is_valid(disp(x,y+1)) ) {
            outDisp(i,j) = 2*disp(x,y+1);
          }
        }
      }

      else {
        if ( is_valid(disp(x,y)) && is_valid(disp(x,y+1)) &&
             is_valid(disp(x+1,y)) && is_valid(disp(x+1,y+1)) ) {

          // All good pixels
          float normx = float(i)/2.0-x, normy = float(j)/2.0-y, norm1mx = 1.0-normx, norm1my = 1.0-normy;
          outDisp(i,j) = PixelMask<Vector2f>( 2 * (disp(x  ,y  )[0] * norm1mx*norm1my +
                                                   disp(x+1,y  )[0] * normx*norm1my +
                                                   disp(x  ,y+1)[0] * norm1mx*normy +
                                                   disp(x+1,y+1)[0] * normx*normy),
                                              2 * (disp(x  ,y  )[1] * norm1mx*norm1my +
                                                   disp(x+1,y  )[1] * normx*norm1my +
                                                   disp(x  ,y+1)[1] * norm1mx*normy +
                                                   disp(x+1,y+1)[1] * normx*normy) );

        } else if ( !is_valid(disp(x,y)) && !is_valid(disp(x,y+1)) &&
                    !is_valid(disp(x+1,y)) && !is_valid(disp(x+1,y+1)) ) {
          // no good pixels
          invalidate( outDisp(i,j) );

        } else {
          // some good & some bad pixels
          //
          // We fudge things a little bit here by picking the
          // first good pixel.  This isn't exactly correct, but
          // it's close enough in the handful of cases where we
          // are near some missing pixels, and we just need an
          // approximately valid value.  The subpixel refinement
          // will correct any minor mistakes we introduce here.
          if ( is_valid(disp(x,y)) ) {
            outDisp(i,j) = 2*disp(x,y);
          } else if ( is_valid(disp(x+1,y)) ) {
            outDisp(i,j) = 2*disp(x+1,y);
          } else if ( is_valid(disp(x,y+1)) ) {
            outDisp(i,j) = 2*disp(x,y+1);
          } else if ( is_valid(disp(x+1,y+1)) ) {
            outDisp(i,j) = 2*disp(x+1,y+1);
          }
        }
      }
    }
  }

  return outDisp;
}

}}} // namespace vw::stereo::detail
