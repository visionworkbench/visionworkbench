// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#include <vw/Stereo/Correlate.h>

#include <vw/Image/ImageView.h>
#include <vw/Image/ImageMath.h>
#include <vw/Image/PixelMask.h>

namespace vw {
namespace stereo {

  /// This routine cross checks L2R and R2L, placing the final version
  /// of the disparity map in L2R.
  void cross_corr_consistency_check(ImageView<PixelMask<Vector2f> > &L2R,
				    ImageView<PixelMask<Vector2f> > const& R2L,
				    double cross_corr_threshold, bool verbose) {
    int32 xx,yy;
    int32 count = 0, match_count = 0;

    if (verbose)
      vw_out(InfoMessage, "stereo") << "\tCrosscorr threshold: "
				    << cross_corr_threshold << "\n";
    if (cross_corr_threshold < 0)
      vw_throw( vw::ArgumentErr() << "CrossCorrConsistencyCheck2D: the crosscorr threshold was less than 0." );

    for(xx = 0; xx < L2R.cols(); xx++) {
      for(yy = 0; yy < L2R.rows(); yy++) {

	int32 xOffset = (int32)L2R(xx,yy)[0];
	int32 yOffset = (int32)L2R(xx,yy)[1];

	// Check to make sure we are within the image bounds
	if(xx+xOffset < 0 || yy+yOffset < 0 ||
	   xx+xOffset >= R2L.cols() || yy+yOffset >= R2L.rows()) {
	  invalidate( L2R(xx,yy) );
	}

	// Check for missing pixels
	else if ( !is_valid(L2R(xx,yy)) ||
		  !is_valid(R2L(xx+xOffset, yy+yOffset)) ) {
	  invalidate(L2R(xx,yy));
	}

	// Check for correlation consistency
	//
	// Since the hdisp for the R2L and L2R buffers will be opposite
	// in sign, we determine their similarity by *summing* them, rather
	// than differencing them as you might expect.
	else if (cross_corr_threshold >= fabs(L2R(xx,yy)[0] + R2L(xx+xOffset,yy+yOffset)[0]) &&
		 cross_corr_threshold >= fabs(L2R(xx,yy)[1] + R2L(xx+xOffset,yy+yOffset)[1])) {
	  count++;
	  match_count++;
	}

	// Otherwise, the pixel is bad.
	else {
	  match_count++;
	  invalidate(L2R(xx,yy));
	}
      }
    }
    if (verbose)
      vw_out(InfoMessage, "stereo") << "\tCross-correlation retained " << count << " / " << match_count << " matches (" << ((float)count/match_count*100) <<" percent).\n";
  }

}} // namespace vw::stereo
