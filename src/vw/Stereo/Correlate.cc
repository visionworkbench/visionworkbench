// __BEGIN_LICENSE__
// 
// Copyright (C) 2006 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration
// (NASA).  All Rights Reserved.
// 
// Copyright 2006 Carnegie Mellon University. All rights reserved.
// 
// This software is distributed under the NASA Open Source Agreement
// (NOSA), version 1.3.  The NOSA has been approved by the Open Source
// Initiative.  See the file COPYING at the top of the distribution
// directory tree for the complete NOSA document.
// 
// THE SUBJECT SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY WARRANTY OF ANY
// KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT
// LIMITED TO, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL CONFORM TO
// SPECIFICATIONS, ANY IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR
// A PARTICULAR PURPOSE, OR FREEDOM FROM INFRINGEMENT, ANY WARRANTY THAT
// THE SUBJECT SOFTWARE WILL BE ERROR FREE, OR ANY WARRANTY THAT
// DOCUMENTATION, IF PROVIDED, WILL CONFORM TO THE SUBJECT SOFTWARE.
// 
// __END_LICENSE__
#include <vw/Stereo/Correlate.h>

namespace vw {
namespace stereo {  

  /// This routine cross checks L2R and R2L, placing the final version
  /// of the disparity map in L2R.
  void cross_corr_consistency_check(ImageView<PixelDisparity<float> > &L2R, 
                                    ImageView<PixelDisparity<float> > &R2L,
                                    double cross_corr_threshold, bool verbose) {
    int32 xx,yy;
    int count = 0, match_count = 0;
  
    if (verbose)
      vw_out(InfoMessage, "stereo") << "\tCrosscorr threshold: " << cross_corr_threshold << "\n";
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
      vw_out(InfoMessage, "stereo") << "\tCross-correlation retained " << count << " / " << match_count << " matches (" << ((float)count/match_count*100) <<" percent).\n";
  }

}} // namespace vw::stereo
