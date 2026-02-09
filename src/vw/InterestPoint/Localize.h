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


/// \file Localize.h
///
/// Functions for subpixel localization .
///
#ifndef __INTERESTPOINT_LOCALIZE_H__
#define __INTERESTPOINT_LOCALIZE_H__

#include <vw/InterestPoint/InterestPoint.h>
#include <vw/InterestPoint/ImageOctave.h>
#include <vw/Math/Matrix.h>
#include <vw/Math/Vector.h>

namespace vw {
namespace ip {

  /// Fits a quadratic curve to three samples to compute the
  /// subpixel shift to the fitted peak.
  template <class ElmtT>
  int fit_peak_1D(Vector<ElmtT,3> const& y, ElmtT& x_sub,
                  ElmtT* coef = NULL ) {
    // Fit the data in the vector y to the parabola a*x^2 + b*x + c
    // where x={-1,0,1}.  This is a simple 3x3 linear system of the form
    // [x_0^2, x_0, 1] [a]
    // [x_1^2, x_1, 1] [b] = y
    // [x_2^2, x_2, 1] [c]
    // which we can easily solve by inverting the coefficient matrix X.
    // Since the x coordinates are always the same, we can even
    // precompute inverse(X)
    Matrix<ElmtT,3,3> invX;
    Vector<ElmtT,3> abc;

    // We are always fitting samples at x={-1,0,1}, to compute an offset
    // which can be added to the actual x location of the peak later.
    // This means we can precompute the inverse here, and also ensures
    // better numerical stability
    invX(0,0) =  0.5; invX(0,1) = -1.0; invX(0,2) = 0.5;
    invX(1,0) = -0.5; invX(1,1) =  0.0; invX(1,2) = 0.5;
    invX(2,0) =  0.0; invX(2,1) =  1.0; invX(2,2) = 0.0;

    // Solve for polynomial coefficients (a,b,c)
    abc = invX * y;

    // Solve for subpixel shift:
    x_sub = -abc(1) / (2*abc(0));

    if (coef) *coef = abc(0);

    return 0;
  }

  /// Fit the subpixel location of the peak in the plane.
  template <class ViewT>
  bool fit_peak(ImageViewBase<ViewT> const& _interest, InterestPoint& pt,
                float* x2coef=NULL, float* y2coef=NULL) {
    ViewT const& interest = _interest.impl();
    Vector<float,3> f_vals;
    float di, dj;

    // Fit x location
    f_vals(0) = interest(pt.ix-1,pt.iy);
    f_vals(1) = interest(pt.ix,pt.iy);
    f_vals(2) = interest(pt.ix+1,pt.iy);
    fit_peak_1D( f_vals, di, x2coef );

    // Fit y location
    f_vals(0) = interest(pt.ix,pt.iy-1);
    f_vals(1) = interest(pt.ix,pt.iy);
    f_vals(2) = interest(pt.ix,pt.iy+1);
    fit_peak_1D( f_vals, dj, y2coef );

    if ( (fabs(di)>1) || (fabs(dj)>1)) {
      vw::vw_out(VerboseDebugMessage, "interest_point") << "\tLocalization correction too large\n";
      return false;
    }

    pt.x = (float)pt.ix + di;
    pt.y = (float)pt.iy + dj;

    return true;
  }

  /// Fit the subpixel location of the peak.  Stuff it into the
  /// InterestPoint object provided and return "true" if the peak is a
  /// strong one.  Return "false" if it is not.
  template <class DataT>
  bool fit_peak( std::vector<DataT> const& data,
                 InterestPoint& pt,
                 ImageOctave<typename DataT::source_type> const& octave,
                 float* x2coef=NULL, float* y2coef=NULL,
                 float* s2coef=NULL ) {
    // Fit peak by fitting only along each of the three principal
    // directions.  Besides being a lot faster, this also worked
    // better than fitting the whole 3x3x3 neighborhood using a
    // weighted regression.

    Vector<float,3> f_vals;
    float dp;
    int k0 = octave.scale_to_plane_index(pt.scale);

    // Fit (x,y) position
    if ( !fit_peak(data[k0].interest(), pt, x2coef, y2coef) )
      return false;

    // Fit scale plane location
    f_vals(0) = data[k0-1].interest()(pt.ix,pt.iy);
    f_vals(1) = data[k0].interest()(pt.ix,pt.iy);
    f_vals(2) = data[k0+1].interest()(pt.ix,pt.iy);
    fit_peak_1D( f_vals, dp, s2coef );

    if ( fabs(dp)>1) {
      vw::vw_out(VerboseDebugMessage, "interest_point") << "\tLocalization correction too large\n";
      return false;
    }

    pt.scale = octave.plane_index_to_scale((float)k0 + dp);

    // Return value indicates the peak is good enough to keep
    return true;
  }

}} // namespace vw::ip

#endif // __INTERESTPOINT_LOCALIZE_H__
