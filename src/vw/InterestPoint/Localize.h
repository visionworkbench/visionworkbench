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

/// \file Localize.h
/// 
/// Functions for subpixel localization .
/// 
#ifndef __INTERESTPOINT_LOCALIZE_H__
#define __INTERESTPOINT_LOCALIZE_H__

#include <vw/InterestPoint/InterestData.h>
#include <vw/InterestPoint/ImageOctave.h>
#include <vw/Math/Matrix.h>
#include <vw/Math/Vector.h>

namespace vw {
namespace ip {

  // Fit a quadratic to three samples.
  template <class ElmtT>
  int fit_peak_1D( const Vector<ElmtT,3>& y, ElmtT& x_sub,
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

    // Estimated parameters
    //cout << "fit_peak_1D:" << endl;
    //cout << "y = " << y << endl;
    //cout << "a = " << abc(0) << endl;
    //cout << "b = " << abc(1) << endl;
    //cout << "c = " << abc(2) << endl;
  
    // Solve for subpixel shift:
    x_sub = -abc(1) / (2*abc(0));
    //cout << "x_sub = " << x_sub << endl;
    if (NULL!=coef)
      *coef = abc(0);

    return 0;
  }

  // Fit the subpixel location of the peak in the plane.
  template <class T>
  bool fit_peak( const ImageView<T>& interest, InterestPoint& pt,
		 T* x2coef=NULL, T* y2coef=NULL) {
    Vector<T,3> f_vals;
    T di, dj;

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
      vw::vw_out(VerboseDebugMessage) << "Correction too large\n";
      return false;
    }

    pt.x += di;
    pt.y += dj;

    return true;
  }

  // Fit the subpixel location of the peak.  Stuff it into the
  // InterestPoint object provided and return "true" if the peak is a
  // strong one.  Return "false" if it is not.
  template <class T>
  bool fit_peak( std::vector<ImageInterestData<T> >& data,
		 InterestPoint& pt,
		 ImageOctave<T> const& octave,
		 T* x2coef=NULL, T* y2coef=NULL,
		 T* s2coef=NULL ) {
    // Fit peak by fitting only along each of the three principal
    // directions.  Besides being a lot faster, this also worked
    // better than fitting the whole 3x3x3 neighborhood using a
    // weighted regression.

    Vector<T,3> f_vals;
    T dp;
    int k0 = octave.scale_to_plane_index(pt.scale);

    // Fit (x,y) position
    if ( !fit_peak(data[k0].interest, pt,
		   x2coef, y2coef) )
      return false;

    // Fit scale plane location
    f_vals(0) = data[k0-1].interest(pt.ix,pt.iy);
    f_vals(1) = data[k0].interest(pt.ix,pt.iy);
    f_vals(2) = data[k0+1].interest(pt.ix,pt.iy);
    fit_peak_1D( f_vals, dp, s2coef );

    if ( fabs(dp)>1) {
      vw::vw_out(VerboseDebugMessage) << "Correction too large\n";
      return false;
    }
    
    pt.scale = octave.plane_index_to_scale((T)k0 + dp);

    // Return value indicates the peak is good enough to keep
    return true;
  }


}} // namspace vw::ip 

#endif // __INTERESTPOINT_LOCALIZE_H__
