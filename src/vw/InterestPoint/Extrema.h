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

/// \file Extrema.h
/// 
/// Functions for finding local extrema.
/// 
#ifndef __INTERESTPOINT_EXTREMA_H__
#define __INTERESTPOINT_EXTREMA_H__

#include <vw/InterestPoint/Detector.h>
#include <vw/InterestPoint/ImageOctave.h>

#include <vw/FileIO.h>

namespace vw {
namespace ip {

#define IP_BORDER_WIDTH (20)

  enum PeakType { IP_MAX, IP_MIN, IP_MINMAX };

  template <class T>
  struct ImageInterestData;

  /* Local min/max functions. This is a bit verbose. */

  // Search for spatial extrema only
  template <class T>
  inline bool is_local_max(const ImageView<T>& interest, int i, int j) {
    // 4-connected neighborhood in the plane
    if (interest(i,j) <= interest(i,j-1)) return false;// up
    if (interest(i,j) <= interest(i,j+1)) return false;// down

    if (interest(i,j) <= interest(i-1,j)) return false;// left
    if (interest(i,j) <= interest(i+1,j)) return false;// right

    // The rest of the 8-connected neighborhood in the plane
    if (interest(i,j) <= interest(i-1,j-1)) return false;
    if (interest(i,j) <= interest(i-1,j+1)) return false;
    if (interest(i,j) <= interest(i+1,j-1)) return false;
    if (interest(i,j) <= interest(i+1,j+1)) return false;

    return true;
  }

  template <class T>
  inline bool is_local_min(const ImageView<T>& interest, int i, int j) {
    // 4-connected neighborhood in the plane
    if (interest(i,j) >= interest(i,j-1)) return false;// up
    if (interest(i,j) >= interest(i,j+1)) return false;// down

    if (interest(i,j) >= interest(i-1,j)) return false;// left
    if (interest(i,j) >= interest(i+1,j)) return false;// right

    // The rest of the 8-connected neighborhood in the plane
    if (interest(i,j) >= interest(i-1,j-1)) return false;
    if (interest(i,j) >= interest(i-1,j+1)) return false;
    if (interest(i,j) >= interest(i+1,j-1)) return false;
    if (interest(i,j) >= interest(i+1,j+1)) return false;

    return true;
  }

  template <class T>
  inline bool is_local_minmax(const ImageView<T>& interest, int i, int j) {
    return (is_local_max(interest, i, j) || is_local_min(interest, i, j));
  }

  // Search for spacial and scale space extrema
  template <class T>
  inline bool is_local_max(std::vector<ImageInterestData<T> >& data, int i, int j, int k) {
    // These checks are in order such that hopefully the average
    // number of checks done is small.

    // 4-connected neighborhood in the plane
    if (data[k].interest(i,j) <= data[k].interest(i,j-1)) return false;// up
    if (data[k].interest(i,j) <= data[k].interest(i,j+1)) return false;// down

    if (data[k].interest(i,j) <= data[k].interest(i-1,j)) return false;// left
    if (data[k].interest(i,j) <= data[k].interest(i+1,j)) return false;// right

    // Same pixel location in planes above and below
    if (data[k].interest(i,j) <= data[k-1].interest(i,j)) return false;// below
    if (data[k].interest(i,j) <= data[k+1].interest(i,j)) return false;// above

    // The rest of the 8-connected neighborhood in the plane
    if (data[k].interest(i,j) <= data[k].interest(i-1,j-1)) return false;
    if (data[k].interest(i,j) <= data[k].interest(i-1,j+1)) return false;
    if (data[k].interest(i,j) <= data[k].interest(i+1,j-1)) return false;
    if (data[k].interest(i,j) <= data[k].interest(i+1,j+1)) return false;

    // TODO: Add 8-connected check in adjacent planes?

    return true;
  }

  template <class T>
  inline bool is_local_min(std::vector<ImageInterestData<T> >& data, int i, int j, int k) {
    // These checks are in order such that hopefully the average
    // number of checks done is small.

    // 4-connected neighborhood in the plane
    if (data[k].interest(i,j) >= data[k].interest(i,j-1)) return false;// up
    if (data[k].interest(i,j) >= data[k].interest(i,j+1)) return false;// down

    if (data[k].interest(i,j) >= data[k].interest(i-1,j)) return false;// left
    if (data[k].interest(i,j) >= data[k].interest(i+1,j)) return false;// right

    // Same pixel location in planes above and below
    if (data[k].interest(i,j) >= data[k-1].interest(i,j)) return false;// below
    if (data[k].interest(i,j) >= data[k+1].interest(i,j)) return false;// above

    // The rest of the 8-connected neighborhood in the plane
    if (data[k].interest(i,j) >= data[k].interest(i-1,j-1)) return false;
    if (data[k].interest(i,j) >= data[k].interest(i-1,j+1)) return false;
    if (data[k].interest(i,j) >= data[k].interest(i+1,j-1)) return false;
    if (data[k].interest(i,j) >= data[k].interest(i+1,j+1)) return false;

    // TODO: Add 8-connected check in adjacent planes?

    return true;
  }

  template <class T>
  inline bool is_local_minmax(std::vector<ImageInterestData<T> >& data,
			      int i, int j, int k) {
    return (is_local_max(data, i, j, k) || is_local_min(data, i, j, k));
  }

  // Find peaks in the image
  template <class T>
  int find_peaks( std::vector<InterestPoint>& interest_points,
		  const ImageView<T>& interest, T min_interest = 0,
		  PeakType type = IP_MAX) {
    unsigned ncols = interest.cols();
    unsigned nrows = interest.rows();

    // Find local extrema
    // TODO: this really needs to be sped up
    for (unsigned j=IP_BORDER_WIDTH; j<nrows-IP_BORDER_WIDTH; j++) {   // row j
      for (unsigned i=IP_BORDER_WIDTH; i<ncols-IP_BORDER_WIDTH; i++) { // col i
	// check if it is a local extremum
	if (interest(i,j)>min_interest) {
	  if (((type == IP_MAX) && is_local_max(interest, i, j)) ||
	      ((type == IP_MIN) && is_local_min(interest, i, j)) ||
	      ((type == IP_MINMAX) && is_local_minmax(interest, i, j))) {
	    vw::vw_out(DebugMessage) << "Found a local max at [" << i << ", " << j
				     << "]" << "    Interest: " << interest(i,j) << "\n";
	  
	    InterestPoint pt;
	    pt.x = i;
	    pt.y = j;
	    pt.scale = 1.0;
	    pt.interest = interest(i,j);
	    interest_points.push_back( pt );
	  }
	}
      } // col j
    } // row i
    
    return 0;
  }


  // Find peaks in the image octave
  template <class T>
  int find_peaks( std::vector<InterestPoint>& interest_points, 
		  std::vector<ImageInterestData<T> >& data,
		  const ImageOctave<T>& octave,
		  T min_interest = 0, PeakType type = IP_MAX) {
    // Check that we have a few planes of corner response function
    // (interest) and that all planes are the same size.
    assert( data.size() > 0 );
    unsigned ncols = data[0].interest.cols();
    unsigned nrows = data[0].interest.rows();
    unsigned nplanes = data.size();
    //printf( "interest is %d by %d by %d\n", ncols, nrows, nplanes );

    // Make sure all planes are the same size
    for (unsigned k=0; k<nplanes; k++){
      assert( ncols == data[k].interest.cols() );
      assert( nrows == data[k].interest.rows() );
    }
    
    // In this implementation, we don't want to compare plane 0 to
    // plane 1, or compare plane N to plane N-1.  Each of these has a
    // scale space plane only on one side, so we cannot bound the
    // scale from both sides.  Below, search all of the internal
    // pixels on all of the internal planes for local maxima in
    // (x,y,scale)

    // Find local maxima
    // TODO: this really needs to be sped up
    for (unsigned k=1; k<nplanes-1; k++) {     // plane k
      for (unsigned j=IP_BORDER_WIDTH; j<nrows-IP_BORDER_WIDTH; j++) {   // row j
	for (unsigned i=IP_BORDER_WIDTH; i<ncols-IP_BORDER_WIDTH; i++) { // col i
	  // check if it is a local max
	  if (data[k].interest(i,j)>min_interest) {
	    if (((type == IP_MAX) && is_local_max(data, i, j, k)) ||
		((type == IP_MIN) && is_local_min(data, i, j, k)) ||
		((type == IP_MINMAX) && is_local_minmax(data, i, j, k))) {
	      vw::vw_out(DebugMessage) << "Found a local max at [" << i << ", " << j << ", "
				       << k << "]" << "    Interest: "
				       << data[k].interest(i,j) << "\n";
	      //printf("Found one at (%i, %i, %i)\n", i, j, k);
	      InterestPoint pt;
	      pt.x = i;
	      pt.y = j;
	      pt.scale = octave.plane_index_to_scale(k);
	      pt.interest = data[k].interest(i,j);
	      interest_points.push_back( pt );
	    }
	  }
	} // col j
      } // row i
    } // plane k
    
    return 0;
  }

}} // namespace vw::ip 

#endif // __INTERESTPOINT_EXTREMA_H__
