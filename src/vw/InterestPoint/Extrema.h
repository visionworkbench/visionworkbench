// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// \file Extrema.h
///
/// Functions for finding local extrema.
///
#ifndef __INTERESTPOINT_EXTREMA_H__
#define __INTERESTPOINT_EXTREMA_H__

#include <vw/InterestPoint/InterestData.h>
#include <vw/InterestPoint/ImageOctave.h>

namespace vw {
namespace ip {

// Do not search too close to image boundary.
#define IP_BORDER_WIDTH (1)

  // Local min/max functions.

  /* Functions to search for spatial extrema */

  /// Find spatial local max.
  template <class ViewT>
  inline bool is_local_max(ImageViewBase<ViewT> const& _interest, int i, int j) {
    ViewT const& interest = _interest.impl();
    typename ViewT::result_type interest_ij = interest(i,j);

    // 4-connected neighborhood in the plane
    if (interest_ij <= interest(i,j-1)) return false;// up
    if (interest_ij <= interest(i,j+1)) return false;// down

    if (interest_ij <= interest(i-1,j)) return false;// left
    if (interest_ij <= interest(i+1,j)) return false;// right

    // The rest of the 8-connected neighborhood in the plane
    if (interest_ij <= interest(i-1,j-1)) return false;
    if (interest_ij <= interest(i-1,j+1)) return false;
    if (interest_ij <= interest(i+1,j-1)) return false;
    if (interest_ij <= interest(i+1,j+1)) return false;

    return true;
  }

  /// Find spatial local min.
  template <class ViewT>
  inline bool is_local_min(ImageViewBase<ViewT> const& _interest, int i, int j) {
    ViewT const& interest = _interest.impl();
    typename ViewT::result_type interest_ij = interest(i,j);

    // 4-connected neighborhood in the plane
    if (interest_ij >= interest(i,j-1)) return false;// up
    if (interest_ij >= interest(i,j+1)) return false;// down

    if (interest_ij >= interest(i-1,j)) return false;// left
    if (interest_ij >= interest(i+1,j)) return false;// right

    // The rest of the 8-connected neighborhood in the plane
    if (interest_ij >= interest(i-1,j-1)) return false;
    if (interest_ij >= interest(i-1,j+1)) return false;
    if (interest_ij >= interest(i+1,j-1)) return false;
    if (interest_ij >= interest(i+1,j+1)) return false;

    return true;
  }

  /// Find spatial local extremum.
  template <class ViewT>
  inline bool is_local_minmax(ImageViewBase<ViewT> const& interest, int i, int j) {
    return (is_local_max(interest, i, j) || is_local_min(interest, i, j));
  }

  /// Find spatial local extremum of the specified type.
  template <class ViewT>
  inline bool is_extremum(ImageViewBase<ViewT> const& interest,
                          int i, int j, int type) {
    if (type == IP_MAX) return is_local_max(interest, i, j);
    if (type == IP_MIN) return is_local_min(interest, i, j);
    return is_local_minmax(interest, i, j);
  }

  /// Compile-time equivalent of spatial is_extremum
  template <int PeakT> struct IsSpatialExtremum {};

  template <> struct IsSpatialExtremum <IP_MAX> {
    template <class ViewT>
    static inline bool check(ImageViewBase<ViewT> const& interest, int i, int j) {
      return is_local_max(interest, i, j);
    }
  };

  template <> struct IsSpatialExtremum <IP_MIN> {
    template <class ViewT>
    static inline bool check(ImageViewBase<ViewT> const& interest, int i, int j) {
      return is_local_min(interest, i, j);
    }
  };

  template <> struct IsSpatialExtremum <IP_MINMAX> {
    template <class ViewT>
    static inline bool check(ImageViewBase<ViewT> const& interest, int i, int j) {
      return is_local_minmax(interest, i, j);
    }
  };

  /* Functions to search for scale space extrema */

  /// Find local max in space and scale.
  template <class DataT>
  inline bool is_local_max(std::vector<DataT> const& data,
                           int i, int j, int k) {
    typedef typename DataT::interest_type interest_type;
    interest_type const& interest_k = data[k].interest();
    interest_type const& interest_k_p = data[k+1].interest();
    interest_type const& interest_k_m = data[k-1].interest();
    typename interest_type::pixel_type interest_ijk = interest_k(i,j);

    // These checks are in order such that hopefully the average
    // number of checks done is small.

    // 4-connected neighborhood in the plane
    if (interest_ijk <= interest_k(i,j-1)) return false;// up
    if (interest_ijk <= interest_k(i,j+1)) return false;// down

    if (interest_ijk <= interest_k(i-1,j)) return false;// left
    if (interest_ijk <= interest_k(i+1,j)) return false;// right

    // Same pixel location in planes above and below
    if (interest_ijk <= interest_k_m(i,j)) return false;// below
    if (interest_ijk <= interest_k_p(i,j)) return false;// above

    // The rest of the 8-connected neighborhood in the plane
    if (interest_ijk <= interest_k(i-1,j-1)) return false;
    if (interest_ijk <= interest_k(i-1,j+1)) return false;
    if (interest_ijk <= interest_k(i+1,j-1)) return false;
    if (interest_ijk <= interest_k(i+1,j+1)) return false;

    // TODO: Add 8-connected check in adjacent planes?

    return true;
  }

  /// Find local min in space and scale.
  template <class DataT>
  inline bool is_local_min(std::vector<DataT> const& data,
                           int i, int j, int k) {
    typedef typename DataT::interest_type interest_type;
    interest_type const& interest_k = data[k].interest();
    interest_type const& interest_k_p = data[k+1].interest();
    interest_type const& interest_k_m = data[k-1].interest();
    typename interest_type::pixel_type interest_ijk = interest_k(i,j);

    // These checks are in order such that hopefully the average
    // number of checks done is small.

    // 4-connected neighborhood in the plane
    if (interest_ijk >= interest_k(i,j-1)) return false;// up
    if (interest_ijk >= interest_k(i,j+1)) return false;// down

    if (interest_ijk >= interest_k(i-1,j)) return false;// left
    if (interest_ijk >= interest_k(i+1,j)) return false;// right

    // Same pixel location in planes above and below
    if (interest_ijk >= interest_k_m(i,j)) return false;// below
    if (interest_ijk >= interest_k_p(i,j)) return false;// above

    // The rest of the 8-connected neighborhood in the plane
    if (interest_ijk >= interest_k(i-1,j-1)) return false;
    if (interest_ijk >= interest_k(i-1,j+1)) return false;
    if (interest_ijk >= interest_k(i+1,j-1)) return false;
    if (interest_ijk >= interest_k(i+1,j+1)) return false;

    // TODO: Add 8-connected check in adjacent planes?

    return true;
  }

  /// Find local extremum in space and scale.
  template <class DataT>
  inline bool is_local_minmax(std::vector<DataT> const& data,
                              int i, int j, int k) {
    return (is_local_max(data, i, j, k) || is_local_min(data, i, j, k));
  }

  /// Find local extremum of the specified type in space and scale.
  template <class DataT>
  inline bool is_extremum(std::vector<DataT> const& data,
                          int i, int j, int k, int type) {
    if (type == IP_MAX) return is_local_max(data, i, j, k);
    if (type == IP_MIN) return is_local_min(data, i, j, k);
    return is_local_minmax(data, i, j, k);
  }

  /// Compile-time equivalent of scale-space is_extremum
  template <int PeakT> struct IsScaleExtremum {};

  template <> struct IsScaleExtremum <IP_MAX> {
    template <class DataT>
    static inline bool check(std::vector<DataT> const& data, int i, int j, int k) {
      return is_local_max(data, i, j, k);
    }
  };

  template <> struct IsScaleExtremum <IP_MIN> {
    template <class DataT>
    static inline bool check(std::vector<DataT> const& data, int i, int j, int k) {
      return is_local_min(data, i, j, k);
    }
  };

  template <> struct IsScaleExtremum <IP_MINMAX> {
    template <class DataT>
    static inline bool check(std::vector<DataT> const& data, int i, int j, int k) {
      return is_local_minmax(data, i, j, k);
    }
  };

  /// Find spatial peaks of a certain type in the interest image.
  template <class DataT>
  int find_peaks(InterestPointList& interest_points, DataT const& data) {
    typename DataT::interest_type const& interest = data.interest();
    int32 ncols = interest.cols();
    int32 nrows = interest.rows();

    // Find local extrema
    // TODO: speed up further?
    for (int32 j=IP_BORDER_WIDTH; j<nrows-IP_BORDER_WIDTH; j++) {   // row j
      for (int32 i=IP_BORDER_WIDTH; i<ncols-IP_BORDER_WIDTH; i++) { // col i
        // check if it is a local extremum
        if (IsSpatialExtremum<DataT::peak_type>::check(interest, i, j)) {
          interest_points.push_back(InterestPoint(i, j, 1.0, interest(i,j)));
        }
      } // col j
    } // row i

    return 0;
  }

  /// Find spatial/scale peaks of some type in the image octave.
  template <class DataT, class ViewT>
  int find_peaks( InterestPointList& interest_points,
                  std::vector<DataT> const& data,
                  ImageOctave<ViewT> const& octave) {
    // Check that we have a few planes of corner response function
    // (interest) and that all planes are the same size.
    VW_ASSERT( !data.empty(), ArgumentErr() << "Data vector is empty.\n" );
    int32 ncols = data[0].interest().cols();
    int32 nrows = data[0].interest().rows();
    unsigned nplanes = data.size();

    // Make sure all planes are the same size
    for (unsigned k=0; k < nplanes; k++){
      VW_ASSERT( ncols == data[k].interest().cols(),
                 ArgumentErr() << "Planes must be the same size.\n" );
      VW_ASSERT( nrows == data[k].interest().rows(),
                 ArgumentErr() << "Planes must be the same size.\n" );
    }

    // In this implementation, we don't want to compare plane 0 to
    // plane 1, or compare plane N to plane N-1.  Each of these has a
    // scale space plane only on one side, so we cannot bound the
    // scale from both sides.  Below, search all of the internal
    // pixels on all of the internal planes for local extrema in
    // (x,y,scale)

    // Find local extrema
    // TODO: speed up further?
    for (unsigned k = 1; k < nplanes-1; k++) {     // plane k
      for (int32 j = IP_BORDER_WIDTH; j < nrows-IP_BORDER_WIDTH; j++) {   // row j
        for (int32 i = IP_BORDER_WIDTH; i < ncols-IP_BORDER_WIDTH; i++) { // col i
          // check if it is a local extremum
          if (IsScaleExtremum<DataT::peak_type>::check(data, i, j, k)) {
            interest_points.push_back(InterestPoint(i, j, octave.plane_index_to_scale(k),
                                                    data[k].interest()(i,j)));
          }
        } // col j
      } // row i
    } // plane k

    return 0;
  }

}} // namespace vw::ip

#endif // __INTERESTPOINT_EXTREMA_H__
