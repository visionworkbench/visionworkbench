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

/// \file InterestData.h
///
/// Basic classes and structures for storing image interest points.
///
#ifndef __INTEREST_DATA_H__
#define __INTEREST_DATA_H__

#include <vector>
#include <list>
#include <algorithm>
#include <sstream>

#include <vw/Math/Vector.h>
#include <vw/Math/Matrix.h>
#include <vw/Math/Functors.h>
#include <vw/Image/ImageViewBase.h>
#include <vw/Image/Algorithms.h>
#include <vw/InterestPoint/InterestTraits.h>

#if defined(VW_HAVE_PKG_OPENCV) && (VW_HAVE_PKG_OPENCV == 1)
#include "opencv2/core/core.hpp"
#include "opencv2/features2d/features2d.hpp"
#endif

namespace vw { namespace ip {

  /// A class for storing information about an interest point.
  struct InterestPoint {
    typedef vw::Vector<float>               descriptor_type;
    typedef descriptor_type::iterator       iterator;
    typedef descriptor_type::const_iterator const_iterator;

    // The best way of creating interest points is:
    // ip = InterestPoint(x, y). But in the worst case, when the user
    // chooses to simply create a blank interest point, like
    // InterestPoint ip;
    // at least make sure that its members are always initialized,
    // as seen in the constructor below.
    // TODO: There is no way to enforce that ix be in sync with x or
    // iy with y.
    InterestPoint(float x = 0, float y = 0, float scale=1.0, float interest=0.0, float ori=0.0,
                  bool pol=false, uint32 octave = 0, uint32 scale_lvl = 0)
      : x(x), y(y), scale(scale), ix(int32(x)), iy(int32(y)), orientation(ori), interest(interest),
        polarity(pol), octave(octave), scale_lvl(scale_lvl) {}

    /// Subpixel (col,row) location of point
    float x,y;

    /// Scale of point.  This may come from the pyramid level, from
    /// interpolating the interest function between levels, or from some
    /// other scale detector like the Laplace scale used by Mikolajczyk & Schmid
    float scale;

    /// Integer location, mainly for internal use.
    int32 ix;
    int32 iy;

    /// Since the orientation is not necessarily unique, we may have more
    /// than one hypothesis for the orientation of an interest point.  I
    /// considered making a vector of orientations for a single point.
    /// However, it is probably better to make more than one interest
    /// point with the same (x,y,s) since the descriptor will be unique
    /// for a given orientation anyway.
    float orientation;

    /// The interest measure (could be Harris, LoG, etc.).
    float interest;

    /// These are some extras for SURF-like implementations
    bool polarity;
    /// This is the integer location in scale space (used for indexing
    /// a vector of interest images)
    uint32 octave, scale_lvl;

    /// And finally the descriptor for the interest point.  For example,
    /// PCA descriptors would have a vector of floats or doubles...
    descriptor_type descriptor;

    const_iterator begin() const { return descriptor.begin(); }
          iterator begin()       { return descriptor.begin(); }
    const_iterator end()   const { return descriptor.end(); }
          iterator end()         { return descriptor.end(); }

    size_t size() const { return descriptor.size(); }
    float operator[] (int index) { return descriptor[index]; }

    /// std::sort can be used to sort InterestPoints in descending
    /// order of interest.
    bool operator< (const InterestPoint& other) const {
      return (other.interest < interest);
    }

    /// Generates a human readable string
    std::string to_string() const {
      std::stringstream s;
      s << "IP: ("<<x<<","<<y<<")  scale: "<<scale<<"  orientation: "<<orientation<<"  interest: "<<interest<<
           "  polarity: "<<polarity<<"  octave: "<<octave<<"  scale_lvl: "<<scale_lvl<<"\n";
      s << "  descriptor: ";
      for (size_t i=0; i<descriptor.size(); ++i)
        s << descriptor[i] << "  ";
      s << std::endl;
      return s.str();
    }

#if defined(VW_HAVE_PKG_OPENCV) && VW_HAVE_PKG_OPENCV == 1
    // TODO: Move the definitions to the cc file!

    /// Copy IP information from an OpenCV KeyPoint object.
    void setFromCvKeypoint(cv::KeyPoint const& cvKey) {
      x  = cvKey.pt.x;
      y  = cvKey.pt.y;
      ix = round(x);
      iy = round(y);
      interest    = cvKey.response;
      octave      = cvKey.octave;
      scale_lvl   = cvKey.octave;
      scale       = cvKey.size;
      orientation = cvKey.angle;
      polarity    = false;
    }

    /// Create an OpenCV KeyPoint object from this IP.
    cv::KeyPoint makeOpenCvKeypoint() const {
      cv::KeyPoint cvKey;
      cvKey.pt.x     = x;
      cvKey.pt.y     = y;
      cvKey.response = interest;
      cvKey.octave   = octave;
      cvKey.size     = scale;
      cvKey.angle    = orientation;
      return cvKey;
    }

#endif
  }; // End class InterestPoint

inline bool InterestPointLessThan (InterestPoint P1, InterestPoint P2){
    if (P1.x           < P2.x           ) return true; if (P1.x           > P2.x           ) return false;
    if (P1.y           < P2.y           ) return true; if (P1.y           > P2.y           ) return false;
    if (P1.scale       < P2.scale       ) return true; if (P1.scale       > P2.scale       ) return false;
    if (P1.orientation < P2.orientation ) return true; if (P1.orientation > P2.orientation ) return false;
    if (P1.interest    < P2.interest    ) return true; if (P1.interest    > P2.interest    ) return false;
    if (P1.polarity    < P2.polarity    ) return true; if (P1.polarity    > P2.polarity    ) return false;
    if (P1.octave      < P2.octave      ) return true; if (P1.octave      > P2.octave      ) return false;
    if (P1.scale_lvl   < P2.scale_lvl   ) return true; if (P1.scale_lvl   > P2.scale_lvl   ) return false;
    return false;
}
  
  // Need to use list instead of vector for efficient thresholding.
  typedef std::list<InterestPoint> InterestPointList;

  // Utility function converts from a list of interest points to a
  // vector of interest point locations.  (Useful when preparing data for RANSAC.)
  std::vector<Vector3> iplist_to_vectorlist(std::vector<InterestPoint> const& iplist);
  //std::vector<InterestPoint> vectorlist_to_iplist(std::vector<Vector3> const& veclist); // Avoid using, info is lost!

  /// Select only the interest points that fall within the specified bounding box.
  template <class RealT>
  InterestPointList crop(InterestPointList const& interest_points, BBox<RealT,2> const& bbox) {
    InterestPointList return_val;
    for (InterestPointList::const_iterator i = interest_points.begin(); i != interest_points.end(); ++i) {
      if (bbox.contains(Vector<RealT,2>(RealT((*i).x), RealT((*i).y))))
        return_val.push_back(*i);
    }
    return return_val;
  }

  /// Helpful functors
  void remove_descriptor( InterestPoint & ip );

  
  /// Convert a an InterestPointList into a dense matrix of IP descriptors, one per row.
  template <class LIST_T, typename T>
  void ip_list_to_matrix(LIST_T const& ip_list, math::Matrix<T> &ip_matrix) {

    // Construct a plain matrix
    const size_t num_points        = ip_list.size();
    const size_t descriptor_length = ip_list.begin()->size();
    ip_matrix.set_size( num_points, descriptor_length );

    typename LIST_T::const_iterator ip_iter;
    size_t matrix_row = 0;
    for (ip_iter=ip_list.begin(); ip_iter != ip_list.end(); ++ip_iter) {
      std::copy( ip_iter->begin(), ip_iter->end(), ip_matrix[matrix_row].begin() );
      ++matrix_row;
    }
  }

  /// ImageInterestData
  ///
  /// This struct encapsulates some basic and widely useful processed
  /// views of a source image: the horizontal and vertical gradients,
  /// the orientation image, the gradient magnitude image, and the
  /// interest image. This is useful to ensure that these images are
  /// not redundantly calculated by different steps of the feature
  /// detection algorithm.
  ///
  /// The interest type is used to determine at compile-time which processed
  /// views should be fully rasterized. For speed in feature detection, the
  /// source type should be ImageView<T> or a simple manipulation of it.
  /// For memory efficiency, the source type should be ImageViewRef<T>.
  ///
  /// If some other sort of shared data is needed or any of the temporaries
  /// should be calculated in a different fashion, ImageInterestData can be
  /// partially specialized on InterestT.
  template <class SrcT, class InterestT>
  class ImageInterestData {
  public:
    /// The image types defined by InterestTraits control whether each processed
    /// view is fully rasterized or not. Only those used in calculating each
    /// pixel's interest measure should be fully rasterized. Later operations
    /// (thresholding, orientation assignment, etc.) require at most support
    /// regions around the interest points.
    typedef          SrcT             source_type;
    typedef typename SrcT::pixel_type pixel_type;

    typedef typename InterestOperatorTraits<SrcT, InterestT>::rasterize_type rasterize_type;
    typedef typename InterestOperatorTraits<SrcT, InterestT>::gradient_type  gradient_type;
    typedef typename InterestOperatorTraits<SrcT, InterestT>::mag_type       mag_type;
    typedef typename InterestOperatorTraits<SrcT, InterestT>::ori_type       ori_type;
    typedef typename InterestOperatorTraits<SrcT, InterestT>::interest_type  interest_type;
    typedef typename InterestOperatorTraits<SrcT, InterestT>::integral_type  integral_type;

    static const int peak_type = InterestPeakType<InterestT>::peak_type;

    /// Constructor which sets the source image and creates the processed views.
    /// This is a generic constructor that assumes the requires grad x/y.
    template <class ViewT>
    ImageInterestData(ImageViewBase<ViewT> const& img) :
      m_src(img.impl()),
      m_grad_x(derivative_filter(m_src, 1, 0).impl()),
      m_grad_y(derivative_filter(m_src, 0, 1).impl()),
      m_mag(hypot(m_grad_x, m_grad_y).impl()),
      m_ori(atan2(m_grad_y, m_grad_x).impl()),
      m_interest(NULL),
      m_integral(NULL) {}

    template <class ViewT>
    ImageInterestData(ImageViewBase<ViewT> const& img,
                      ImageViewBase<integral_type> const& integral ) :
      m_src(img.impl()),
      m_interest(NULL),
      m_integral(&integral.impl()) {}

    ~ImageInterestData() {
      if (m_interest) delete m_interest;
    }

    /// Accessors to immutable processed views.
    inline source_type   const& source     () const { return m_src;    }
    inline gradient_type const& gradient_x () const { return m_grad_x; }
    inline gradient_type const& gradient_y () const { return m_grad_y; }
    inline mag_type      const& magnitude  () const { return m_mag;    }
    inline ori_type      const& orientation() const { return m_ori;    }

    /// Accessors to mutable interest image.
    inline interest_type& interest() const {
      if (!m_interest) 
        vw_throw(LogicErr() << "ImageInterestData::interest() Interest image has not yet been computed.");
      return *m_interest;
    }

    template <class ViewT>
    inline void set_interest(ImageViewBase<ViewT> const& interest) {
      if (m_interest) 
        delete m_interest;
      m_interest = new interest_type(interest.impl());
    }

    /// Accessors to mutable integral image.
    inline integral_type const& integral() const {
      if (!m_integral) 
        vw_throw(LogicErr() << "ImageInterestData::integral() Integral image has not yet been computed.");
      return *m_integral;
    }

    template <class ViewT>
    inline void set_integral(ImageViewBase<ViewT> const& integral) {
      // I don't really recommend using this -ZMM
      m_integral = &integral;
    }

  protected:
    /// Cached processed data
    source_type   m_src;
    gradient_type m_grad_x, m_grad_y;
    mag_type      m_mag;
    ori_type      m_ori;
          interest_type *m_interest;
    const integral_type *m_integral;
  }; // End class ImageInterestData

}} // namespace vw::ip

#endif //__INTEREST_DATA_H__
