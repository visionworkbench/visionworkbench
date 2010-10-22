// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// \file InterestData.h
///
/// Basic classes and structures for storing image interest points.
///
#ifndef __INTEREST_DATA_H__
#define __INTEREST_DATA_H__

#include <vw/Math/Vector.h>
#include <vw/Math/Functors.h>
#include <vw/Image/ImageViewBase.h>
#include <vw/InterestPoint/InterestTraits.h>

#include <vector>
#include <list>
#include <algorithm>

namespace vw {
namespace ip {

  /// A class for storing information about an interest point.
  struct InterestPoint {
    typedef vw::Vector<float> descriptor_type;
    typedef descriptor_type::iterator iterator;
    typedef descriptor_type::const_iterator const_iterator;

    InterestPoint() {}

    InterestPoint(float x, float y, float scale=1.0, float interest=0.0, float ori=0.0,
		  bool pol=false, uint32 octave = 0, uint32 scale_lvl = 0)
      : x(x), y(y), scale(scale), ix(int32(x)), iy(int32(y)), orientation(ori), interest(interest),
	polarity(pol), octave(octave), scale_lvl(scale_lvl) {}


    /// Subpixel (col,row) location of point
    float x,y;

    /// Scale of point.  This may come from the pyramid level, from
    /// interpolating the interest function between levels, or from some
    /// other scale detector like the Laplace scale used by Mikolajczyk
    /// & Schmid
    float scale;

    /// Integer location (unnormalized), mainly for internal use.
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
    iterator begin() { return descriptor.begin(); }
    const_iterator end() const { return descriptor.end(); }
    iterator end() { return descriptor.end(); }

    size_t size() const { return descriptor.size(); }
    float operator[] (int index) { return descriptor[index]; }

    /// std::sort can be used to sort InterestPoints in descending
    /// order of interest.
    bool operator< (const InterestPoint& other) const {
      return (other.interest < interest);
    }

  };

  // Need to use list instead of vector for efficient thresholding.
  typedef std::list<InterestPoint> InterestPointList;

  // Utility function converts from a list of interest points to a
  // vector of interest point locations.  (Useful when preping data
  // far RANSAC...)
  std::vector<Vector3> iplist_to_vectorlist(std::vector<InterestPoint> const& iplist);
  std::vector<InterestPoint> vectorlist_to_iplist(std::vector<Vector3> const& veclist);

  // Routines for reading & writing interest point data files
  void write_lowe_ascii_ip_file(std::string ip_file, InterestPointList ip);
  void write_binary_ip_file(std::string ip_file, InterestPointList ip);
  std::vector<InterestPoint> read_binary_ip_file(std::string ip_file);

  // Routines for reading & writing interest point match files
  void write_binary_match_file(std::string match_file, std::vector<InterestPoint> const& ip1,
                               std::vector<InterestPoint> const& ip2);
  void read_binary_match_file(std::string match_file, std::vector<InterestPoint> &ip1,
                              std::vector<InterestPoint> &ip2);

  /// Select only the interest points that fall within the specified
  /// bounding box.
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
    typedef SrcT source_type;
    typedef typename SrcT::pixel_type pixel_type;

    typedef typename InterestOperatorTraits<SrcT, InterestT>::rasterize_type rasterize_type;
    typedef typename InterestOperatorTraits<SrcT, InterestT>::gradient_type gradient_type;
    typedef typename InterestOperatorTraits<SrcT, InterestT>::mag_type mag_type;
    typedef typename InterestOperatorTraits<SrcT, InterestT>::ori_type ori_type;
    typedef typename InterestOperatorTraits<SrcT, InterestT>::interest_type interest_type;
    typedef typename InterestOperatorTraits<SrcT, InterestT>::integral_type integral_type;

    static const int peak_type = InterestPeakType<InterestT>::peak_type;

    /// Constructor which sets the source image and creates the processed views.
    /// This is a generic constructor that assumes the requires grad
    /// x/y.
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
    inline source_type const& source() const { return m_src; }
    inline gradient_type const& gradient_x() const { return m_grad_x; }
    inline gradient_type const& gradient_y() const { return m_grad_y; }
    inline mag_type const& magnitude() const { return m_mag; }
    inline ori_type const& orientation() const { return m_ori; }

    /// Accessors to mutable interest image.
    inline interest_type& interest() const {
      if (!m_interest) vw_throw(LogicErr() << "ImageInterestData::interest() Interest image has not yet been computed.");
      return *m_interest;
    }

    template <class ViewT>
    inline void set_interest(ImageViewBase<ViewT> const& interest) {
      if (m_interest) delete m_interest;
      m_interest = new interest_type(interest.impl());
    }

    /// Accessors to mutable integral image.
    inline integral_type const& integral() const {
      if (!m_integral) vw_throw(LogicErr() << "ImageInterestData::integral() Integral image  has not yet been computed.");
      return *m_integral;
    }

    template <class ViewT>
    inline void set_integral(ImageViewBase<ViewT> const& integral) {
      // I don't really recommend using this -ZMM
      m_integral = &integral;
    }

  protected:
    /// Cached processed data
    source_type m_src;
    gradient_type m_grad_x, m_grad_y;
    mag_type m_mag;
    ori_type m_ori;
    interest_type *m_interest;
    const integral_type *m_integral;
  };


}} // namespace vw::ip

#endif //__INTEREST_DATA_H__
