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


#ifndef __VW_STEREO_STEREOVIEW_H__
#define __VW_STEREO_STEREOVIEW_H__

#include <vw/Image/ImageViewBase.h>
#include <vw/Image/PixelMask.h>
#include <vw/Image/PixelTypes.h>
#include <vw/Image/PixelAccessors.h>
#include <vw/Stereo/StereoModel.h>
#include <limits>

namespace vw {

  namespace camera {
    class CameraModel;
  }

namespace stereo {

  template <class T>
  inline typename boost::enable_if<IsScalar<T>,Vector2>::type
  DispHelper(T const& disparity) {
    return Vector2((double)disparity, 0);
  }
  
  template <class T>
  inline typename boost::enable_if_c<IsCompound<T>::value && (CompoundNumChannels<typename UnmaskedPixelType<T>::type>::value == 1),Vector2>::type
  DispHelper(T const& disparity) {
    return Vector2((double)disparity, 0);
  }

  template <class T>
  inline typename boost::enable_if_c<IsCompound<T>::value && (CompoundNumChannels<typename UnmaskedPixelType<T>::type>::value != 1),Vector2>::type
  DispHelper(T const& disparity) {
    return Vector2((double)disparity[0], (double)disparity[1]);
  }

  // Class definition
  template <class DisparityImageT>
  class StereoView : public ImageViewBase<StereoView<DisparityImageT> >
  {
    DisparityImageT m_disparity_map;
    StereoModel m_stereo_model;
    typedef typename DisparityImageT::pixel_type dpixel_type;

    template <class PixelT>
    struct NotSingleChannel {
      static const bool value = (1 != CompoundNumChannels<typename UnmaskedPixelType<PixelT>::type>::value);
    };

  public:

    typedef Vector3 pixel_type;
    typedef const Vector3 result_type;
    typedef ProceduralPixelAccessor<StereoView> pixel_accessor;

    StereoView( DisparityImageT const& disparity_map,
                vw::camera::CameraModel const* camera_model1,
                vw::camera::CameraModel const* camera_model2) :
      m_disparity_map(disparity_map),
      m_stereo_model(camera_model1, camera_model2) {}

    StereoView( DisparityImageT const& disparity_map,
                StereoModel const& stereo_model) :
      m_disparity_map(disparity_map),
      m_stereo_model(stereo_model) {}

    inline int32 cols() const { return m_disparity_map.cols(); }
    inline int32 rows() const { return m_disparity_map.rows(); }
    inline int32 planes() const { return 1; }

    inline pixel_accessor origin() const { return pixel_accessor(*this); }

    inline result_type operator()( size_t i, size_t j, size_t p=0 ) const {
      double error;
      if ( is_valid(m_disparity_map(i,j,p)) )
        return m_stereo_model(Vector2(i,j),
                              Vector2(i,j) + DispHelper(m_disparity_map(i,j,p)),
                              error);
      // For missing pixels in the disparity map, we return a null 3D position.
      return Vector3();
    }

    // Returns the distance between the two stereo rays at their
    // closest point of intersection.  A value of zero is returned if
    // the requested pixel is missing in the disparity map or if the
    // stereo geometry returns a bad match (e.g. if the rays
    // diverged).
    inline double error( int32 i, int32 j, int32 p=0 ) const {
      double error = 1e-10;
      if ( is_valid(m_disparity_map(i,j,p)) )
        StereoModelHelper( m_stereo_model, Vector2(i,j),
                           m_disparity_map(i,j,p), error );
      if ( error < 0 )
        return 0;
      return error;
    }

    DisparityImageT const& disparity_map() const { return m_disparity_map; }

    /// \cond INTERNAL
    typedef StereoView<typename DisparityImageT::prerasterize_type> prerasterize_type;
    inline prerasterize_type prerasterize( BBox2i const& bbox ) const { return prerasterize_type( m_disparity_map.prerasterize(bbox), m_stereo_model ); }
    template <class DestT> inline void rasterize( DestT const& dest, BBox2i const& bbox ) const { vw::rasterize( prerasterize(bbox), dest, bbox ); }
    /// \endcond
  };

  template <class ImageT>
  StereoView<ImageT> stereo_triangulate( ImageViewBase<ImageT> const& v,
                                         vw::camera::CameraModel const* camera1,
                                         vw::camera::CameraModel const* camera2 ) {
    return StereoView<ImageT>( v.impl(), camera1, camera2 );
  }

  // This per pixel functor applies a universe radius to a point
  // image.  Points that fall outside of the annulus specified with
  // the near and far radii are set to the "missing" pixel value
  // of (0,0,0).
  //
  // You can disable either the near or far universe radius
  // computation by setting either near_radius or far_radius to 0.0.
  class UniverseRadiusFunc: public vw::UnaryReturnSameType {
  private:

    // This small subclass gives us the wiggle room we need to update
    // the state of this object from within the PerPixelView, where
    // the UniverseRadius object itself is treated as a const copy of
    // the original functor. By maintaining a smart pointer to this
    // small status class, we can change state that is shared between
    // any copies of the UniverseRadius object and the original.
    struct UniverseRadiusState {
      int32 rejected_points, total_points;
    };

    Vector3 m_origin;
    double m_near_radius, m_far_radius;
    boost::shared_ptr<UniverseRadiusState> m_state;

  public:
    UniverseRadiusFunc(Vector3 universe_origin, double near_radius = 0, double far_radius = std::numeric_limits<double>::max()):
      m_origin(universe_origin), m_near_radius(near_radius), m_far_radius(far_radius),
      m_state( new UniverseRadiusState() ) {
      VW_ASSERT(m_near_radius >= 0 && m_far_radius >= 0,
                vw::ArgumentErr() << "UniverseRadius: radii must be >= 0.");
      VW_ASSERT(m_near_radius <= m_far_radius,
                vw::ArgumentErr() << "UniverseRadius: near radius must be <= far radius.");
      m_state->rejected_points = m_state->total_points = 0;
    }

    double near_radius() const { return m_near_radius; }
    double far_radius() const { return m_far_radius; }
    int32 rejected_points() const { return m_state->rejected_points; }
    int32 total_points() const { return m_state->total_points; }

    template <class ElemT>
    Vector<ElemT,3> operator() (Vector<ElemT,3> const& pix) const {
      m_state->total_points++;
      if (pix != Vector<ElemT,3>() ) {
        double dist = norm_2(pix - m_origin);
        if ((m_near_radius != 0 && dist < m_near_radius) ||
            (m_far_radius != 0 && dist > m_far_radius)) {
          m_state->rejected_points++;
          return Vector<ElemT,3>();
        } else {
          return pix;
        }
      }
      return Vector<ElemT,3>();
    }

    // A version that is carrying error
    template <class ElemT>
    Vector<ElemT,4> operator() (Vector<ElemT,4> const& pix) const {
      m_state->total_points++;
      if (subvector(pix,0,3) != Vector<ElemT,3>() ) {
        double dist = norm_2(subvector(pix,0,3) - m_origin);
        if ((m_near_radius != 0 && dist < m_near_radius) ||
            (m_far_radius != 0 && dist > m_far_radius)) {
          m_state->rejected_points++;
          return Vector<ElemT,4>();
        } else {
          return pix;
        }
      }
      return Vector<ElemT,4>();
    }

    // A version that is carrying 3D error
    template <class ElemT>
    Vector<ElemT,6> operator() (Vector<ElemT,6> const& pix) const {
      m_state->total_points++;
      if (subvector(pix,0,3) != Vector<ElemT,3>() ) {
        double dist = norm_2(subvector(pix,0,3) - m_origin);
        if ((m_near_radius != 0 && dist < m_near_radius) ||
            (m_far_radius != 0 && dist > m_far_radius)) {
          m_state->rejected_points++;
          return Vector<ElemT,6>();
        } else {
          return pix;
        }
      }
      return Vector<ElemT,6>();
    }
    
  };

  // Useful routine for printing how many points have been rejected
  // using a particular UniverseRadius functor.
  inline std::ostream& operator<<(std::ostream& os, UniverseRadiusFunc const& u) {
    os << "Universe Radius Limits: [ " << u.near_radius() << ", ";
    if (u.far_radius() == 0)
      os << "Inf ]\n" ;
    else
      os << u.far_radius() << "]\n";
    os << "\tRejected " << u.rejected_points() << "/" << u.total_points() << " vertices (" << double(u.rejected_points())/u.total_points()*100 << "%).\n";
    return os;
  }


}} // namespace vw::stereo

#endif // __VW_STEREO_STEREOVIEW_H__
