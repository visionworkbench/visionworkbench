// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifndef __VW_STEREO_STEREOVIEW_H__
#define __VW_STEREO_STEREOVIEW_H__

#include <vw/Image/ImageViewBase.h>
#include <vw/Stereo/StereoModel.h>
#include <vw/Camera/CameraModel.h>
#include <limits>

namespace vw {

  // Registering point cloud for FileIO
  template<> struct PixelFormatID<Vector3> { static const PixelFormatEnum value = VW_PIXEL_GENERIC_3_CHANNEL; };
  template<> struct PixelFormatID<Vector3f> { static const PixelFormatEnum value = VW_PIXEL_GENERIC_3_CHANNEL; };

namespace stereo {

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

    template <class T>
    inline typename boost::enable_if<IsScalar<T>,Vector3>::type
    StereoModelHelper( StereoModel const& model, Vector2 const& index,
                       T const& disparity, double& error ) const {
      return model( index, Vector2( index[0] + disparity, index[1] ), error );
    }

    template <class T>
    inline typename boost::enable_if_c<IsCompound<T>::value && (CompoundNumChannels<typename UnmaskedPixelType<T>::type>::value == 1),Vector3>::type
    StereoModelHelper( StereoModel const& model, Vector2 const& index,
                       T const& disparity, double& error ) const {
      return model( index, Vector2( index[0] + disparity, index[1] ), error );
    }

    template <class T>
    inline typename boost::enable_if_c<IsCompound<T>::value && (CompoundNumChannels<typename UnmaskedPixelType<T>::type>::value != 1),Vector3>::type
    StereoModelHelper( StereoModel const& model, Vector2 const& index,
                       T const& disparity, double& error ) const {
      return model( index, Vector2( index[0] + disparity[0],
                                    index[1] + disparity[1] ), error );
    }

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
        return StereoModelHelper( m_stereo_model, Vector2(i,j),
                                  m_disparity_map(i,j,p), error );
      // For missing pixels in the disparity map, we return a null 3D position.
      return Vector3();
    }

    // Returns the distance between the two stereo rays at their
    // closest point of intersection.  A value of zero is returned if
    // the requested pixel is missing in the disparity map or if the
    // stereo geometry returns a bad match (e.g. if the rays
    // diverged).
    inline double error( int i, int j, int p=0 ) const {
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
      int rejected_points, total_points;
    };

    Vector3 m_origin;
    double m_near_radius, m_far_radius;
    boost::shared_ptr<UniverseRadiusState> m_state;

  public:
    UniverseRadiusFunc(Vector3 universe_origin, double near_radius = 0, double far_radius = std::numeric_limits<double>::max()):
      m_origin(universe_origin), m_near_radius(near_radius), m_far_radius(far_radius),
      m_state( new UniverseRadiusState() ) {
      VW_ASSERT(m_near_radius >= 0 && m_far_radius >= 0,
                vw::ArgumentErr() << "UniverseRadius: radii must be >= to zero.");
      VW_ASSERT(m_near_radius <= m_far_radius,
                vw::ArgumentErr() << "UniverseRadius: near radius must be <= to far radius.");
      m_state->rejected_points = m_state->total_points = 0;
    }

    double near_radius() const { return m_near_radius; }
    double far_radius() const { return m_far_radius; }
    int rejected_points() const { return m_state->rejected_points; }
    int total_points() const { return m_state->total_points; }

    Vector3 operator() (Vector3 const& pix) const {
      m_state->total_points++;
      if (pix != Vector3() ) {
        double dist = norm_2(pix - m_origin);
        if ((m_near_radius != 0 && dist < m_near_radius) ||
            (m_far_radius != 0 && dist > m_far_radius)) {
          m_state->rejected_points++;
          return Vector3();
        } else {
          return pix;
        }
      }
      return Vector3();
    }
  };

  // Useful routine for printing how many points have been rejected
  // using a particular UniverseRadius funtor.
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
