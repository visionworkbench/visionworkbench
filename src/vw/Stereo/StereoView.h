#ifndef __VW_STEREO_STEREOVIEW_H__
#define __VW_STEREO_STEREOVIEW_H__

#include <vw/Image/ImageViewBase.h>
#include <vw/Stereo/DisparityMap.h>
#include <vw/Stereo/StereoModel.h>
#include <vw/Camera/CameraModel.h>
#include <vw/Math/Vector.h>

namespace vw { 
namespace stereo {

  // Class definition
  template <class DisparityImageT>
  class StereoView : public ImageViewBase<StereoView<DisparityImageT> >
  {
  private:
    DisparityImageT m_disparity_map;
    
    StereoModel m_stereo_model;
    
  public:

    typedef Vector3 pixel_type;
    typedef const Vector3 result_type;
    typedef ProceduralPixelAccessor<StereoView> pixel_accessor;

    StereoView( DisparityImageT const& disparity_map, 
                vw::camera::CameraModel const& camera_model1, 
                vw::camera::CameraModel const& camera_model2) : 
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

    inline result_type operator()( int i, int j, int p=0 ) const { 
      double error;
      if ( !m_disparity_map(i,j).missing() ) {
        return m_stereo_model(Vector2( i, j ),
                              Vector2( i + m_disparity_map(i,j).h(),
                                       j + m_disparity_map(i,j).v() ), 
                              error );
      } else {
        // For missing pixels in the disparity map, we return a null 3D position.
        return Vector3();
      }
    }

    // Returns the distance between the two stereo rays at their
    // closest point of intersection.  A value of zero is returned if
    // the requested pixel is missing in the disparity map or if the
    // stereo geometry returns a bad match (e.g. if the rays
    // diverged).
    inline double error( int i, int j, int p=0 ) const { 
      double error;     
      if ( !disparity_map(i,j).missing() ) {
        m_disparity_map(Vector2( i, j ),
                        Vector2( i + disparity_map(i,j).h(),
                                 j + disparity_map(i,j).v() ), 
                           error );
        if (error >= 0) 
          return error;
      }

      // If we reach here, it means the disparity pixel or the 3D
      // point was invalid, and we return zero for error.
      return 0;
    }

    DisparityImageT const& disparity_map() const { return m_disparity_map; }

    /// \cond INTERNAL
    typedef StereoView<typename DisparityImageT::prerasterize_type> prerasterize_type;
    inline prerasterize_type prerasterize( BBox2i const& bbox ) const { return prerasterize_type( m_disparity_map.prerasterize(bbox), m_stereo_model ); }
    template <class DestT> inline void rasterize( DestT const& dest, BBox2i const& bbox ) const { vw::rasterize( prerasterize(bbox), dest, bbox ); }
    /// \endcond
  };

}} // namespace vw::stereo

#endif // __VW_STEREO_STEREOVIEW_H__
