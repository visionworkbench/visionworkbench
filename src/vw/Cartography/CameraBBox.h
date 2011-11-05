// __BEGIN_LICENSE__
// Copyright (C) 2006-2011 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifndef __VW_CARTOGRAPHY_CAMERABBOX_H__
#define __VW_CARTOGRAPHY_CAMERABBOX_H__

#include <vw/config.h>
#if defined(VW_HAVE_PKG_CAMERA) && (VW_HAVE_PKG_CAMERA==1)

#include <vw/Image/ImageViewRef.h>
#include <vw/Image/Algorithms.h>
#include <vw/Image/Interpolation.h>
#include <vw/Math/LevenbergMarquardt.h>
#include <vw/Camera/CameraModel.h>
#include <vw/Cartography/SimplePointImageManipulation.h>
#include <vw/Cartography/GeoReference.h>
#include <vw/Cartography/detail/BresenhamLine.h>

#include <boost/shared_ptr.hpp>

namespace vw {
namespace cartography {

  // Return map projected point location (the intermediate between LLA
  // and Pixel)
  Vector2 geospatial_intersect( Vector2 pix,
                                GeoReference const& georef,
                                boost::shared_ptr<camera::CameraModel> camera_model,
                                double z_scale, bool& did_intersect );

  // Define an LMA model to solve for an intersection ...
  template <class DEMImageT>
  class DEMIntersectionLMA : public math::LeastSquaresModelBase< DEMIntersectionLMA< DEMImageT > > {
    InterpolationView<EdgeExtensionView<DEMImageT, ConstantEdgeExtension>, BilinearInterpolation> m_dem;
    GeoReference m_georef;
    boost::shared_ptr<camera::CameraModel> m_camera_model;
    BBox2i m_dem_bbox;

    // Provide safe interaction with DEMs that are scalar or compound
    template <class PixelT>
    typename boost::enable_if< IsScalar<PixelT>, double >::type
    inline Helper( double const& x, double const& y ) const {
      return m_dem(x,y);
    }

    template <class PixelT>
    typename boost::enable_if< IsCompound<PixelT>, double>::type
    inline Helper( double const& x, double const& y ) const {
      return m_dem(x,y)[0];
    }

  public:
    // What is returned by evaluating the functor. In this case it is
    // the projection into the camera.
    typedef Vector<double> result_type;
    // Defines the search space. In this case it is the point location
    // on the DEM.
    typedef Vector<double> domain_type;
    // Jacobian form. Auto.
    typedef Matrix<double> jacobian_type;

    // Constructor
    DEMIntersectionLMA( ImageViewBase<DEMImageT> const& dem_image,
                        GeoReference const& georef,
                        boost::shared_ptr<camera::CameraModel> camera_model ) :
    m_dem(interpolate(dem_image)), m_georef(georef), m_camera_model(camera_model) {
      m_dem_bbox = bounding_box( m_dem );
    }

    // Evaluator
    inline result_type operator()( domain_type const& x ) const {
      Vector2 dem_pixel_loc = m_georef.point_to_pixel( x );
      if ( !m_dem_bbox.contains( dem_pixel_loc ) )
        return Vector2(-100,-100);
      Vector2 dem_ll_loc = m_georef.point_to_lonlat( x );
      Vector3 dem_cart_loc = m_georef.datum().geodetic_to_cartesian( Vector3( dem_ll_loc.x(), dem_ll_loc.y(), Helper<typename DEMImageT::pixel_type >(dem_pixel_loc.x(),dem_pixel_loc.y())) );
      return m_camera_model->point_to_pixel(dem_cart_loc);
    }
  };

  namespace detail {

    template <class FunctionT>
    void bresenham_apply( BresenhamLine line, size_t step,
                          FunctionT& f ) {
      while ( line.is_good() ) {
        f( *line );
        for ( size_t i = 0; i < step; i++ )
          ++line;
      }
    }

    class CameraDatumBBoxHelper {
      GeoReference m_georef;
      boost::shared_ptr<camera::CameraModel> m_camera;
      double m_z_scale;
      Vector2 m_last_intersect;

    public:
      bool last_valid, center_on_zero;
      BBox2 box;
      double scale;

      CameraDatumBBoxHelper( GeoReference const& georef,
                             boost::shared_ptr<camera::CameraModel> camera,
                             bool center=false) : m_georef(georef), m_camera(camera), last_valid(false), center_on_zero(center), scale( std::numeric_limits<double>::max() ) {
        m_z_scale = m_georef.datum().semi_major_axis() / m_georef.datum().semi_minor_axis();
      }

      void operator() ( Vector2 pixel ) {
        bool test_intersect;
        Vector2 geospatial_point =
          geospatial_intersect( pixel, m_georef, m_camera,
                                m_z_scale, test_intersect );
        if ( !test_intersect ) {
          last_valid = false;
          return;
        }

        if ( center_on_zero && geospatial_point[0] > 180 )
          geospatial_point[0] -= 360.0;

        if ( last_valid ) {
          double current_scale =
            norm_2( geospatial_point - m_last_intersect );
          if ( current_scale < scale )
            scale = current_scale;
        }
        m_last_intersect = geospatial_point;
        box.grow( geospatial_point );
        last_valid = true;
      }
    };

    template <class DEMImageT>
    class CameraDEMBBoxHelper {
      GeoReference m_georef;
      boost::shared_ptr<camera::CameraModel> m_camera;
      DEMImageT m_dem;
      double m_z_scale;
      Vector2 m_last_intersect;

    public:
      bool last_valid, center_on_zero;
      BBox2 box;
      double scale;

      CameraDEMBBoxHelper( ImageViewBase<DEMImageT> const& dem_image,
                           GeoReference const& georef,
                           boost::shared_ptr<camera::CameraModel> camera,
                           bool center=false ) : m_georef(georef), m_camera(camera), m_dem(dem_image.impl()), last_valid(false), center_on_zero(center), scale( std::numeric_limits<double>::max() ) {
        m_z_scale = m_georef.datum().semi_major_axis() / m_georef.datum().semi_minor_axis();
      }

      void operator() ( Vector2 pixel ) {
        bool test_intersect;
        Vector2 geospatial_point =
          geospatial_intersect( pixel, m_georef, m_camera,
                                m_z_scale, test_intersect );
        if ( !test_intersect ) {
          last_valid = false;
          return;
        }

        // Refining with more accurate intersection
        DEMIntersectionLMA<DEMImageT> model( m_dem, m_georef, m_camera );
        int status = 0;
        geospatial_point = math::levenberg_marquardt( model, geospatial_point,
                                                      pixel, status );
        if ( status < 0 ) {
          last_valid = false;
          return;
        }

        if ( center_on_zero && geospatial_point[0] > 180 )
          geospatial_point[0] -= 360.0;
        else if ( center_on_zero && geospatial_point[0] < -180 )
          geospatial_point[0] += 360.0;
        else if ( !center_on_zero && geospatial_point[0] < 0 )
          geospatial_point[0] += 360.0;
        else if ( !center_on_zero && geospatial_point[0] > 360 )
          geospatial_point[0] -= 360.0;

        if ( last_valid ) {
          double current_scale =
            norm_2( geospatial_point - m_last_intersect );
          if ( current_scale < scale )
            scale = current_scale;
        }
        m_last_intersect = geospatial_point;
        box.grow( geospatial_point );
        last_valid = true;
      }
    };
  }

  // Functions for Users
  //////////////////////////////////////////////////////

  // Simple Intersection interfaces
  BBox2 camera_bbox( GeoReference const& georef,
                     boost::shared_ptr<vw::camera::CameraModel> camera_model,
                     int32 cols, int32 rows, float &scale );

  inline BBox2 camera_bbox( GeoReference const& georef,
                            boost::shared_ptr<vw::camera::CameraModel> camera_model,
                            int32 cols, int32 rows ) {
    float scale;
    return camera_bbox( georef, camera_model, cols, rows, scale );
  }

  // Intersections that take in account DEM topography
  template< class DEMImageT >
  BBox2 camera_bbox( ImageViewBase<DEMImageT> const& dem_image,
                     GeoReference const& georef,
                     boost::shared_ptr<vw::camera::CameraModel> camera_model,
                     int32 cols, int32 rows, float &scale ) {

    // Testing to see if we should be centering on zero
    bool center_on_zero = true;
    Vector3 camera_llr =
      XYZtoLonLatRadFunctor::apply(camera_model->camera_center(Vector2()));
    if ( camera_llr[0] < -90 ||
         camera_llr[0] > 90 )
      center_on_zero = false;

    int32 step_amount = (2*cols+2*rows)/100;
    detail::CameraDEMBBoxHelper<DEMImageT> functor( dem_image, georef, camera_model,
                                                    center_on_zero );

    // Running the edges
    bresenham_apply( BresenhamLine(0,0,cols,0),
                     step_amount, functor );
    functor.last_valid = false;
    bresenham_apply( BresenhamLine(cols-1,0,cols-1,rows),
                     step_amount, functor );
    functor.last_valid = false;
    bresenham_apply( BresenhamLine(cols-1,rows-1,0,rows-1),
                     step_amount, functor );
    functor.last_valid = false;
    bresenham_apply( BresenhamLine(0,rows-1,0,0),
                     step_amount, functor );
    functor.last_valid = false;

    // Running once through the center
    bresenham_apply( BresenhamLine(0,0,cols,rows),
                     step_amount, functor );

    scale = functor.scale/double(step_amount);
    return functor.box;
  }

  template< class DEMImageT >
  inline BBox2 camera_bbox( ImageViewBase<DEMImageT> const& dem_image,
                            GeoReference const& georef,
                            boost::shared_ptr<vw::camera::CameraModel> camera_model,
                            int32 cols, int32 rows ) {
    float scale;
    return camera_bbox<DEMImageT>( dem_image.impl(), georef, camera_model, cols, rows, scale );
  }

} // namespace cartography
} // namespace vw

#endif // VW_HAVE_PKG_CAMERA

#endif // __VW_CARTOGRAPHY_CAMERABBOX_H__
