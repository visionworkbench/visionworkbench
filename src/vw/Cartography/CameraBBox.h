// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
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
#include <vw/Cartography/GeoReference.h>
#include <vw/Camera/CameraModel.h>
#include <vw/Math/LevenbergMarquardt.h>

#include <boost/shared_ptr.hpp>

namespace vw {
namespace cartography {

  // Intersection Devices
  //////////////////////////////////////////////////////

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

    double semi_major_axis = georef.datum().semi_major_axis();
    double semi_minor_axis = georef.datum().semi_minor_axis();
    double z_scale = semi_major_axis / semi_minor_axis;

    BBox2 georeference_space_bbox;
    bool last_valid=false;
    Vector2 last_geospatial_point;
    scale = -1;

    // Top row
    for( int x=0; x<cols; ++x ) {
      Vector2 pix(x,0);
      bool test_intersect;
      Vector2 geospatial_point = geospatial_intersect( pix, georef,
                                                       camera_model,
                                                       z_scale,
                                                       test_intersect );
      if ( !test_intersect ) {
        last_valid = false;
        continue;
      }

      // Refining with more accurate intersection
      DEMIntersectionLMA<DEMImageT> model( dem_image, georef, camera_model );
      int status = 0;
      geospatial_point = math::levenberg_marquardt( model, geospatial_point,
                                                    pix, status );
      if ( status < 0 ) {
        last_valid = false;
        continue;
      }

      if( last_valid ) {
        double current_scale = norm_2( geospatial_point - last_geospatial_point );
        if ( current_scale < 0 ||
             current_scale < scale )
          scale = current_scale;
      }
      last_geospatial_point = geospatial_point;
      georeference_space_bbox.grow( geospatial_point );
      last_valid = true;
    }
    // Bottom row
    for( int x=cols-1; x>=0; --x ) {
      Vector2 pix(x,rows-1);

      bool test_intersect;
      Vector2 geospatial_point = geospatial_intersect( pix, georef,
                                                       camera_model,
                                                       z_scale,
                                                       test_intersect );
      if ( !test_intersect ) {
        last_valid = false;
        continue;
      }

      // Refining with more accurate intersection
      DEMIntersectionLMA<DEMImageT> model( dem_image, georef, camera_model );
      int status = 0;
      geospatial_point = math::levenberg_marquardt( model, geospatial_point,
                                                    pix, status );
      if ( status < 0 ) {
        last_valid = false;
        continue;
      }

      if( last_valid ) {
        double current_scale = norm_2( geospatial_point - last_geospatial_point );
        if ( current_scale < 0 ||
             current_scale < scale )
          scale = current_scale;
      }
      last_geospatial_point = geospatial_point;
      georeference_space_bbox.grow( geospatial_point );
      last_valid = true;
    }
    // Left side
    for( int y=rows-2; y>0; --y ) {
      Vector2 pix(0,y);

      bool test_intersect;
      Vector2 geospatial_point = geospatial_intersect( pix, georef,
                                                       camera_model,
                                                       z_scale,
                                                       test_intersect );
      if ( !test_intersect ) {
        last_valid = false;
        continue;
      }

      // Refining with more accurate intersection
      DEMIntersectionLMA<DEMImageT> model( dem_image, georef, camera_model );
      int status = 0;
      geospatial_point = math::levenberg_marquardt( model, geospatial_point,
                                                    pix, status );
      if ( status < 0 ) {
        last_valid = false;
        continue;
      }

      if( last_valid ) {
        double current_scale = norm_2( geospatial_point - last_geospatial_point );
        if ( current_scale < 0 ||
             current_scale < scale )
          scale = current_scale;
      }
      last_geospatial_point = geospatial_point;
      georeference_space_bbox.grow( geospatial_point );
      last_valid = true;
    }
    // Right side
    last_valid = false;
    for( int y=1; y<rows-1; ++y ) {
      Vector2 pix(cols-1,y);

      bool test_intersect;
      Vector2 geospatial_point = geospatial_intersect( pix, georef,
                                                       camera_model,
                                                       z_scale,
                                                       test_intersect );
      if ( !test_intersect ) {
        last_valid = false;
        continue;
      }

      // Refining with more accurate intersection
      DEMIntersectionLMA<DEMImageT> model( dem_image, georef, camera_model );
      int status = 0;
      geospatial_point = math::levenberg_marquardt( model, geospatial_point,
                                                    pix, status );
      if ( status < 0 ) {
        last_valid = false;
        continue;
      }

      if( last_valid ) {
        double current_scale = norm_2( geospatial_point - last_geospatial_point );
        if ( scale < 0 ||
             current_scale < scale )
          scale = current_scale;
      }
      last_geospatial_point = geospatial_point;
      georeference_space_bbox.grow( geospatial_point );
      last_valid = true;
    }

    BBox2 bbox = georeference_space_bbox;

    // Did we fail to find scale?
    if ( scale == -1 ) {
      Vector2 pix(cols,rows);
      pix = pix / 2;
      Vector2 pix2 = pix + Vector2(1,1);

      bool test_intersect, test_intersect2;
      Vector2 geospatial_point = geospatial_intersect( pix, georef,
                                                       camera_model,
                                                       z_scale,
                                                       test_intersect );
      Vector2 geospatial_point2 = geospatial_intersect( pix2, georef,
                                                        camera_model,
                                                        z_scale,
                                                        test_intersect2 );
      if ( (!test_intersect) || (!test_intersect2) )
        return bbox;
      scale = norm_2( geospatial_point - geospatial_point2 );
    }

    return bbox;
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
