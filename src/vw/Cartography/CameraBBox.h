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

  // Intersect the ray back-projected from the camera with the datum.
  Vector3 datum_intersection( Datum const& datum,
                              camera::CameraModel const* model,
                              Vector2 const& pix );

  Vector3 datum_intersection( Datum const& datum,
                              Vector3 camera_ctr, Vector3 camera_vec );

  // Return the intersection between the ray emanating from the
  // current camera pixel with the datum ellipsoid. The return value
  // is a map-projected point location (the intermediate between
  // lon-lat-altitude and pixel).
  Vector2 geospatial_intersect( GeoReference const& georef,
                                Vector3 const& camera_ctr, Vector3 const& camera_vec,
                                bool& has_intersection );

  // Define an LMA model to solve for a DEM intersecting a ray.
  template <class DEMImageT>
  class DEMIntersectionLMA : public math::LeastSquaresModelBase< DEMIntersectionLMA< DEMImageT > > {
    InterpolationView<EdgeExtensionView<DEMImageT, ConstantEdgeExtension>, BilinearInterpolation> m_dem;
    GeoReference m_georef;
    Vector3 m_camera_ctr;
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
      // Note: We ignore the pixel mask, in effect treating invalid
      // pixels as having 0 value.
      return m_dem(x,y)[0];
    }

  public:
    // What is returned by evaluating the functor. In this case it is
    // the unit vector from the camera to the current xyz.
    typedef Vector3 result_type;
    // Defines the search space. In this case it is the point location
    // on the DEM.
    typedef Vector2 domain_type;
    // Jacobian form. Auto.
    typedef Matrix<double> jacobian_type;

    // Constructor
    DEMIntersectionLMA( ImageViewBase<DEMImageT> const& dem_image,
                        GeoReference const& georef,
                        Vector3 const& camera_ctr): m_dem(interpolate(dem_image)), m_georef(georef),
                                                    m_camera_ctr(camera_ctr){
      m_dem_bbox = bounding_box( m_dem );
    }

    inline Vector3 point_to_xyz_helper( domain_type const& projected_point ) const {
      // Convert a point in the projected space to a Cartesian point
      // by interpolating into the DEM.
      Vector2 dem_pixel = m_georef.point_to_pixel( projected_point );
      if ( !m_dem_bbox.contains( dem_pixel ) )
        return Vector3();
      Vector2 dem_lonlat = m_georef.point_to_lonlat( projected_point );
      Vector3 dem_xyz
        = m_georef.datum().geodetic_to_cartesian( Vector3( dem_lonlat.x(),
                                                           dem_lonlat.y(),
                                                           Helper<typename DEMImageT::pixel_type >(dem_pixel.x(),dem_pixel.y()))
                                                  );
      return dem_xyz;
    }

    // Evaluator
    inline result_type operator()( domain_type const& projected_point ) const {
      Vector3 xyz = point_to_xyz_helper(projected_point);
      if (xyz == Vector3()) return Vector3(-100,-100,-100);
      return normalize( xyz - m_camera_ctr );
    }
  };

  // Intersect the ray going from the given camera pixel with the DEM.
  // The return value is a point in the projected space. In this
  // function we treat no-data DEM values as 0.
  template <class DEMImageT>
  Vector2 camera_pixel_to_dem_point(Vector2 const& camera_pixel,
                                    ImageViewBase<DEMImageT> const& dem_image,
                                    GeoReference const& georef,
                                    boost::shared_ptr<camera::CameraModel> camera_model,
                                    bool & has_intersection
                                    ){
    has_intersection = true;

    double max_abs_tol = 1e-14, max_rel_tol = 1e-14;
    int num_max_iter = 100;

    Vector3 camera_ctr = camera_model->camera_center(camera_pixel);
    Vector3 camera_vec = camera_model->pixel_to_vector(camera_pixel);

    // Get an initial guess from intersecting the ray with the datum.
    Vector2 projected_point = geospatial_intersect(georef, camera_ctr, camera_vec, has_intersection);
    if ( !has_intersection ) {
      return Vector2();
    }

    // Refining the intersection using Levenberg-Marquardt
    DEMIntersectionLMA<DEMImageT> model(dem_image, georef, camera_ctr);
    int status = 0;
    projected_point = math::levenberg_marquardt(model, projected_point,
                                                camera_vec, status,
                                                max_abs_tol, max_rel_tol,
                                                num_max_iter
                                                );

    if ( status < 0 ) {
      has_intersection = false;
      return Vector2();
    }

    has_intersection = true;
    return projected_point;
  }

  // Define an LMA model to solve for a DEM intersecting a ray. The
  // variable of optimization is position on the ray. The cost
  // function is difference between datum height and DEM height at
  // current point on the ray.
  template <class DEMImageT>
  class RayDEMIntersectionLMA : public math::LeastSquaresModelBase< RayDEMIntersectionLMA< DEMImageT > > {
    InterpolationView<EdgeExtensionView<DEMImageT, ConstantEdgeExtension>, BilinearInterpolation> m_dem;
    GeoReference m_georef;
    Vector3 m_camera_ctr;
    Vector3 m_camera_vec;
    bool m_treat_nodata_as_zero;
    // Provide safe interaction with DEMs that are scalar or compound
    template <class PixelT>
    typename boost::enable_if< IsScalar<PixelT>, double >::type
    inline Helper( double x, double y ) const {
      if ( 0 <= x && x <= m_dem.cols() - 1 && // for interpolation
           0 <= y && y <= m_dem.rows() - 1 ){
        PixelT val = m_dem(x, y);
        if (is_valid(val)) return val;
      }
      if (m_treat_nodata_as_zero) return 0;
      return big_val();
    }

    template <class PixelT>
    typename boost::enable_if< IsCompound<PixelT>, double>::type
    inline Helper( double x, double y ) const {
      if ( 0 <= x && x <= m_dem.cols() - 1 && // for interpolation
           0 <= y && y <= m_dem.rows() - 1 ){
        PixelT val = m_dem(x, y);
        if (is_valid(val)) return val[0];
      }
      if (m_treat_nodata_as_zero) return 0;
      return big_val();
    }

  public:
    typedef Vector<double, 1> result_type;
    typedef Vector<double, 1> domain_type;
    // Jacobian form. Auto.
    typedef Matrix<double> jacobian_type;

    inline double big_val() const {
      // Don't make this too big as in the LMA algorithm it may get
      // squared and may cause overflow.
      return 1.0e+50;
    }

    // Constructor
    RayDEMIntersectionLMA(ImageViewBase<DEMImageT> const& dem_image,
                          GeoReference const& georef,
                          Vector3 const& camera_ctr,
                          Vector3 const& camera_vec,
                          bool treat_nodata_as_zero
                          )
      : m_dem(interpolate(dem_image)), m_georef(georef),
        m_camera_ctr(camera_ctr), m_camera_vec(camera_vec),
        m_treat_nodata_as_zero(treat_nodata_as_zero){}

    // Evaluator. See description above.
    inline result_type operator()( domain_type const& len ) const {
      Vector3 xyz = m_camera_ctr + len[0]*m_camera_vec;
      Vector3 llh = m_georef.datum().cartesian_to_geodetic( xyz );
      Vector2 pix = m_georef.lonlat_to_pixel( Vector2( llh.x(), llh.y() ) );
      result_type result;
      result[0] = Helper<typename DEMImageT::pixel_type >(pix.x(),pix.y()) - llh[2];
      return result;
    }
  };

  // Intersect the ray going from the given camera pixel with a DEM.
  // The return value is a Cartesian point. If the ray goes through a
  // hole in the DEM where there is no data, we return no-intersection
  // or intersection with the datum, depending on whether the variable
  // treat_nodata_as_zero is false or true.
  template <class DEMImageT>
  Vector3 camera_pixel_to_dem_xyz(Vector3 const& camera_ctr, Vector3 const& camera_vec,
                                  ImageViewBase<DEMImageT> const& dem_image,
                                  GeoReference const& georef,
                                  bool treat_nodata_as_zero,
                                  bool & has_intersection,
                                  double height_error_tol = 1e-1,  // error in DEM height
                                  double max_abs_tol      = 1e-14, // abs cost function change b/w iters
                                  double max_rel_tol      = 1e-14,
                                  int num_max_iter        = 100,
                                  Vector3 xyz_guess       = Vector3()
                                  ){

    has_intersection = false;
    RayDEMIntersectionLMA<DEMImageT> model(dem_image, georef,
                                           camera_ctr, camera_vec, treat_nodata_as_zero);

    Vector3 xyz;
    if ( xyz_guess == Vector3() ){
      // Intersect the ray with the datum, this is a good initial
      // guess.
      xyz = datum_intersection(georef.datum(), camera_ctr, camera_vec);

      if ( xyz == Vector3() ) {
        has_intersection = false;
        return Vector3();
      }
    }else{
      xyz = xyz_guess;
    }

    // Length along the ray from camera center to intersection point
    Vector<double, 1> len0, len;
    len0[0] = norm_2(xyz - camera_ctr);

    // If the ray intersects the datum at a point which does not
    // correspond to a valid location in the DEM, wiggle that point
    // along the ray until hopefully it does.
    double radius = norm_2(xyz);
    int n = 10;
    double small = radius*0.02/( 1 << (n-1) ); // wiggle up to 0.02*radius
    for (int i = 0; i <= n; i++){
      double delta = 0;
      if (i > 0) delta = small*( 1 << (i-1) );
      for (int k = -1; k <= 1; k += 2){
        len[0] = len0[0] + k*delta;
        Vector<double, 1> height_diff = model(len);
        if ( std::abs(height_diff[0]) < model.big_val()/10.0 ){
          has_intersection = true;
          break;
        }
      }
      if (has_intersection) break;
    }

    if ( !has_intersection ) {
      return Vector3();
    }

    // Refining the intersection using Levenberg-Marquardt
    int status = 0;
    Vector<double, 1> observation; observation[0] = 0;
    len = math::levenberg_marquardt(model, len, observation, status,
                                    max_abs_tol, max_rel_tol,
                                    num_max_iter
                                    );

    Vector<double, 1> dem_height = model(len);

    if ( status < 0 || std::abs(dem_height[0]) > height_error_tol ){
      has_intersection = false;
      return Vector3();
    }

    has_intersection = true;
    xyz = camera_ctr + len[0]*camera_vec;
    return xyz;
  }

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
      Vector2 m_last_intersect;

    public:
      bool last_valid, center_on_zero;
      BBox2 box;
      double scale;

      CameraDatumBBoxHelper( GeoReference const& georef,
                             boost::shared_ptr<camera::CameraModel> camera,
                             bool center=false) : m_georef(georef), m_camera(camera), last_valid(false), center_on_zero(center), scale( std::numeric_limits<double>::max() ) {}

      void operator() ( Vector2 const& pixel );
    };

    template <class DEMImageT>
    class CameraDEMBBoxHelper {
      GeoReference m_georef;
      boost::shared_ptr<camera::CameraModel> m_camera;
      DEMImageT m_dem;
      Vector2 m_last_intersect;

    public:
      bool last_valid, center_on_zero;
      BBox2 box;
      double scale;

      CameraDEMBBoxHelper( ImageViewBase<DEMImageT> const& dem_image,
                           GeoReference const& georef,
                           boost::shared_ptr<camera::CameraModel> camera,
                           bool center=false ) : m_georef(georef), m_camera(camera), m_dem(dem_image.impl()), last_valid(false), center_on_zero(center), scale( std::numeric_limits<double>::max() ) {}

      void operator() ( Vector2 const& pixel ) {

        bool has_intersection;
        Vector2 point
          = camera_pixel_to_dem_point(pixel, m_dem, m_georef, m_camera, has_intersection);

        if ( !has_intersection ) {
          last_valid = false;
          return;
        }

        if (!m_georef.is_projected()){
          // If we don't use a projected coordinate system, then the
          // coordinates of this point are simply lon and lat.
          if ( center_on_zero && point[0] > 180 )
            point[0] -= 360.0;
          else if ( center_on_zero && point[0] < -180 )
            point[0] += 360.0;
          else if ( !center_on_zero && point[0] < 0 )
            point[0] += 360.0;
          else if ( !center_on_zero && point[0] > 360 )
            point[0] -= 360.0;
        }

        if ( last_valid ) {
          double current_scale =
            norm_2( point - m_last_intersect );
          if ( current_scale < scale )
            scale = current_scale;
        }
        m_last_intersect = point;
        box.grow( point );
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

    // To do: Integrate the almost identical functions camera_bbox in
    // CameraBBox.h and CameraBBox.cc. One of them uses a DEM and the
    // second one does not.

    // Testing to see if we should be centering on zero
    bool center_on_zero = true;
    Vector3 camera_llr =
      XYZtoLonLatRadFunctor::apply(camera_model->camera_center(Vector2()));
    if ( camera_llr[0] < -90 ||
         camera_llr[0] > 90 )
      center_on_zero = false;

    int32 step_amount = (2*cols+2*rows)/100;
    step_amount = std::min(step_amount, cols/4); // must have at least several points per column
    step_amount = std::min(step_amount, rows/4); // must have at least several points per row
    step_amount = std::max(step_amount, 1);      // step amount must be > 0

    detail::CameraDEMBBoxHelper<DEMImageT> functor( dem_image, georef, camera_model,
                                                    center_on_zero );

    // Running the edges. Note: The last valid point on a
    // BresenhamLine is the last point before the endpoint.
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
