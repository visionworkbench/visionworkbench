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

/// \file CameraBBox.h Contains bounding box, pixel interesection, and misc utilities.

#include <vw/config.h>
#if defined(VW_HAVE_PKG_CAMERA) && (VW_HAVE_PKG_CAMERA==1)

#include <vw/Image/ImageViewRef.h>
#include <vw/Image/Algorithms.h>
#include <vw/Image/Interpolation.h>
#include <vw/Math/LevenbergMarquardt.h>
#include <vw/Math/BresenhamLine.h>
#include <vw/Camera/CameraModel.h>
#include <vw/Cartography/PointImageManipulation.h>
#include <vw/Cartography/GeoReference.h>


#include <boost/shared_ptr.hpp>


namespace vw {
namespace cartography {

  // Intersect the ray back-projected from the camera with the datum.
  Vector3 datum_intersection( Datum const& datum,
                              camera::CameraModel const* model,
                              Vector2 const& pix );

  // Return the intersection between the ray emanating from the
  // current camera pixel with the datum ellipsoid. The return value
  // is a map-projected point location (the intermediate between
  // lon-lat-altitude and pixel).
  Vector2 geospatial_intersect( GeoReference const& georef,
                                Vector3 const& camera_ctr, Vector3 const& camera_vec,
                                bool& has_intersection );


  // Define an LMA model to solve for a DEM intersecting a ray. The
  // variable of optimization is position on the ray. The cost
  // function is difference between datum height and DEM height at
  // current point on the ray.
  template <class DEMImageT>
  class RayDEMIntersectionLMA : public math::LeastSquaresModelBase< RayDEMIntersectionLMA< DEMImageT > > {

    // TODO: Why does this use EdgeExtension if Helper() restricts access to the bounds?
    InterpolationView<EdgeExtensionView<DEMImageT, ConstantEdgeExtension>,
                      BilinearInterpolation> m_dem;
    GeoReference m_georef;
    Vector3      m_camera_ctr;
    Vector3      m_camera_vec;
    bool         m_treat_nodata_as_zero;

    /// Provide safe interaction with DEMs that are scalar
    /// - If m_dem(x,y) is in bounds, return the interpolated value.
    /// - Otherwise return 0 or big_val()
    template <class PixelT>
    typename boost::enable_if< IsScalar<PixelT>, double >::type
    inline Helper( double x, double y ) const {
      if ( (0 <= x) && (x <= m_dem.cols() - 1) && // for interpolation
           (0 <= y) && (y <= m_dem.rows() - 1) ){
        PixelT val = m_dem(x, y);
        if (is_valid(val)) return val;
      }
      if (m_treat_nodata_as_zero) return 0;
      return big_val();
    }

    /// Provide safe interaction with DEMs that are compound
    template <class PixelT>
    typename boost::enable_if< IsCompound<PixelT>, double>::type
    inline Helper( double x, double y ) const {
      if ( (0 <= x) && (x <= m_dem.cols() - 1) && // for interpolation
           (0 <= y) && (y <= m_dem.rows() - 1) ){
        PixelT val = m_dem(x, y);
        if (is_valid(val)) return val[0];
      }
      if (m_treat_nodata_as_zero) return 0;
      return big_val();
    }

  public:
    typedef Vector<double, 1> result_type;
    typedef Vector<double, 1> domain_type;
    typedef Matrix<double>    jacobian_type; ///< Jacobian form. Auto.

    /// Return a very large error to penalize locations that fall off the edge of the DEM.
    inline double big_val() const {
      // Don't make this too big as in the LMA algorithm it may get squared and may cause overflow.
      return 1.0e+50;
    }

    /// Constructor
    RayDEMIntersectionLMA(ImageViewBase<DEMImageT> const& dem_image,
                          GeoReference const& georef,
                          Vector3 const& camera_ctr,
                          Vector3 const& camera_vec,
                          bool treat_nodata_as_zero
                          )
      : m_dem(interpolate(dem_image)), m_georef(georef),
        m_camera_ctr(camera_ctr), m_camera_vec(camera_vec),
        m_treat_nodata_as_zero(treat_nodata_as_zero){}

    /// Evaluator. See description above.
    inline result_type operator()( domain_type const& len ) const {
      // The proposed intersection point
      Vector3 xyz = m_camera_ctr + len[0]*m_camera_vec;

      // Convert to geodetic coordinates, then to DEM pixel coordinates
      Vector3 llh = m_georef.datum().cartesian_to_geodetic( xyz );
      Vector2 pix = m_georef.lonlat_to_pixel( Vector2( llh.x(), llh.y() ) );

      // Return a measure of the elevation difference between the DEM and the guess
      // at its current location.
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
    RayDEMIntersectionLMA<DEMImageT> model(dem_image, georef, camera_ctr,
                                           camera_vec, treat_nodata_as_zero);

    Vector3 xyz;
    if ( xyz_guess == Vector3() ){ // If no guess provided
      // Intersect the ray with the datum, this is a good initial guess.
      xyz = datum_intersection(georef.datum(), camera_ctr, camera_vec);

      if ( xyz == Vector3() ) { // If we failed to intersect the datum, give up!
        has_intersection = false;
        return Vector3();
      }
    }else{ // User provided guess
      xyz = xyz_guess;
    }

    // Length along the ray from camera center to intersection point
    Vector<double, 1> len0, len;
    len0[0] = norm_2(xyz - camera_ctr);

    // If the ray intersects the datum at a point which does not
    // correspond to a valid location in the DEM, wiggle that point
    // along the ray until hopefully it does.
    const double radius     = norm_2(xyz); // Radius from XZY coordinate center
    const int    ITER_LIMIT = 10; // There are two solver attempts per iteration
    const double small      = radius*0.02/( 1 << (ITER_LIMIT-1) ); // Wiggle
    for (int i = 0; i <= ITER_LIMIT; i++){
      // Gradually expand delta until on final iteration it is == radius*0.02
      double delta = 0;
      if (i > 0)
        delta = small*( 1 << (i-1) );

      for (int k = -1; k <= 1; k += 2){ // For k==-1, k==1
        len[0] = len0[0] + k*delta; // Ray guess length +/- 2% planetary radius
        // Use our model to compute the height diff at this length
        Vector<double, 1> height_diff = model(len);
        // TODO: This is an EXTREMELY lenient threshold! big_val()/10.0 == 1.0e+49!!!
        //       The effect of this may be to just stop this loop when we get over valid DEM terrain.
        if ( std::abs(height_diff[0]) < (model.big_val()/10.0) ){
          has_intersection = true;
          break;
        }
        //if (i == 0) break; // When k*delta==0, no reason to do both + and -!

      } // End k loop
      if (has_intersection)
        break;
    } // End i loop

    // Failed to compute an intersection in the hard coded iteration limit!
    if ( !has_intersection ) {
      return Vector3();
    }

    // Refining the intersection using Levenberg-Marquardt
    // - This will actually use the L-M solver to play around with the len
    //   value to minimize the height difference from the DEM.
    int status = 0;
    Vector<double, 1> observation; observation[0] = 0;
    len = math::levenberg_marquardt(model, len, observation, status,
                                    max_abs_tol, max_rel_tol,
                                    num_max_iter
                                    );

    Vector<double, 1> dem_height = model(len);

    if ( (status < 0) || (std::abs(dem_height[0]) > height_error_tol) ){
      has_intersection = false;
      return Vector3();
    }

    has_intersection = true;
    xyz = camera_ctr + len[0]*camera_vec;
    return xyz;
  }

  namespace detail {

    /// Apply a function to evenly spaced locations along a line of pixels
    template <class FunctionT>
    void bresenham_apply( math::BresenhamLine line, size_t step,
                          FunctionT& f ) {
      while ( line.is_good() ) { // Run until the end of the line
        f( *line ); // Execute the function on this pixel location
        for ( size_t i = 0; i < step; i++ )
          ++line; // Advance "step" pixels along the line
      }
    }

    ///
    class CameraDatumBBoxHelper {
      GeoReference m_georef;
      boost::shared_ptr<camera::CameraModel> m_camera;
      Vector2 m_last_intersect;

    public:
      bool   last_valid, center_on_zero;
      BBox2  box;
      double scale;

      CameraDatumBBoxHelper( GeoReference const& georef,
                             boost::shared_ptr<camera::CameraModel> camera,
                             bool center=false) : m_georef(georef), m_camera(camera), last_valid(false), center_on_zero(center), scale( std::numeric_limits<double>::max() ) {}

      void operator() ( Vector2 const& pixel );
    }; // End class CameraDatumBBoxHelper

    /// Class to accumulate some information about a series of DEM intersections
    template <class DEMImageT>
    class CameraDEMBBoxHelper {
      GeoReference m_georef;
      boost::shared_ptr<camera::CameraModel> m_camera;
      DEMImageT    m_dem;
      Vector2      m_last_intersect;

    public:
      bool   last_valid, center_on_zero;
      BBox2  box;   ///< Bounding box containing all intersections so far.
      double scale; ///< Closest distance between two sequential intersections.
                    //   This value is in units of the georeference, either projected coords
                    //   or screwy lat/lon degrees!

      /// Constructor initializes class with DEM, camera model, etc.
      CameraDEMBBoxHelper( ImageViewBase<DEMImageT> const& dem_image,
                           GeoReference const& georef,
                           boost::shared_ptr<camera::CameraModel> camera,
                           bool center=false )
        : m_georef(georef),  m_camera(camera), m_dem(dem_image.impl()),
          last_valid(false), center_on_zero(center),
          scale( std::numeric_limits<double>::max() ) {}

      /// Intersect this pixel with the DEM and record some information about the intersection
      void operator() ( Vector2 const& pixel ) {

        bool   has_intersection;
        // TODO: Make this an input option!  Whether or not this goes outside the dem is IMPORTANT
        bool   treat_nodata_as_zero = true; // Intersect with datum if no dem
        double height_error_tol = 1e-3;   // error in DEM height
        double max_abs_tol      = 1e-14;  // abs cost function change b/w iters
        double max_rel_tol      = 1e-14;
        int    num_max_iter     = 100;
        Vector3 xyz_guess       = Vector3();
        Vector3 camera_ctr = m_camera->camera_center(pixel);  // Get ray from this pixel
        Vector3 camera_vec = m_camera->pixel_to_vector(pixel);
        Vector3 xyz // Use iterative solver call to compute an intersection of the pixel with the DEM
          = camera_pixel_to_dem_xyz(camera_ctr, camera_vec,
                                    m_dem, m_georef,
                                    treat_nodata_as_zero,
                                    has_intersection,
                                    height_error_tol, max_abs_tol, max_rel_tol,
                                    num_max_iter, xyz_guess
                                   );
        // Quit if we did not find an intersection
        if ( !has_intersection ) {
          last_valid = false;
          return;
        }
        // Use the datum to convert GCC coordinate to lon/lat/height and to a projected coordinate system
        Vector3 llh   = m_georef.datum().cartesian_to_geodetic(xyz);
        Vector2 point = m_georef.lonlat_to_point( Vector2(llh.x(), llh.y()) );

        if (!m_georef.is_projected()){
          // If we don't use a projected coordinate system, then the
          // coordinates of this point are simply lon and lat.
          // - Normalize the longitude coordinate.
          if ( center_on_zero && point[0] > 180 )
            point[0] -= 360.0;
          else if ( center_on_zero && point[0] < -180 )
            point[0] += 360.0;
          else if ( !center_on_zero && point[0] < 0 )
            point[0] += 360.0;
          else if ( !center_on_zero && point[0] > 360 )
            point[0] -= 360.0;
        }

        if ( last_valid ) { // If the last call successfully intersected
          // Compute distance from last intersection
          // - This is either in projected coordinate system units (probably meters) or screwy lat/lon degrees
          double current_scale = norm_2( point - m_last_intersect );
          if ( current_scale < scale ) // Record this distance if less than last distance
            scale = current_scale;
        }
        m_last_intersect = point; // Record this intersection
        box.grow( point ); // Expand a bounding box of all points intersected so far
        last_valid = true; // Record intersection success
      }
    }; // End class CameraDEMBBoxHelper
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

  /// Intersections that take in account DEM topography
  /// - Returns a bounding box in Georeference coordinate system (projected if available)
  ///    containing everything visible in the camera image.
  /// - Computes "scale" which is the estimated ground resolution of the camera.
  ///    This is in GeoReference measurement units (not necessarily meters!)
  template< class DEMImageT >
  BBox2 camera_bbox( ImageViewBase<DEMImageT> const& dem_image,
                     GeoReference const& georef,
                     boost::shared_ptr<vw::camera::CameraModel> camera_model,
                     int32 cols, int32 rows, float &scale ) {

    // To do: Integrate the almost identical functions camera_bbox() in
    // CameraBBox.h and CameraBBox.cc. One of them uses a DEM and the
    // second one does not.

    // Testing to see if we should be centering on zero
    bool center_on_zero = true;
    Vector3 camera_llr = // Compute lon/lat/radius of camera center
      georef.datum().cartesian_to_geodetic(camera_model->camera_center(Vector2()));
    if ( camera_llr[0] < -90 ||
         camera_llr[0] > 90 )
      center_on_zero = false;

    int32 step_amount = (2*cols+2*rows)/100;
    step_amount = std::min(step_amount, cols/4); // must have at least several points per column
    step_amount = std::min(step_amount, rows/4); // must have at least several points per row
    step_amount = std::max(step_amount, 1);      // step amount must be > 0

    // Construct helper class with DEM and camera information.
    detail::CameraDEMBBoxHelper<DEMImageT> functor( dem_image, georef, camera_model,
                                                    center_on_zero );

    // Running the edges. Note: The last valid point on a
    // BresenhamLine is the last point before the endpoint.

    // Left to right across the top side
    bresenham_apply( math::BresenhamLine(0,0,cols,0),
                     step_amount, functor );
    functor.last_valid = false;
    // Top to bottom down the right side
    bresenham_apply( math::BresenhamLine(cols-1,0,cols-1,rows),
                     step_amount, functor );
    functor.last_valid = false;
    // Right to left across the bottom side
    bresenham_apply( math::BresenhamLine(cols-1,rows-1,0,rows-1),
                     step_amount, functor );
    functor.last_valid = false;
    // Bottom to top up the left side
    bresenham_apply( math::BresenhamLine(0,rows-1,0,0),
                     step_amount, functor );
    functor.last_valid = false;

    // Do the x pattern
    bresenham_apply( math::BresenhamLine(0,0,cols-1,rows-1),
                     step_amount, functor );
    bresenham_apply( math::BresenhamLine(0,rows-1,cols-1,0),
		     step_amount, functor );
    
    // Estimate the smallest distance between adjacent points on the bounding box edges
    // - This is perhaps the finest resolution of the image.
    // - The units for this value are defined by the GeoReference and can be something weird!
    scale = functor.scale/double(step_amount);
    return functor.box;
  }

  /// Overload of camera_bbox when we don't care about getting the scale back.
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
