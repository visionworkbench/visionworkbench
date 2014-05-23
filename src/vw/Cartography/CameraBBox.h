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
#if defined(VW_HAVE_PKG_CAMERA)

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
                                  double max_abs_tol      = 1e-14, // abs cost fun change b/w iters
                                  double max_rel_tol      = 1e-14,
                                  int num_max_iter        = 100,
                                  Vector3 xyz_guess       = Vector3()
                                  ){
    
    // This is a very fragile function and things can easily go wrong. 
    try {
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
    }catch(...){
      has_intersection = false;
    }
    return Vector3();
  }

  namespace detail {

    // TODO: This should be done by default!
    // Normalize the coordinate if lonlat
    void recenter_point(bool center_on_zero, GeoReference const& georef, Vector2 & point);
                        
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

    /// Class to accumulate some information about a series of DEM intersections
    template <class DEMImageT>
    class CameraDEMBBoxHelper {
      GeoReference m_dem_georef, m_target_georef;
      boost::shared_ptr<camera::CameraModel> m_camera;
      DEMImageT    m_dem;
      Vector2      m_last_intersect;
      std::vector<Vector3> *m_coords;

    public:
      bool   m_last_valid, m_center_on_zero;
      BBox2  box;   ///< Bounding box containing all intersections so far.
      double scale; ///< Closest distance between two sequential intersections, in projected coords.
      std::vector<Vector2> cam_pixels; // Collect sampled pixels here
      
      /// Constructor initializes class with DEM, camera model, etc.
      CameraDEMBBoxHelper( ImageViewBase<DEMImageT> const& dem,
                           GeoReference const& dem_georef,
                           GeoReference const& target_georef, // return box in this projection
                           boost::shared_ptr<camera::CameraModel> camera,
                           bool center_on_zero,
                           std::vector<Vector3> *coords=0)
        : m_dem_georef(dem_georef),
          m_target_georef(target_georef), 
          m_camera(camera), m_dem(dem.impl()), m_coords(coords),
          m_last_valid(false), m_center_on_zero(center_on_zero),
          scale( std::numeric_limits<double>::max() ) {
        if (m_coords)
          m_coords->clear();
      }

      // This is a function we don't want exposed outside the logic of this file,
      // as we make too many particular choices.
      static bool camera_pixel_to_dem_point(Vector2 const& pixel,
                                            ImageViewBase<DEMImageT> const& dem,
                                            GeoReference const& dem_georef,
                                            GeoReference const& target_georef, 
                                            boost::shared_ptr<camera::CameraModel> camera,
                                            bool center_on_zero,
                                            Vector2 & point, // output
                                            Vector3 & xyz){

        // This is a very fragile function, and at many steps something can fail.
        try {
          // TODO: Make this an input option!  Whether or not this goes outside the dem is IMPORTANT
          bool   treat_nodata_as_zero = false; // Intersect with datum if no dem
          
          bool   has_intersection = false;
          double height_error_tol = 1e-3;   // error in DEM height
          double max_abs_tol      = 1e-14;  // abs cost function change b/w iters
          double max_rel_tol      = 1e-14;
          int    num_max_iter     = 100;
          Vector3 xyz_guess       = Vector3();
          Vector3 camera_ctr = camera->camera_center(pixel);  // Get ray from this pixel
          Vector3 camera_vec = camera->pixel_to_vector(pixel);
          
          // Use iterative solver call to compute an intersection of the pixel with the DEM	
          xyz = camera_pixel_to_dem_xyz(camera_ctr, camera_vec,
                                        dem, dem_georef,
                                        treat_nodata_as_zero,
                                        has_intersection,
                                        height_error_tol, max_abs_tol, max_rel_tol,
                                        num_max_iter, xyz_guess
                                       );
          // Quit if we did not find an intersection
          if (!has_intersection)
            return false;
          
          // Use the datum to convert GCC coordinate to lon/lat/height
          // and to a projected coordinate system
          Vector3 llh = target_georef.datum().cartesian_to_geodetic(xyz);
          point = target_georef.lonlat_to_point( Vector2(llh.x(), llh.y()) );
          recenter_point(center_on_zero, target_georef, point);
                    
          return has_intersection;
        }catch(...){
          return false;
        }
      }
      
      /// Intersect this pixel with the DEM and record some information about the intersection
      void operator() ( Vector2 const& pixel ) {

        Vector2 point;
        Vector3 xyz;
        bool has_intersection = camera_pixel_to_dem_point(pixel, m_dem, m_dem_georef,
                                                          m_target_georef,
                                                          m_camera, m_center_on_zero,  
                                                          point, // output
                                                          xyz);
        // Quit if we did not find an intersection
        if ( !has_intersection ) {
          m_last_valid = false;
          return;
        }
        
        if ( m_last_valid ) {
          // If the call before this successfully intersected, compute
          // distance from last intersection.
          double current_scale = norm_2( point - m_last_intersect );
          if ( current_scale < scale ) // Record this distance if less than last distance
            scale = current_scale;
        }
        
        m_last_intersect = point; // Record this intersection
        
        if (m_coords)
          m_coords->push_back(xyz);
        
        box.grow( point ); // Expand a bounding box of all points intersected so far
        m_last_valid = true; // Record intersection success
        cam_pixels.push_back(pixel);
        
      }
    }; // End class CameraDEMBBoxHelper

    // Collect valid pixel coordinates on the perimeter of the DEM,
    // and also inside using an X pattern. Some of these points may be duplicated. 
    template <class DEMImageT>
    void sample_points_on_dem(DEMImageT const& dem, int dem_step, 
                              std::vector<Vector2> & dem_pixels){
      
      dem_pixels.clear();
      
      int dem_cols = dem.impl().cols(); 
      int dem_rows = dem.impl().rows();
      if (dem_cols <= 0 || dem_rows <= 0) return; // nothing to do

      // top and bottom edge
      for (int col0 = 0; col0 < dem_cols + dem_step; // add dem_step to ensure we catch last col
           col0 += dem_step) {
      
        int col = std::min(col0, dem_cols - 1); // the last will be dem_cols-1

        // go down
        for (int row = 0; row < dem_rows; row++) {
          if (is_valid(dem.impl()(col, row))) {
            // Stop at the first valid pixel
            dem_pixels.push_back(Vector2(col, row));
            break;
          }
        }

        // go up
        for (int row = dem_rows-1; row >= 0; row--) {
          if (is_valid(dem.impl()(col, row))) {
            dem_pixels.push_back(Vector2(col, row));
            break;
          }
        }
      
      }

      // left and right edge
      for (int row0 = 0; row0 < dem_rows + dem_step; // add dem_step to ensure we catch last row
           row0 += dem_step) {
        
        int row = std::min(row0, dem_rows - 1); // the last will be dem_rows-1

        // go right
        for (int col = 0; col < dem_cols; col++) {
          if (is_valid(dem.impl()(col, row))) {
            // Hit first valid pixel from the left
            dem_pixels.push_back(Vector2(col, row));
            break;
          }
        }

        // go left
        for (int col = dem_cols-1; col >= 0; col--) {
          if (is_valid(dem.impl()(col, row))) {
            // Hit first valid pixel from the left
            dem_pixels.push_back(Vector2(col, row));
            break;
          }
        }
      
      }

      // Go on the diagonals. 
      double diag = sqrt( double(dem_cols)*dem_cols + double(dem_rows)*dem_rows);
      for (double val = 0; val < diag + dem_step; val += dem_step) { // add dem_step to go past last

        int col = std::min( int(round((val/diag)*dem_cols)), dem_cols-1);
        int row = std::min( int(round((val/diag)*dem_rows)), dem_rows-1);

        // Main diagonal
        if (is_valid(dem.impl()(col, row))) {
          dem_pixels.push_back(Vector2(col, row));
        }

        // Other diagonal
        col = dem_cols - 1 - col;
        if (is_valid(dem.impl()(col, row))) {
          dem_pixels.push_back(Vector2(col, row));
        }
      
      }
    }
    
  } // end namespace detail

  // Functions for Users
  //////////////////////////////////////////////////////

  // Simple Intersection interfaces
  
  /// Compute the bounding box in points (georeference space) that is
  /// defined by georef. Scale is MPP as georeference space is in meters.
  /// - If coords is provided the intersection coordinates will be stored there.
  BBox2 camera_bbox( GeoReference const& georef,
                     boost::shared_ptr<vw::camera::CameraModel> camera_model,
                     int32 cols, int32 rows, float &scale,
                     std::vector<Vector2> *coords=0 );

  /// Overload with no scale return
  inline BBox2 camera_bbox( GeoReference const& dem_georef,
                            boost::shared_ptr<vw::camera::CameraModel> camera_model,
                            int32 cols, int32 rows ) {
    float scale;
    return camera_bbox( dem_georef, camera_model, cols, rows, scale );
  }

  /// Intersections that take into account DEM topography
  /// - Returns a bounding box in Georeference coordinate system (projected if available)
  ///    containing everything visible in the camera image.
  /// - Computes mean_gsd which is the estimated mean ground resolution of the camera.
  ///   Note that ground resolution in row and col directions can be different for LRO NAC.
  ///   This will just return a mean of the two. 
  ///   The mean_gsd is in GeoReference measurement units (not necessarily meters!)
  /// - If the quick option is enabled, only rays along the image borders will be used
  ///   to perform the computation.
  /// - If coords is provided the intersection coordinates will be stored there.
  template< class DEMImageT >
  BBox2 camera_bbox( ImageViewBase<DEMImageT> const& dem,
                     GeoReference const& dem_georef,
                     GeoReference const& target_georef, // return box in this projection
                     boost::shared_ptr<vw::camera::CameraModel> camera_model,
                     int32 cols, int32 rows, float &mean_gsd,
                     bool quick=false,
                     std::vector<Vector3> *coords=0 ) {

    // Testing to see if we should be centering on zero
    bool center_on_zero = true;
    Vector3 camera_llr = // Compute lon/lat/radius of camera center
      target_georef.datum().cartesian_to_geodetic(camera_model->camera_center(Vector2()));
    if ( camera_llr[0] < -90 || camera_llr[0] > 90 )
      center_on_zero = false;

    int dem_cols = dem.impl().cols(); 
    int dem_rows = dem.impl().rows();
    if (dem_cols <= 0 || dem_rows <= 0)
      return BBox2(); // nothing to do

    int NUM_SAMPLES = 1000; // increased from 100, which was cutting corners
    
    // Image sampling
    int32 image_step = (2*cols+2*rows)/NUM_SAMPLES;
    image_step = std::min(image_step, cols/4); // must have at least several points per col
    image_step = std::min(image_step, rows/4); // must have at least several points per row
    image_step = std::max(image_step, 1);      // step amount must be > 0

    // DEM sampling
    int32 dem_step = (2*dem_cols+2*dem_rows)/NUM_SAMPLES;
    dem_step = std::min(dem_step, dem_cols/4); // must have at least several points per col
    dem_step = std::min(dem_step, dem_rows/4); // must have at least several points per row
    dem_step = std::max(dem_step, 1);          // step amount must be > 0

    // Construct helper class with DEM and camera information.
    detail::CameraDEMBBoxHelper<DEMImageT> functor( dem, dem_georef, target_georef,
                                                    camera_model, center_on_zero, coords );

    // Running the edges. Note: The last valid point on a
    // BresenhamLine is the last point before the endpoint.

    // Left to right across the top side
    functor.m_last_valid = false;
    bresenham_apply( math::BresenhamLine(0,0,cols,0), image_step, functor );
    
    // Top to bottom down the right side
    functor.m_last_valid = false;
    bresenham_apply( math::BresenhamLine(cols-1,0,cols-1,rows), image_step, functor );

    // Right to left across the bottom side
    functor.m_last_valid = false;
    bresenham_apply( math::BresenhamLine(cols-1,rows-1,0,rows-1), image_step, functor );

    // Bottom to top up the left side
    functor.m_last_valid = false;
    bresenham_apply( math::BresenhamLine(0,rows-1,0,0), image_step, functor );

    if (!quick) {
      // Do the x pattern
      functor.m_last_valid = false;
      bresenham_apply( math::BresenhamLine(0,0,cols-1,rows-1), image_step, functor );

      functor.m_last_valid = false;
      bresenham_apply( math::BresenhamLine(0,rows-1,cols-1,0), image_step, functor );
      functor.m_last_valid = false;
    }
    
    // Estimate the smallest distance between adjacent points on the bounding box edges
    // We in fact want the average distance, so we will re-estimate that below.
    mean_gsd = functor.scale/double(image_step);

    // The bounding box collected so far. 
    BBox2 cam_bbox = functor.box;

    // Sampled camera pixels collected so far 
    std::vector<Vector2> cam_pixels = functor.cam_pixels; 

    if (!quick) {

      //vw_out() << "Computed image to DEM bbox: " << cam_bbox << std::endl;

      // Bugfix. Traversing the bbox of the image and drawing an X on
      // its diagonals is not enough sometimes to accurately determine
      // where the map-projected image overlaps with the DEM. It fails
      // if the DEM is small. Therefore, also do the reverse, from the
      // DEM project points in the camera, traversing the bbox of the
      // DEM and doing an X pattern, and see which fall inside.
      std::vector<Vector2> dem_pixels;
      detail::sample_points_on_dem(dem, dem_step, dem_pixels);
        
      // Project the sampled points into the camera
      for (size_t it = 0; it < dem_pixels.size(); it++) {

        Vector2 lonlat, point, dem_pix, cam_pix;
        double  height;
        Vector3 llh, xyz;

        try {
          // Get the point for this DEM pixel and convert it to GCC coords
          dem_pix = dem_pixels[it];
          if (!is_valid(dem.impl()(dem_pix[0], dem_pix[1])))
            continue; // redundant
          lonlat = dem_georef.pixel_to_lonlat(dem_pix);
          height = dem.impl()(dem_pix[0], dem_pix[1]);

          point = target_georef.lonlat_to_point(lonlat);
          detail::recenter_point(center_on_zero, target_georef, point);
          
          llh[0] = lonlat[0]; llh[1] = lonlat[1]; llh[2] = height;

          xyz = dem_georef.datum().geodetic_to_cartesian(llh);
          if (xyz == Vector3() || xyz != xyz) // watch for invalid values
            continue;

          cam_pix = camera_model->point_to_pixel(xyz);
          if (cam_pix != cam_pix)
            continue; // watch for nan
	
          if (cam_pix[0] >= 0 && cam_pix[0] <= cols-1 &&
              cam_pix[1] >= 0 && cam_pix[1] <= rows-1 ) {

            // Finally a good point we can accept
            cam_bbox.grow(point);
            //vw_out() << "cam_pix: " << cam_pix << std::endl;
            //vw_out() << "point: " << point << std::endl;
            //vw_out() << "llh: " << llh << std::endl;

            // Add to cam_pixels from this different way of sampling
            cam_pixels.push_back(cam_pix);
          }
        }
        catch(...) {
          // It is possible to hit exceptions in here from coordinate transformation and such which
          //  do not cause further problems, for example with points on large DEMs that do not fit
          //  well into the target georef.  We can safely skip these since they probably don't intersect
          //  the image anyways.
          continue;
        }  
      } // End loop through points on the DEM
      
      //vw_out() << "Expanded bbox with DEM to image: " << cam_bbox << std::endl;
    } // End if (!quick)

    // Now estimate the gsd, in point units, by projecting onto the ground neighboring points
    std::vector<double> gsd;
    BBox2i image_box(0, 0, cols, rows);
    for (size_t it = 0; it < cam_pixels.size(); it++) {

      Vector2i ctr_pix = cam_pixels[it];
      if (!image_box.contains(ctr_pix))
        continue;

      Vector2 ctr_point;
      Vector3 xyz;
      bool has_intersection = functor.camera_pixel_to_dem_point(ctr_pix, dem, dem_georef,
                                                                target_georef,  
                                                                camera_model, center_on_zero,  
                                                                ctr_point, // output
                                                                xyz
                                                               );
      if ( !has_intersection )
        continue; 
      
      // Four neighboring pixels
      for (int j = 0; j < 4; j++) {
        Vector2i off_pix = ctr_pix;
        if (j == 0) off_pix += Vector2i(1, 0);
        if (j == 1) off_pix += Vector2i(0, 1);
        if (j == 2) off_pix += Vector2i(-1, 0);
        if (j == 3) off_pix += Vector2i(0, -1);
        if (!image_box.contains(off_pix))
          continue;
        
        Vector2 off_point;
        bool has_intersection = functor.camera_pixel_to_dem_point(off_pix, dem, dem_georef,
                                                                  target_georef,  
                                                                  camera_model, center_on_zero,  
                                                                  off_point, // output
                                                                  xyz
                                                                 );
        if ( !has_intersection )
          continue; 

        gsd.push_back(norm_2(ctr_point-off_point));
      }
    }

    VW_ASSERT(!gsd.empty(), ArgumentErr() << "Could not sample correctly the image.");

    // Note that, at least for LRO NAC, the GSD in row and column
    // direction can be wildly different (not true for WV
    // though). Hence we should do an average, not a median. But first
    // trimming some outliers.
    std::sort(gsd.begin(), gsd.end()); // in order
    int gsd_len = gsd.size();
    int beg = int(0.1*gsd_len);
    int end = int(0.9*gsd_len);
    VW_ASSERT(beg < end, ArgumentErr() << "Could not sample correctly the image.");

    mean_gsd = 0;
    int num = 0;
    for (int it = beg; it < end; it++) {
      double val = gsd[it];
      if (val <= 0 || val != val)
        continue;
      mean_gsd += val;
      num   += 1;
    }

    if (num == 0)
      VW_ASSERT(beg < end, ArgumentErr() << "Could not sample correctly the image.");

    mean_gsd /= num;
    
    return cam_bbox;
  }

  /// Overload of camera_bbox when we don't care about getting the mean_gsd back.
  template< class DEMImageT >
  inline BBox2 camera_bbox( ImageViewBase<DEMImageT> const& dem,
                            GeoReference const& dem_georef,
                            GeoReference const& target_georef,
                            boost::shared_ptr<vw::camera::CameraModel> camera_model,
                            int32 cols, int32 rows) {
    float mean_gsd;
    return camera_bbox<DEMImageT>( dem.impl(), dem_georef, target_georef,
                                    camera_model, cols, rows, mean_gsd );
  }

} // namespace cartography
} // namespace vw

#endif // VW_HAVE_PKG_CAMERA

#endif // __VW_CARTOGRAPHY_CAMERABBOX_H__
