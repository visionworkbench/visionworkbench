// __BEGIN_LICENSE__
//  Copyright (c) 2006-2025, United States Government as represented by the
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

#include <vw/Cartography/CameraBBox.h>
#include <vw/Math/BresenhamLine.h>
#include <vw/Image/Algorithms.h>
#include <vw/Image/Interpolation.h>
#include <vw/Math/LevenbergMarquardt.h>
#include <vw/Cartography/PointImageManipulation.h>

using vw::math::BresenhamLine;

namespace vw { namespace cartography {

Vector3 datum_intersection(Datum const& datum,
                           camera::CameraModel const* model,
                           Vector2 const& pix) {
  return datum_intersection(datum, model->camera_center(pix), model->pixel_to_vector(pix));
}

// Return the intersection between the ray emanating from the
// current camera pixel with the datum ellipsoid. The return value
// is a map projected point location (the intermediate between
// lon-lat-altitude and pixel).
Vector2 geospatial_intersect(vw::cartography::GeoReference const& georef,
                             Vector3 const& camera_ctr, Vector3 const& camera_vec,
                             bool& has_intersection) {

  Vector3 intersection = datum_intersection(georef.datum(), camera_ctr, camera_vec);
  if (intersection == Vector3()) {
    has_intersection = false;
    return Vector2();
  }

  has_intersection = true;
  Vector3 llh = georef.datum().cartesian_to_geodetic(intersection);
  Vector2 projected_point = georef.lonlat_to_point(Vector2(llh.x(), llh.y()));

  return projected_point;
}

namespace detail {

  /// A class to help identify the extent of an image when
  /// projected onto a datum.
  class CameraDatumHelper {
    vw::cartography::GeoReference const& m_georef; // alias
    boost::shared_ptr<camera::CameraModel> m_camera;
    Vector2      m_last_intersect;
    std::vector<Vector2> *m_coords;

  public:
    bool   last_valid, m_center_on_zero;
    BBox2  box;
    double scale;

    CameraDatumHelper(GeoReference const& georef,
                      boost::shared_ptr<camera::CameraModel> camera,
                      bool center=false,
                      std::vector<Vector2> *coords = NULL):
      m_georef(georef), m_camera(camera), m_coords(coords), last_valid(false),
      m_center_on_zero(center), scale(std::numeric_limits<double>::max()) {
      if (m_coords)
        m_coords->clear();
    }

    void operator()(Vector2 const& pixel) {
      bool has_intersection;
      Vector2 point = geospatial_intersect(m_georef,
			                                     m_camera->camera_center(pixel),
			                                     m_camera->pixel_to_vector(pixel),
			                                     has_intersection);
      if (!has_intersection) {
        last_valid = false;
        return;
      }

      if (!m_georef.is_projected()) {
        // If we don't use a projected coordinate system, then the
        // coordinates of this point are simply lon and lat.
        if (m_center_on_zero && point[0] > 180)
          point[0] -= 360.0;
      }
      if (m_coords)
        m_coords->push_back(point);

      if (last_valid) {
        double current_scale = norm_2(point - m_last_intersect);
        if (current_scale < scale)
          scale = current_scale;
      }
      m_last_intersect = point;
      box.grow(point);
      last_valid = true;
    }

  }; // End class CameraDatumHelper

  /// Class to accumulate some information about a series of DEM intersections
  class CameraDemHelper {
    GeoReference m_dem_georef, m_target_georef;
    boost::shared_ptr<camera::CameraModel> m_camera;
    vw::ImageViewRef<vw::PixelMask<float>> m_dem;
    double m_height_guess;
    Vector2 m_last_intersect;
    std::vector<Vector3> *m_coords;

  public:
    bool   m_last_valid;
    BBox2  box;   ///< Bounding box containing all intersections so far.
    double scale; ///< Closest distance between two sequential intersections, in projected coords.
    std::vector<Vector2> cam_pixels; // Collect sampled pixels here

    /// Constructor initializes class with DEM, camera model, etc.
    CameraDemHelper(vw::ImageViewRef<vw::PixelMask<float>> const& dem,
                    GeoReference const& dem_georef,
                    GeoReference const& target_georef, // return box in this projection
                    boost::shared_ptr<camera::CameraModel> camera,
                    double height_guess,
                    std::vector<Vector3> *coords = NULL): 
        m_dem_georef(dem_georef),
        m_target_georef(target_georef),
        m_camera(camera), m_dem(dem.impl()), 
        m_height_guess(height_guess),
        m_coords(coords),
        m_last_valid(false), 
        scale(std::numeric_limits<double>::max()) {
      if (m_coords)
        m_coords->clear();
    }

    // This is a function we don't want exposed outside the logic of this file,
    // as we make too many particular choices.
    static bool pix_to_pt_aux(Vector2 const& pixel,
                              vw::ImageViewRef<vw::PixelMask<float>> const& dem,
                              GeoReference const& dem_georef,
                              GeoReference const& target_georef,
                              boost::shared_ptr<camera::CameraModel> camera,
                              vw::Vector3 const& xyz_guess, double height_guess,
                              // Outputs
                              Vector2 & point, Vector3 & xyz) {

      // This is a very fragile function, and at many steps something can fail.
      try {
        // TODO: Make this an input option!  Whether or not this goes outside the dem is IMPORTANT
        bool   treat_nodata_as_zero = false; // Intersect with datum if no dem
        bool   has_intersection = false;
        double height_error_tol = 1e-3;   // error in DEM height
        double max_abs_tol      = 1e-14;  // abs cost function change b/w iters
        double max_rel_tol      = 1e-14;
        int    num_max_iter     = 100;
        Vector3 camera_ctr = camera->camera_center(pixel);  // Get ray from this pixel
        Vector3 camera_vec = camera->pixel_to_vector(pixel);

        // Use iterative solver call to compute an intersection of the pixel with the DEM
        xyz = camera_pixel_to_dem_xyz(camera_ctr, camera_vec,
                                      dem, dem_georef,
                                      treat_nodata_as_zero,
                                      has_intersection,
                                      height_error_tol, max_abs_tol, max_rel_tol,
                                      num_max_iter, xyz_guess, height_guess);
        // Quit if we did not find an intersection
        if (!has_intersection)
          return false;

        // Use the datum to convert GCC coordinate to lon/lat/height
        // and to a projected coordinate system
        Vector3 llh = target_georef.datum().cartesian_to_geodetic(xyz);
        point = target_georef.lonlat_to_point(Vector2(llh.x(), llh.y()));

        return has_intersection;
      } catch(...) {
        return false;
      }
    }

    /// Intersect this pixel with the DEM and record some information about the intersection
    void operator()(Vector2 const& pixel) {
      
      Vector2 point;
      Vector3 xyz_guess, xyz;
      bool has_intersection
        = pix_to_pt_aux(pixel, m_dem, m_dem_georef, m_target_georef,
                        m_camera, xyz_guess, m_height_guess,
                        point, xyz); // outputs

      // Quit if we did not find an intersection
      if (!has_intersection) {
        m_last_valid = false;
        return;
      }

      if (m_last_valid) {
        // If the call before this successfully intersected, compute
        // distance from last intersection.
        double current_scale = norm_2(point - m_last_intersect);
        if (current_scale < scale) // Record this distance if less than last distance
          scale = current_scale;
      }

      m_last_intersect = point; // Record this intersection

      if (m_coords)
        m_coords->push_back(xyz);
        
      box.grow(point); // Expand a bounding box of all points intersected so far
      
      m_last_valid = true; // Record intersection success
      cam_pixels.push_back(pixel);

    }
  }; // End class CameraDemHelper

  /// Apply a function to evenly spaced locations along a line of pixels
  template <class FunctionT>
  void bresenham_apply(math::BresenhamLine line, size_t step,
                        FunctionT& f) {
    while (line.is_good()) { // Run until the end of the line
      f(*line); // Execute the function on this pixel location
      for (size_t i = 0; i < step; i++)
        ++line; // Advance "step" pixels along the line
    }
  }

  // Collect valid pixel coordinates on the perimeter of the DEM,
  // and also inside using an X pattern. Some of these points may be duplicated.
  void sample_points_on_dem(vw::ImageViewRef<vw::PixelMask<float>> const& dem,
                            int dem_step,
                            std::vector<Vector2> & dem_pixels) {

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

    // Left and right edge. // Add dem_step to ensure we catch last row.
    for (int row0 = 0; row0 < dem_rows + dem_step; row0 += dem_step) {

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
    double diag = sqrt(double(dem_cols)*dem_cols + double(dem_rows)*dem_rows);
    for (double val = 0; val < diag + dem_step; val += dem_step) { // add dem_step to go past last

      int col = std::min(int(round((val/diag)*dem_cols)), dem_cols-1);
      int row = std::min(int(round((val/diag)*dem_rows)), dem_rows-1);

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

// Find a handful of valid DEM values and average them. It helps later when
// intersecting with the DEM, especially for Mars, where the DEM heights ca be
// very far from the datum.
double demHeightGuess(vw::ImageViewRef<vw::PixelMask<float>> const& dem) {

  double height_guess = 0.0;
  bool found = false;
  double sum = 0.0, num = 0.0;
  for (double row = 0; row < dem.rows(); row += dem.rows()/10.0) {
    for (double col = 0; col < dem.cols(); col += dem.cols()/10.0) {
      if (is_valid(dem(col, row))) {
        sum += dem(col, row).child();
        num++;
        if (num > 20) {
          // Those are enough, as going on for too long may take too much time
          found = true;
          break;
        }
      }
    }
    if (found) break;
  }
  if (num > 0)
    height_guess = sum/num;

  return height_guess;
} // End function demHeightGuess()

// Define an LMA model to solve for a DEM intersecting a ray. The
// variable of optimization is position on the ray. The cost
// function is difference between datum height and DEM height at
// current point on the ray.
class RayDEMIntersectionLMA:
  public math::LeastSquaresModelBase<RayDEMIntersectionLMA> {

  // TODO: Why does this use EdgeExtension if Helper() restricts access to the bounds?
  InterpolationView<EdgeExtensionView<vw::ImageViewRef<vw::PixelMask<float>>, ConstantEdgeExtension>, BilinearInterpolation> m_dem;
  vw::cartography::GeoReference const& m_georef; // alias
  Vector3      m_camera_ctr;
  Vector3      m_camera_vec;
  bool         m_treat_nodata_as_zero;

  inline vw::PixelMask<float> Helper(double x, double y) const {
    if ((0 <= x) && (x <= m_dem.cols() - 1) && // for interpolation
        (0 <= y) && (y <= m_dem.rows() - 1)) {
      vw::PixelMask<float> val = m_dem(x, y);
      if (is_valid(val))
        return val[0];
    }
    if (m_treat_nodata_as_zero)
      return 0;

    return big_val();
  }

public:
  typedef Vector<double, 1> result_type;
  typedef Vector<double, 1> domain_type;
  typedef Matrix<double>    jacobian_type; ///< Jacobian form. Auto.

  /// Return a very large error to penalize locations that fall off the edge of the DEM.
  inline double big_val() const {
    // Don't make this too big as in the LMA algorithm it may get
    // squared and may cause overflow.
    return 1.0e+50;
  }

  /// Constructor
  RayDEMIntersectionLMA(vw::ImageViewRef<vw::PixelMask<float>> const& dem_image,
                        vw::cartography::GeoReference const& georef,
                        Vector3 const& camera_ctr,
                        Vector3 const& camera_vec,
                        bool treat_nodata_as_zero):
    m_dem(interpolate(dem_image)), // create an interpolation object
    m_georef(georef),
    m_camera_ctr(camera_ctr), m_camera_vec(camera_vec),
    m_treat_nodata_as_zero(treat_nodata_as_zero) {}

  /// Evaluator. See description above.
  inline result_type operator()(domain_type const& len) const {
    
    // Check for NaN. This is a bugfix.
    if (len != len) 
      return len[0]; 
    
    // The proposed intersection point
    Vector3 xyz = m_camera_ctr + len[0]*m_camera_vec;

    // Convert to geodetic coordinates, then to DEM pixel coordinates
    Vector3 llh = m_georef.datum().cartesian_to_geodetic(xyz);
    Vector2 pix = m_georef.lonlat_to_pixel(Vector2(llh.x(), llh.y()));

    // Return a measure of the elevation difference between the DEM and the guess
    // at its current location.
    result_type result;
    result[0] = Helper(pix.x(), pix.y()) - llh[2];
    return result;
  }
};

// A very customized secant method function, for this specific application.
// Make several attempts to improve robustness. Return the solved length along
// the ray (which is also passed in as an initial guess) and the
// has_intersection flag.
void secantMethod(RayDEMIntersectionLMA & model, Vector3 const& camera_ctr,
                  Vector3 const& camera_vec, double height_error_tol,
                  // Outputs
                  bool & has_intersection,  Vector<double, 1> & len) {

  // Initialize the output
  has_intersection = false;

  // Do several attempts with different starting values for x1. The value of j
  // will control the step size.
  int num_j = 100; // will use this variable in two places below
  Vector<double, 1> len1, len2;

  for (int j = 0; j <= num_j; j++) {

    double x0 = len[0];
    double f0 = model(len)[0];

    if (std::abs(f0) < height_error_tol) {
      // As a robustness check, return early if we are already close enough.
      // Otherwise may end up with a denominator close to 0 below.
      len[0] = x0;
      has_intersection = true;
      return; // Done
    }

    double x1 = len[0] + 10.0 * (j + 1);
    len1[0] = x1;
    double f1 = model(len1)[0];

    if (std::abs(f0 - f1) < 1e-6 && j < num_j)
      continue; // Try next j, as the f values are too close

    // Do 100 iterations of secant method
    for (int i = 0; i < 100; i++) {
      double x2 = x1 - f1 * (x1 - x0) / (f1 - f0);
      len2[0] = x2;

      if (std::isnan(x2) || std::isinf(x2))
        break; // bad value, give up

      double f2 = model(len2)[0];
      if (std::abs(f2) < height_error_tol) {
        x0 = x2;
        f0 = f2;
        break;
      }
      x0 = x1; f0 = f1;
      x1 = x2; f1 = f2;
    }

    if (std::abs(f0) < height_error_tol) {
      len[0] = x0;
      has_intersection = true;
      return; // Done
    } else {
      has_intersection = false;
      // Try again
    }
  }

  return;
}

// The first step in finding where a ray intersects the DEM. Find a position
// on the ray where one is at least above the DEM. Do that by wiggling the
// initial guess along the ray.
void findInitPositionAboveDEM(RayDEMIntersectionLMA & model,
                              Vector3 const& camera_ctr,
                              Vector3 const& xyz,
                              // outputs
                              bool & has_intersection,
                              Vector<double, 1> & len) {

  Vector<double, 1> len0;
  len0[0] = norm_2(xyz - camera_ctr);

  // Initialize the outputs
  len = len0;
  has_intersection = false;

  // If the ray intersects the datum at a point which does not
  // correspond to a valid location in the DEM, wiggle that point
  // along the ray until hopefully it does.
  const double radius     = norm_2(xyz); // Radius from XZY coordinate center
  const int    ITER_LIMIT = 10; // There are two solver attempts per iteration
  const double small      = radius*0.02/(1 << ITER_LIMIT); // Wiggle
  for (int i = 0; i <= ITER_LIMIT; i++) {

    // Gradually expand delta magnitude until on final iteration it is == radius*0.02.
    // We flip between positive and negative values of ever-increasing magnitude.
    double delta = small*(1 << i);
    if (i == 0)
      delta = 0.0; // In first try, start at the initial guess

    for (int k = -1; k <= 1; k += 2) { // For k==-1, k==1
      // Use below -k because we want len to increase first time
      len[0] = len0[0] - k*delta; // Ray guess length +/- 2% planetary radius

      // Use our model to compute the height diff at this length
      Vector<double, 1> height_diff = model(len);
      // TODO: This is a very lenient threshold. big_val()/10.0 == 1.0e+49.
      // The effect of this is to stop this loop when we get over valid DEM terrain.
      if (std::abs(height_diff[0]) < (model.big_val()/10.0)) {
        has_intersection = true;
        break;
      }
    } // End k loop

    if (has_intersection)
      break;

  } // End i loop
} // End function

// Intersect the ray going from the given camera pixel with a DEM.
// The return value is a Cartesian point. If the ray goes through a
// hole in the DEM where there is no data, we return no-intersection
// or intersection with the datum, depending on whether the variable
// treat_nodata_as_zero is false or true.
Vector3 camera_pixel_to_dem_xyz(Vector3 const& camera_ctr, Vector3 const& camera_vec,
                                vw::ImageViewRef<vw::PixelMask<float>> const& dem_image,
                                vw::cartography::GeoReference const& georef,
                                bool treat_nodata_as_zero,
                                bool & has_intersection,
                                double height_error_tol,  // error in DEM height
                                double max_abs_tol, // abs cost fun change b/w iters
                                double max_rel_tol,
                                int num_max_iter,
                                Vector3 xyz_guess,
                                double height_guess) {

  // Must estimate the height guess if not provided, as otherwise the results
  // can be inaccurate.
  if (std::isnan(height_guess))
    height_guess = vw::cartography::demHeightGuess(dem_image);

  try {
    has_intersection = false;
    // This is a very fragile function and things can easily go wrong.
    RayDEMIntersectionLMA model(dem_image, georef, camera_ctr,
                                camera_vec, treat_nodata_as_zero);

    Vector3 xyz;
    if (xyz_guess == Vector3()) { // If no guess provided
      // Intersect the ray with the datum, this is a good initial guess.
      xyz = vw::cartography::datum_intersection
        (georef.datum().semi_major_axis() + height_guess,
          georef.datum().semi_minor_axis() + height_guess,
          camera_ctr, camera_vec);

      if (xyz == Vector3()) { // If we failed to intersect the datum, give up.
        has_intersection = false;
        return Vector3();
      }
    } else { // User provided guess
      xyz = xyz_guess;
    }

    // Now wiggle xyz along the ray until it is somewhere above the DEM.
    // Will return not xyz, but the 'len' along the ray for it.
    Vector<double, 1> len;
    findInitPositionAboveDEM(model, camera_ctr, xyz,
      // outputs
      has_intersection, len);
    
    // Call the secant method function to find the intersection with the
    // ground. Return has_intersection and len. This is 10x faster and more
    // robust than the Levenberg-Marquardt method used below (which used to be
    // the original method).
    Vector<double, 1> len_secant = len; // will change
    secantMethod(model, camera_ctr, camera_vec, height_error_tol,
                  has_intersection, len_secant);
    if (has_intersection) {
      xyz = camera_ctr +  len_secant[0] * camera_vec;
      return xyz;
    }

    // If no luck, fallback to using Levenberg-Marquardt, with the original
    // value of len.
    int status = 0;
    Vector<double, 1> observation; observation[0] = 0;
    len = math::levenberg_marquardt(model, len, observation, status,
      max_abs_tol, max_rel_tol, num_max_iter);
    Vector<double, 1> dem_height = model(len);

#if 0
    // We don't really care of the status of the minimization algorithm, since
    // we use it as a root solver. The good test is whether the height difference
    // is small enough. This test can fail even if we are close enough to the root.
    if (status < 0) {
      has_intersection = false;
      return Vector3();
    }
#endif

    if (std::abs(dem_height[0]) <= height_error_tol) {
        has_intersection = true;
        xyz = camera_ctr + len[0]*camera_vec;
        return xyz;
    } else {
      // Failed
      has_intersection = false;
      return Vector3();
    }
  } catch(...) {
    has_intersection = false;
  }
  return Vector3();
}

// Camera footprint on the datum.
BBox2 camera_bbox(cartography::GeoReference const& georef,
                  boost::shared_ptr<camera::CameraModel> camera_model,
                  int32 cols, int32 rows, float &scale,
                  std::vector<Vector2> *coords) {

  // Testing to see if we should be centering on zero
  bool center_on_zero = true;
  Vector3 cam_llh =
    georef.datum().cartesian_to_geodetic(camera_model->camera_center(Vector2()));
  if (cam_llh[0] < -90 || cam_llh[0] > 90)
    center_on_zero = false;

  int32 step_amount = (2*cols+2*rows)/100;
  step_amount = std::min(step_amount, cols/4); // ensure at least 4 pts/col
  step_amount = std::min(step_amount, rows/4); // ensure at least 4 pts/row
  step_amount = std::max(step_amount, 1);      // step amount must be > 0
  detail::CameraDatumHelper functor(georef, camera_model, center_on_zero, coords);

  // Running the edges. Note: The last valid point on a BresenhamLine is the
  // last point before the endpoint.
  bresenham_apply(BresenhamLine(0, 0, cols, 0), step_amount, functor);
  functor.last_valid = false;
  bresenham_apply(BresenhamLine(cols-1, 0, cols-1, rows), step_amount, functor);
  functor.last_valid = false;
  bresenham_apply(BresenhamLine(cols-1, rows-1, 0, rows-1), step_amount, functor);
  functor.last_valid = false;
  bresenham_apply(BresenhamLine(0, rows-1, 0, 0), step_amount, functor);
  functor.last_valid = false;

  // Do the x pattern
  bresenham_apply(BresenhamLine(0, 0, cols-1, rows-1), step_amount, functor);
  bresenham_apply(BresenhamLine(0, rows-1, cols-1, 0), step_amount, functor);

  scale = functor.scale/double(step_amount);
  return functor.box;
}

// Sample the image boundary. This is a helper function for camera_bbox.
void sampleImageBoundary(int cols, int rows, int num_samples, bool quick,
                         detail::CameraDemHelper& functor) {

  int32 image_step = (2*cols+2*rows) / num_samples;
  image_step = std::min(image_step, cols/4); // must have at least several points per col
  image_step = std::min(image_step, rows/4); // must have at least several points per row
  image_step = std::max(image_step, 1);      // step amount must be > 0
  
  // Running the edges. Note: The last valid point on a
  // BresenhamLine is the last point before the endpoint.

  // Left to right across the top side
  functor.m_last_valid = false;
  bresenham_apply(math::BresenhamLine(0,0,cols,0), image_step, functor);

  // Top to bottom down the right side
  functor.m_last_valid = false;
  bresenham_apply(math::BresenhamLine(cols-1,0,cols-1,rows), image_step, functor);

  // Right to left across the bottom side
  functor.m_last_valid = false;
  bresenham_apply(math::BresenhamLine(cols-1,rows-1,0,rows-1), image_step, functor);

  // Bottom to top up the left side
  functor.m_last_valid = false;
  bresenham_apply(math::BresenhamLine(0,rows-1,0,0), image_step, functor);

  if (!quick) {
    // Do the x pattern
    functor.m_last_valid = false;
    bresenham_apply(math::BresenhamLine(0,0,cols-1,rows-1), image_step, functor);
    
    functor.m_last_valid = false;
    bresenham_apply(math::BresenhamLine(0,rows-1,cols-1,0), image_step, functor);
    functor.m_last_valid = false;
  }

  return;
}

// From the DEM project points in the camera. This is important when the DEM
// is very small. 
void sampleDemBoundary(vw::ImageViewRef<vw::PixelMask<float>> const& dem,
                      GeoReference const& dem_georef,
                      GeoReference const& target_georef,
                      boost::shared_ptr<vw::camera::CameraModel> const& camera_model,
                      // TODO(oalexan1): Must be image_cols, image_rows
                      int image_cols, int image_rows, int num_samples,
                      // Outputs
                      BBox2 & cam_bbox,
                      std::vector<Vector2>& cam_pixels,
                      std::map<std::pair<double, double>, vw::Vector3>& pix2xyz) {
  
  // This output starts empty, unlike the others
  pix2xyz.clear();
  
  int dem_cols = dem.cols();
  int dem_rows = dem.rows();
  if (dem_cols <= 0 || dem_rows <= 0)
    return; // nothing to do

  // DEM sampling
  int32 dem_step = (2*dem_cols+2*dem_rows)/num_samples;
  dem_step = std::min(dem_step, dem_cols/4); // must have at least several points per col
  dem_step = std::min(dem_step, dem_rows/4); // must have at least several points per row
  dem_step = std::max(dem_step, 1);          // step amount must be > 0

  std::vector<Vector2> dem_pixels;
  vw::cartography::detail::sample_points_on_dem(dem, dem_step, dem_pixels);

  // Project the sampled points into the camera
  for (size_t it = 0; it < dem_pixels.size(); it++) {

    Vector2 lonlat, point, dem_pix, cam_pix;
    double  height;
    Vector3 llh, xyz;

    try {
      // Get the point for this DEM pixel and convert it to GCC coords
      dem_pix = dem_pixels[it];
      if (!is_valid(dem.impl()(dem_pix[0], dem_pix[1])))
        continue;

      lonlat = dem_georef.pixel_to_lonlat(dem_pix);
      height = dem.impl()(dem_pix[0], dem_pix[1]);
      point = target_georef.lonlat_to_point(lonlat);

      // Note: This height will be used further down
      llh[0] = lonlat[0]; llh[1] = lonlat[1]; llh[2] = height;

      // Note: This xyz will be used way down
      xyz = dem_georef.datum().geodetic_to_cartesian(llh);

      if (xyz == Vector3() || xyz != xyz) // watch for invalid values
        continue;

      cam_pix = camera_model->point_to_pixel(xyz);

      if (cam_pix != cam_pix)
        continue; // watch for nan

      // Check if the point projects into the camera
      if (!(cam_pix[0] >= 0 && cam_pix[0] <= image_cols-1 &&
            cam_pix[1] >= 0 && cam_pix[1] <= image_rows-1))
        continue;

      // This point looks good. Do more sanity checks

      // Geometric check. If this dot product is non-negative,
      // the point xyz is on the same side of the planet as the
      // camera center. Otherwise throw out this xyz. This is a
      // bugfix.
      vw::Vector3 cam_ctr = camera_model->camera_center(cam_pix);
      vw::Vector3 ray_vec = xyz - cam_ctr;

      double dot = dot_prod(ray_vec, -xyz);
      if (dot < 0.0)
        continue;

      // Normalize the ray vector
      double len = norm_2(ray_vec);
      ray_vec /= len;

      // If the ray from the pixel to the point on the ground
      // does not agree with the camera direction, this point is
      // spurious. The tolerance we use here is too lenient, it
      // is meant to catch only large deviations. At the same
      // time for some camera models the agreement between
      // point_to_pixel and pixel_to_vector may not be too great
      // perhaps if a numerical solver is used, hence the
      // tolerance is not made too tight.
      vw::Vector3 camera_dir = camera_model->pixel_to_vector(cam_pix);
      double DIRECTION_TOLERANCE = 1e-3;
      if (norm_2(camera_dir - ray_vec) > DIRECTION_TOLERANCE)
        continue;

      // Finally a good point we can accept
      cam_bbox.grow(point);

      // Add to cam_pixels from this different way of sampling
      cam_pixels.push_back(cam_pix);

      pix2xyz[std::make_pair(cam_pix.x(), cam_pix.y())] = xyz;
    } catch(...) {
      continue;
    }
  } // End loop through points on the DEM
  
  return;
}

// Estimate the gsd, in point units, by projecting onto the ground neighboring points
double calcMeanGsd(int cols, int rows,
                   detail::CameraDemHelper const& functor,
                   std::vector<Vector2> const& cam_pixels,
                   std::map<std::pair<double, double>, Vector3> const& pix2xyz,
                   vw::ImageViewRef<vw::PixelMask<float>> const& dem,
                   GeoReference const& dem_georef,
                   GeoReference const& target_georef,
                   boost::shared_ptr<vw::camera::CameraModel> const& camera_model,
                   double height_guess) {

  std::vector<double> gsd;
  BBox2i image_box(0, 0, cols, rows);
  for (size_t it = 0; it < cam_pixels.size(); it++) {

    // Note how we cast to int
    Vector2i ctr_pix = cam_pixels[it];
    if (!image_box.contains(ctr_pix))
      continue;

    Vector2 ctr_point;
    Vector3 xyz_guess, xyz;

    // If we ran into into this pixel before, we have an idea of its xyz.
    // This will help with intersections below. Use here the original
    // pixel, not the one cast to int.
    auto coord_pair = std::make_pair(cam_pixels[it].x(), cam_pixels[it].y());
    auto xyz_it = pix2xyz.find(coord_pair);
    if (xyz_it != pix2xyz.end()) {
      xyz_guess = xyz_it->second;
    }

    bool has_intersection
      = functor.pix_to_pt_aux(ctr_pix, dem, dem_georef, target_georef,  camera_model, 
                              xyz_guess, height_guess,
                              ctr_point, xyz); // outputs

    if (!has_intersection)
      continue;

    // Four neighboring pixels. Use the same guess.
    for (int j = 0; j < 4; j++) {
      Vector2i off_pix = ctr_pix;
      if (j == 0) off_pix += Vector2i(1, 0);
      if (j == 1) off_pix += Vector2i(0, 1);
      if (j == 2) off_pix += Vector2i(-1, 0);
      if (j == 3) off_pix += Vector2i(0, -1);
      if (!image_box.contains(off_pix))
        continue;

      Vector2 off_point;
      bool has_intersection
        = functor.pix_to_pt_aux(off_pix, dem, dem_georef, target_georef, camera_model,
                                xyz_guess, height_guess,
                                off_point, xyz); // outputs
      if (!has_intersection)
        continue;

      gsd.push_back(norm_2(ctr_point-off_point));
    }
  }

  VW_ASSERT(!gsd.empty(), ArgumentErr() << "Could not sample correctly the image.");

  // Note that, at least for LRO NAC, the GSD in row and column direction can be
  // wildly different (not true for WV though). Hence we should do an average,
  // not a median. But first trimming some outliers.
  // TODO(oalexan1): Need to do the same for point2dem.
  std::sort(gsd.begin(), gsd.end()); // in order
  int gsd_len = gsd.size();
  int beg = int(0.1*gsd_len);
  int end = int(0.9*gsd_len);
  if (gsd_len > 0 && beg == end) {
    // When there are not enough samples, use all of them
    beg = 0; 
    end = gsd_len;
  }
  VW_ASSERT(beg < end, ArgumentErr() << "Could not sample correctly the image.");

  // Average
  double mean_gsd = 0.0;
  int num = 0;
  for (int it = beg; it < end; it++) {
    double val = gsd[it];
    if (val <= 0 || val != val)
      continue;
    mean_gsd += val;
    num  += 1;
  }

  if (num == 0)
    VW_ASSERT(beg < end, ArgumentErr() << "Could not sample correctly the image.");
  mean_gsd /= num;
  
  return mean_gsd;
}

// Camera footprint on the DEM. See the .h file for details.
// TODO(oalexan1): This needs to be broken up into image-to-ground
// estimation and ground-to-image estimation.
BBox2 camera_bbox(vw::ImageViewRef<vw::PixelMask<float>> const& dem,
                  GeoReference const& dem_georef,
                  GeoReference const& target_georef, // return box in this projection
                  boost::shared_ptr<vw::camera::CameraModel> camera_model,
                  int32 cols, int32 rows, float &mean_gsd,
                  bool quick,
                  std::vector<Vector3> *coords,
                  int num_samples) {
  
  // An estimate of the DEM height above the datum
  double height_guess = vw::cartography::demHeightGuess(dem);

  // Testing to see if we should be centering on zero. The logic here is consistent
  // with point2dem.
  bool center_on_zero = true;
  Vector3 cam_llh = // Compute lon/lat/height of camera center
    target_georef.datum().cartesian_to_geodetic(camera_model->camera_center(Vector2()));
  if (cam_llh[0] < -90 || cam_llh[0] > 90)
    center_on_zero = false;

  // Construct helper class with DEM and camera information.
  detail::CameraDemHelper functor(dem, dem_georef, target_georef,
                                   camera_model, height_guess, coords);

  // Image sampling. About 1000 samples are needed to not cut the corners
  // for complex geometry.
  sampleImageBoundary(cols, rows, num_samples, quick, functor);
  
  // The bounding box collected so far.
  BBox2 cam_bbox = functor.box;

  // Sampled camera pixels collected so far
  std::vector<Vector2> cam_pixels = functor.cam_pixels;

  // From the DEM project points in the camera. This is important when the DEM
  // is very small.
  std::map<std::pair<double, double>, vw::Vector3> pix2xyz;
  if (!quick)
    sampleDemBoundary(dem, dem_georef, target_georef, camera_model, cols, rows,
                      num_samples, cam_bbox, cam_pixels, pix2xyz);

  // Find the mean GSD
  mean_gsd = calcMeanGsd(cols, rows, functor, cam_pixels, pix2xyz,
                         dem, dem_georef, target_georef, camera_model,
                         height_guess);

  return cam_bbox;
}

// Overload with no scale return
BBox2 camera_bbox(GeoReference const& dem_georef,
                  boost::shared_ptr<vw::camera::CameraModel> camera_model,
                  int32 cols, int32 rows) {
  float scale = -1.0; // will be updated
  return camera_bbox(dem_georef, camera_model, cols, rows, scale);
}

// Overload of camera_bbox when we don't care about getting the mean_gsd back.
BBox2 camera_bbox(vw::ImageViewRef<vw::PixelMask<float>> const& dem,
                  GeoReference const& dem_georef,
                  GeoReference const& target_georef,
                  boost::shared_ptr<vw::camera::CameraModel> camera_model,
                  int32 cols, int32 rows) {

  float mean_gsd = -1.0; // will be updated
  return camera_bbox(dem, dem_georef, target_georef, camera_model, cols, rows, mean_gsd);
}

}} // end namespace vw::cartography
