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


#include <vw/Cartography/CameraBBox.h>

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
  if (intersection == Vector3()){
    has_intersection = false;
    return Vector2();
  }
  
  has_intersection = true;
  
  Vector3 llh = georef.datum().cartesian_to_geodetic( intersection );
  Vector2 projected_point = georef.lonlat_to_point( Vector2( llh.x(),
                                                             llh.y() ) );

  return projected_point;
}

namespace detail {

  /// A class to help identify the extent of an image when
  /// projected onto a datum.
  class CameraDatumBBoxHelper {
    vw::cartography::GeoReference const& m_georef; // alias
    boost::shared_ptr<camera::CameraModel> m_camera;
    Vector2      m_last_intersect;
    std::vector<Vector2> *m_coords;

  public:
    bool   last_valid, center_on_zero;
    BBox2  box;
    double scale;

    CameraDatumBBoxHelper(GeoReference const& georef,
                          boost::shared_ptr<camera::CameraModel> camera,
                          bool center=false,
                          std::vector<Vector2> *coords=0):
      m_georef(georef), m_camera(camera), m_coords(coords), last_valid(false),
      center_on_zero(center), scale( std::numeric_limits<double>::max() ) {
      if (m_coords)
        m_coords->clear();
    }
    
    void operator()( Vector2 const& pixel ) {
      bool has_intersection;
      Vector2 point = geospatial_intersect(m_georef,
			                                     m_camera->camera_center(pixel),
			                                     m_camera->pixel_to_vector(pixel),
			                                     has_intersection);
      if ( !has_intersection ) {
        last_valid = false;
        return;
      }
      
      if (!m_georef.is_projected()){
        // If we don't use a projected coordinate system, then the
        // coordinates of this point are simply lon and lat.
        if ( center_on_zero && point[0] > 180 )
          point[0] -= 360.0;
      }
      if (m_coords)
        m_coords->push_back(point);
      
      if ( last_valid ) {
        double current_scale = norm_2( point - m_last_intersect );
        if ( current_scale < scale )
          scale = current_scale;
      }
      m_last_intersect = point;
      box.grow( point );
      last_valid = true;
    }
    
  }; // End class CameraDatumBBoxHelper

} // end namespace detail
  
BBox2 camera_bbox(cartography::GeoReference const& georef,
                  boost::shared_ptr<camera::CameraModel> camera_model,
                  int32 cols, int32 rows, float &scale,
                  std::vector<Vector2> *coords) {

  // Testing to see if we should be centering on zero
  bool center_on_zero = true;
  Vector3 camera_llr =
    georef.datum().cartesian_to_geodetic(camera_model->camera_center(Vector2()));
  if ( camera_llr[0] < -90 || camera_llr[0] > 90 )
    center_on_zero = false;

  int32 step_amount = (2*cols+2*rows)/100;
  step_amount = std::min(step_amount, cols/4); // ensure at least 4 pts/col
  step_amount = std::min(step_amount, rows/4); // ensure at least 4 pts/row
  step_amount = std::max(step_amount, 1);      // step amount must be > 0
  detail::CameraDatumBBoxHelper functor(georef, camera_model, center_on_zero, coords);

  // Running the edges. Note: The last valid point on a
  // BresenhamLine is the last point before the endpoint.
  bresenham_apply( BresenhamLine(0,      0,      cols,   0     ), step_amount, functor );
  functor.last_valid = false;
  bresenham_apply( BresenhamLine(cols-1, 0,      cols-1, rows  ), step_amount, functor );
  functor.last_valid = false;
  bresenham_apply( BresenhamLine(cols-1, rows-1, 0,      rows-1), step_amount, functor );
  functor.last_valid = false;
  bresenham_apply( BresenhamLine(0,      rows-1, 0,      0     ), step_amount, functor );
  functor.last_valid = false;

  // Do the x pattern
  bresenham_apply( BresenhamLine(0, 0,      cols-1, rows-1), step_amount, functor );
  bresenham_apply( BresenhamLine(0, rows-1, cols-1, 0     ), step_amount, functor );

  scale = functor.scale/double(step_amount);
  return functor.box;
}

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
    m_treat_nodata_as_zero(treat_nodata_as_zero){}
  
  /// Evaluator. See description above.
  inline result_type operator()(domain_type const& len) const {
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
  // outputs
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

    if (has_intersection) {
      break;
    }
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
    height_guess 
      = vw::cartography::demHeightGuess(pixel_cast<vw::PixelMask<float>>(dem_image));

  try {
    has_intersection = false;
    // This is a very fragile function and things can easily go wrong. 
    RayDEMIntersectionLMA model(pixel_cast<vw::PixelMask<float>>(dem_image), 
                                georef, camera_ctr,
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
  }catch(...){
    has_intersection = false;
  }
  return Vector3();
}

}} // end namespace vw::cartography
