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


/// \file CameraUtilities.cc
///
/// This file contains miscellaneous functions for working with camera models.
///
#include <vw/Camera/CameraUtilities.h>
#include <vw/Math/Geometry.h>

namespace vw {
namespace camera {

// Unpack a vector into a rotation + translation + scale
void vector_to_transform(Vector<double> const & C, 
			 Matrix3x3            & rotation,
			 Vector3              & translation,
			 double               & scale){

  if (C.size() != 7) 
    vw_throw( LogicErr() << "Expecting a vector of size 7.\n" );
  
  translation = subvector(C, 0, 3);
  rotation    = axis_angle_to_quaternion(subvector(C, 3, 3)).rotation_matrix();
  scale       = C[6];
  
  return;
}

// Pack a rotation + translation + scale into a vector
void transform_to_vector(Vector<double>  & C, 
			 Matrix3x3 const & rotation,
			 Vector3   const & translation,
			 double    const & scale){
  
  Vector3 axis_angle = Quat(rotation).axis_angle();

  C.set_size(7);
  subvector(C, 0, 3) = translation;
  subvector(C, 3, 3) = axis_angle;
  C[6] = scale;
  
  return;
}

/// Find the camera model that best projects given xyz points into given pixels.
/// If a positive camera weight is specified, use that as a constraint to not
/// move the camera center too much during optimization.
template <class CAM>
class CameraSolveLMA : public vw::math::LeastSquaresModelBase<CameraSolveLMA<CAM> > {
  std::vector<vw::Vector3> const& m_xyz;
  vw::Vector3 m_camera_center;
  double m_camera_weight;
  CAM m_camera_model;
  mutable size_t m_iter_count;
  
public:

  typedef vw::Vector<double>    result_type;   // pixel residuals
  typedef vw::Vector<double, 6> domain_type;   // camera parameters (camera center and axis angle)
  typedef vw::Matrix<double> jacobian_type;

  /// Instantiate the solver with a set of xyz to pixel pairs and a pinhole model
  CameraSolveLMA(std::vector<vw::Vector3> const& xyz,
                 CAM const& camera_model, double camera_weight):
    m_xyz(xyz), m_camera_model(camera_model), m_iter_count(0),
    m_camera_center(camera_model.camera_center()), m_camera_weight(camera_weight){}

  /// Given the camera, project xyz into it
  inline result_type operator()(domain_type const& C) const {

    // Create the camera model
    CAM camera_model = m_camera_model; // make a copy local to this function
    vector_to_camera(camera_model, C);      // update its parameters
    
    size_t result_size = 2*m_xyz.size();
    if (m_camera_weight > 0)
      result_size += 3;

    result_type result;
    result.set_size(result_size);
    for (size_t i = 0; i < m_xyz.size(); i++) {
      Vector2 pixel = camera_model.point_to_pixel(m_xyz[i]);
      result[2*i  ] = pixel[0];
      result[2*i+1] = pixel[1];
    }

    if (m_camera_weight > 0){
      for (size_t i = 0; i < 3; i++){
        int j = 2*m_xyz.size() + i;
        result[j] = m_camera_weight*(C[i] - m_camera_center[i]);
      }
    }
    
    ++m_iter_count;
    
    return result;
  }
}; // End class CameraSolveLMA
  
/// Adjust a given camera so that the xyz points project project to
/// the pixel values.
void fit_camera_to_xyz(std::string const& camera_type,
		       bool refine_camera, 
		       std::vector<Vector3> const& xyz_vec,
		       std::vector<double> const& pixel_values,
		       bool verbose, boost::shared_ptr<CameraModel> & out_cam){
  
  // Create fake points in space at given distance from this camera's
  // center and corresponding actual points on the ground.  Use 500
  // km, just some height not too far from actual satellite height.
  const double ht = 500000; 
  vw::Matrix<double> in, out;
  int num_pts = pixel_values.size()/2;
  in.set_size(3, num_pts);
  out.set_size(3, num_pts);
  for (int col = 0; col < in.cols(); col++) {
    Vector3 a = out_cam->camera_center(Vector2(0, 0)) +
      ht*out_cam->pixel_to_vector(Vector2(pixel_values[2*col], pixel_values[2*col+1]));
    for (int row = 0; row < in.rows(); row++) {
      in(row, col)  = a[row];
      out(row, col) = xyz_vec[col][row];
    }
  }

  // Apply a transform to the camera so that the fake points are on top of the real points
  Matrix<double, 3, 3> rotation;
  Vector3 translation;
  double scale;
  find_3D_transform(in, out, rotation, translation, scale);
  if (camera_type == "opticalbar")
    ((vw::camera::OpticalBarModel*)out_cam.get())->apply_transform(rotation,
								    translation, scale);
  else
    ((PinholeModel*)out_cam.get())->apply_transform(rotation, translation, scale);

  // Print out some errors
  if (verbose) {
    vw_out() << "The error between the projection of each ground "
	     << "corner point into the coarse camera and its pixel value:\n";
    for (size_t corner_it = 0; corner_it < num_pts; corner_it++) {
      vw_out () << "Corner and error: ("
		<< pixel_values[2*corner_it] << ' ' << pixel_values[2*corner_it+1]
		<< ") " <<  norm_2(out_cam.get()->point_to_pixel(xyz_vec[corner_it]) -
				 Vector2( pixel_values[2*corner_it],
					  pixel_values[2*corner_it+1]))
		<< std::endl;
    }
  }
  
  // Solve a little optimization problem to make the points on the ground project
  // as much as possible exactly into the image corners.
  if (refine_camera) {
    Vector<double> pixel_vec; // must copy to this structure
    pixel_vec.set_size(pixel_values.size());
    for (size_t corner_it = 0; corner_it < pixel_vec.size(); corner_it++) 
      pixel_vec[corner_it] = pixel_values[corner_it];
    const double abs_tolerance  = 1e-24;
    const double rel_tolerance  = 1e-24;
    const int    max_iterations = 2000;
    int status = 0;
    Vector<double> final_params;
    Vector<double> seed;
    double camera_weight = 0; // if camera_weight > 0, need to add 3 zero values to pixel_vec
    
    if (camera_type == "opticalbar") {
      CameraSolveLMA<vw::camera::OpticalBarModel>
	lma_model(xyz_vec, *((vw::camera::OpticalBarModel*)out_cam.get()), camera_weight);
      camera_to_vector(*((vw::camera::OpticalBarModel*)out_cam.get()), seed);
      final_params = math::levenberg_marquardt(lma_model, seed, pixel_vec,
					       status, abs_tolerance, rel_tolerance,
					       max_iterations);
      vector_to_camera(*((vw::camera::OpticalBarModel*)out_cam.get()), final_params);
    } else {
      CameraSolveLMA<PinholeModel> lma_model(xyz_vec, *((PinholeModel*)out_cam.get()),
                                             camera_weight);
      camera_to_vector(*((PinholeModel*)out_cam.get()), seed);
      final_params = math::levenberg_marquardt(lma_model, seed, pixel_vec,
					       status, abs_tolerance, rel_tolerance,
					       max_iterations);
      vector_to_camera(*((PinholeModel*)out_cam.get()), final_params);
    }
    if (status < 1)
      vw_out() << "The Levenberg-Marquardt solver failed. Results may be inaccurate.\n";

    if (verbose) {
      vw_out() << "The error between the projection of each ground "
	       << "corner point into the refined camera and its pixel value:\n";
      for (size_t corner_it = 0; corner_it < num_pts; corner_it++) {
	vw_out () << "Corner and error: ("
		  << pixel_values[2*corner_it] << ' ' << pixel_values[2*corner_it+1]
		  << ") " <<  norm_2(out_cam.get()->point_to_pixel(xyz_vec[corner_it]) -
				     Vector2( pixel_values[2*corner_it],
					      pixel_values[2*corner_it+1]))
		  << std::endl;
      }
    }
    
  } // End camera refinement case
}
  
/// Load a pinhole camera model of any supported type
boost::shared_ptr<vw::camera::CameraModel>
load_pinhole_camera_model(std::string const& path){

  std::string lcase_file = boost::to_lower_copy(path);
  if (boost::ends_with(lcase_file,".cahvore") ) {
    return boost::shared_ptr<vw::camera::CameraModel>( new vw::camera::CAHVOREModel(path) );
  } else if (boost::ends_with(lcase_file,".cahvor") ||
             boost::ends_with(lcase_file,".cmod"  )   ) {
    return boost::shared_ptr<vw::camera::CameraModel>( new vw::camera::CAHVORModel(path) );
  } else if ( boost::ends_with(lcase_file,".cahv") ||
              boost::ends_with(lcase_file,".pin" )   ) {
    return boost::shared_ptr<vw::camera::CameraModel>( new vw::camera::CAHVModel(path) );
  } else if ( boost::ends_with(lcase_file,".pinhole") ||
              boost::ends_with(lcase_file,".tsai"   )   ) {
    return boost::shared_ptr<vw::camera::CameraModel>( new vw::camera::PinholeModel(path) );
  } else {
    vw::vw_throw(vw::ArgumentErr() << "PinholeStereoSession: unsupported camera file type.\n");
  }
}


/// Load a pinhole, CAHV, CAHVOR, or CAHVORE model and convert to CAHV.
boost::shared_ptr<vw::camera::CAHVModel>
load_cahv_pinhole_camera_model(std::string const& image_path,
                               std::string const& camera_path){
  // Get the image size
  vw::DiskImageView<float> disk_image(image_path);
  vw::Vector2i image_size(disk_image.cols(), disk_image.rows());

  // Load the appropriate camera model object and if necessary
  // convert it to the CAHVModel type.
  std::string lcase_file = boost::to_lower_copy(camera_path);
  boost::shared_ptr<vw::camera::CAHVModel> cahv(new vw::camera::CAHVModel);
  if (boost::ends_with(lcase_file, ".cahvore") ) {
    vw::camera::CAHVOREModel cahvore(camera_path);
    *(cahv.get()) = vw::camera::linearize_camera(cahvore, image_size, image_size);
  } else if (boost::ends_with(lcase_file, ".cahvor")  ||
             boost::ends_with(lcase_file, ".cmod"  )   ) {
    vw::camera::CAHVORModel cahvor(camera_path);
    *(cahv.get()) = vw::camera::linearize_camera(cahvor, image_size, image_size);

  } else if ( boost::ends_with(lcase_file, ".cahv") ||
              boost::ends_with(lcase_file, ".pin" )) {
    *(cahv.get()) = vw::camera::CAHVModel(camera_path);

  } else if ( boost::ends_with(lcase_file, ".pinhole") ||
              boost::ends_with(lcase_file, ".tsai"   )   ) {
    // The CAHV class is constructed from a Pinhole model.
    vw::camera::PinholeModel left_pin(camera_path);
    *(cahv.get()) = vw::camera::strip_lens_distortion(left_pin);

  } else {
    vw_throw(vw::ArgumentErr() << "load_cahv_pinhole_camera_model - unsupported camera file type.\n");
  }

  return cahv;
}

int auto_compute_sample_spacing(Vector2i const image_size) {
  int DEFAULT_SPACING = 50;
  int MIN_SIDE_POINTS = 100;
  int MAX_SIDE_POINTS = 600;
  
  // Use the default spacing unless it gives too many points, then use the min spacing
  //  that hits our point estimate.
  int spacing    = DEFAULT_SPACING;
  int side_total = image_size[0] + image_size[1];
  int num_pts    = side_total / DEFAULT_SPACING;
  
  if (num_pts < MIN_SIDE_POINTS)
    spacing = std::max(side_total / MIN_SIDE_POINTS, 1);
  
  if (num_pts > MAX_SIDE_POINTS)
    spacing = side_total / MAX_SIDE_POINTS;
  
  return spacing;
}

void 
resize_epipolar_cameras_to_fit(PinholeModel const& cam1,      PinholeModel const& cam2,
                               PinholeModel      & epi_cam1,  PinholeModel      & epi_cam2,
                               BBox2i       const& roi1,      BBox2i       const& roi2,
                               Vector2i          & epi_size1, Vector2i          & epi_size2) {
  // Get transforms from input images to epipolar images
  CameraTransform<PinholeModel, PinholeModel> in_to_epi1(cam1, epi_cam1);
  CameraTransform<PinholeModel, PinholeModel> in_to_epi2(cam2, epi_cam2);
  
  // Figure out the bbox needed to contain the transformed image
  // - This just uses the ROIs for the two input images so the cameras can be shifted to align with 0,0.
  BBox2 epi_bbox1 = compute_transformed_bbox_fast(roi1, in_to_epi1);
  BBox2 epi_bbox2 = compute_transformed_bbox_fast(roi2, in_to_epi2);

  // Figure out leftmost and uppermost pixel coordinates in resampled images 
  double min_col = std::min(epi_bbox1.min().x(), epi_bbox2.min().x());
  double min_row = std::min(epi_bbox1.min().y(), epi_bbox2.min().y());
  
  // Compute an adjustment of the camera center point (CCD point below focal point)
  //  such that leftmost and uppermost pixels fall at col 0 and row 0 respectively.
  // - Shift the center by the number of pixels converted to physical CCD units.
  // - We can freely adjust the intrinsic portions of the epipolar cameras as long as
  //   we do the same modification to both cameras.
  Vector2 point_offset  = epi_cam1.point_offset();
  Vector2 center_adjust = Vector2(min_col, min_row)*epi_cam1.pixel_pitch();

  // Apply the adjustments to the input epipolar cameras and recompute the bounding boxes  
  epi_cam1.set_point_offset(point_offset - center_adjust);
  epi_cam2.set_point_offset(point_offset - center_adjust);
  CameraTransform<PinholeModel, PinholeModel> in_to_epi1_new(cam1, epi_cam1);
  CameraTransform<PinholeModel, PinholeModel> in_to_epi2_new(cam2, epi_cam2);
  
  // Recompute the bounding boxes to be sure
  epi_bbox1 = compute_transformed_bbox_fast(roi1, in_to_epi1_new);
  epi_bbox2 = compute_transformed_bbox_fast(roi2, in_to_epi2_new);
 
  // Return the size required to contain all of the transformed image data.
  epi_size1 = epi_bbox1.max();
  epi_size2 = epi_bbox2.max();
} // End resize_epipolar_cameras_to_fit


// Convert an optical model to a pinhole model without distortion.
// (The distortion will be taken care of later.)
PinholeModel opticalbar2pinhole(OpticalBarModel const& opb_model, int sample_spacing,
                                double camera_to_ground_dist){

    double pixel_pitch = opb_model.get_pixel_size();
    double  f = opb_model.get_focal_length(); 
    Vector2 c = opb_model.get_optical_center();
    c *= pixel_pitch; // because OpticalBarModel optical center is in pixels

    // Cook up a pinhole model copying the optical center, focal length,
    // pitch, and orientation. This has no distortion, on purpose.
    Vector2i image_size = opb_model.get_image_size();
    PinholeModel out_model = PinholeModel(opb_model.camera_center(image_size/2.0),
                                          opb_model.camera_pose(image_size/2.0).rotation_matrix(),
                                          f, f, c[0], c[1]);
    out_model.set_pixel_pitch(pixel_pitch);

    // This model is very rough. We will improve it by making sure it
    // projects pixels in the camera as close as possible as the
    // original optical bar model.

    // We will keep only the central region of the image with
    // dimensions this number times the original image dimensions, to
    // avoid areas with most distortion.
    double keep_ratio = 1.0; // 1.0/3.0;

    // Since we will keep a smaller area, sample it more.
    double fine_sample_spacing = sample_spacing * keep_ratio; 

    // Grab point pairs for the solver at every
    // interval of "fine_sample_spacing"
    int num_col_samples = image_size[0]/fine_sample_spacing;
    int num_row_samples = image_size[1]/fine_sample_spacing;
    num_col_samples = std::max(num_col_samples, 2);
    num_row_samples = std::max(num_row_samples, 2);
    
    std::vector<Vector3> input_xyz;
    std::vector<double> input_pixels;
    
    // Create pairs of pixels and xyz for the input camera
    for (int c = 0; c <= num_col_samples - 1; c++) {
      // sample the full interval [0, image_size[0]-1] uniformly
      double col = (image_size[0] - 1.0) * double(c)/(num_col_samples - 1.0);
      
      // See the comment earlier where keep_ratio is defined.
      if (std::abs(col - (image_size[0] - 1.0)/2.0) > (image_size[0] - 1.0)*keep_ratio/2.0)
	continue;
      
      for (int r = 0; r <= num_row_samples - 1; r++) {
        // sample the full interval [0, image_size[1]-1] uniformly
	double row = (image_size[1] - 1.0) * double(r)/(num_row_samples - 1.0); 

	if (std::abs(row - (image_size[1] - 1.0)/2.0) > (image_size[1] - 1.0)*keep_ratio/2.0)
	  continue;
      
        Vector2 pix(col, row);
        Vector3 xyz = opb_model.camera_center(pix)
          + camera_to_ground_dist*opb_model.pixel_to_vector(pix);

        input_pixels.push_back(pix[0]);
        input_pixels.push_back(pix[1]);
        input_xyz.push_back(xyz);
      }
      
    }

    // Copy input_pixels to Vector<double> to please the API
    // A value > 0 will make the camera center move less.
    // This not give good results for optical bar, but it may be worth
    // exploring later. Apparently the rays emanating from the optical
    // camera are different than the ones from a pinhole camera
    // (even with distortion), as the optical bar sensor surface
    // is not planar, hence it is not possible to approximate
    // well an optical bar camera with a pinhole cameras everywhere.
    double camera_weight = 0.0; 
    size_t extra = 0;
    if (camera_weight > 0)
      extra = 3;
    Vector<double> reorg_pixels(input_pixels.size() + extra);
    for (size_t it = 0; it < input_pixels.size(); it++) 
      reorg_pixels[it] = input_pixels[it];
    for (size_t it = input_pixels.size(); it < input_pixels.size() + extra; it++)
      reorg_pixels[it] = 0;

    // Refine the camera position and orientation
    double abs_tolerance  = 1e-24;
    double rel_tolerance  = 1e-24;
    int    max_iterations = 2000;
    int status = 0;
    Vector<double> final_params;
    Vector<double> seed;
    CameraSolveLMA<PinholeModel> lma_model(input_xyz, out_model, camera_weight);
    camera_to_vector(out_model, seed);
    final_params = math::levenberg_marquardt(lma_model, seed, reorg_pixels,
                                             status, abs_tolerance, rel_tolerance,
                                             max_iterations);

    vector_to_camera(out_model, final_params);

    return out_model;
  }

/// Find the rotation + translation + scale that best projects given
/// xyz points into given pixels for the input pinhole cameras (each
/// pinhole camera has a few xyz and pixels).
class CameraSolveRotTransScale: public vw::math::LeastSquaresModelBase<CameraSolveRotTransScale> {

  std::vector< std::vector<Vector3> > const& m_xyz;
  vw::Vector<double> const& m_pixel_vec; // for debugging
  std::vector<vw::camera::PinholeModel> m_cameras;
  
public:

  typedef vw::Vector<double>    result_type;   // pixel residuals
  typedef vw::Vector<double, 7> domain_type;   // axis angle + translation + scale
  typedef vw::Matrix<double> jacobian_type;

  /// Instantiate the solver with a set of xyz to pixel pairs and pinhole models
  CameraSolveRotTransScale(std::vector< std::vector<Vector3> > const& xyz,
		  vw::Vector<double> const& pixel_vec,
		  std::vector<vw::camera::PinholeModel> const& cameras):
    m_xyz(xyz), m_pixel_vec(pixel_vec), m_cameras(cameras){
    
    // Sanity check
    if (m_xyz.size() != m_cameras.size()) 
      vw_throw( ArgumentErr() << "Error in CameraSolveRotTransScale: "
		<< "There must be as many xyz sets as cameras.\n");
  }
  
  /// Given the cameras, project xyz into them
  inline result_type operator()(domain_type const& C, bool verbose = false) const {

    // Create the camera models
    std::vector<vw::camera::PinholeModel> cameras = m_cameras; // make a copy local to this function
    for (size_t it = 0; it < cameras.size(); it++) {
      apply_rot_trans_scale(cameras[it], C); // update its parameters
    }
    
    int result_size = 0;
    for (size_t it = 0; it < m_xyz.size(); it++) {
      bool is_good = (m_xyz[it].size() >= 3);
      if (is_good)
        result_size += m_xyz[it].size() * 2;
    }
    
    result_type result;
    result.set_size(result_size);
    int count = 0;
    for (size_t it = 0; it < m_xyz.size(); it++) {
      
      bool is_good = (m_xyz[it].size() >= 3);
      if (is_good) {
        for (size_t c = 0; c < m_xyz[it].size(); c++) {
          Vector2 pixel = cameras[it].point_to_pixel(m_xyz[it][c]);
          result[2*count  ] = pixel[0];
          result[2*count+1] = pixel[1];
          count++;
        }
      }
    }
    
    if (2*count != result_size) {
      vw_throw( LogicErr() << "Book-keeping failure in CameraSolveRotTransScale.\n");
    }

    if (verbose) {
      vw_out() << "Pixels and pixel errors after optimizing the transform to the ground.";
    }
    double cost = 0;
    for (int it = 0; it < result_size; it++) {
      double diff = result[it] - m_pixel_vec[it];
      cost += diff*diff;
      if (verbose)
        vw_out() << result[it] << ' ' << m_pixel_vec[it] << ' '
                 << std::abs(result[it] - m_pixel_vec[it]) << std::endl;
    }

    return result;
  }
  
}; // End class CameraSolveRotTransScale

// Given original cams in sfm_cams and individually scaled cameras in
// aux_cams, get the median scale change from the first set to the second one.
// It is important to do the median, since scaling the cameras individually
// is a bit of a shaky business.
double find_median_scale_change(std::vector<vw::camera::PinholeModel> const & sfm_cams,
				std::vector<vw::camera::PinholeModel> const & aux_cams,
				std::vector< std::vector<vw::Vector3>> const& xyz){
  
  int num_cams = sfm_cams.size();

  std::vector<double> scales;
  
  for (int it1 = 0; it1 < num_cams; it1++) {

    bool is_good = (xyz[it1].size() >= 3);
    if (!is_good)
      continue;
    
    for (int it2 = it1 + 1; it2 < num_cams; it2++) {
      
      bool is_good = (xyz[it2].size() >= 3);
      if (!is_good)
        continue;
    
      double len1 = norm_2(sfm_cams[it1].camera_center()
			   - sfm_cams[it2].camera_center());
      double len2 = norm_2(aux_cams[it1].camera_center()
			   - aux_cams[it2].camera_center());
      
      double scale = len2/len1;
      scales.push_back(scale);
    }
  }

  if (scales.empty())
    vw_throw( LogicErr() << "Could not find two images with at least 3 GCP each.\n");
    
  double median_scale = vw::math::destructive_median(scales);

  return median_scale;
}

// Given some GCP so that at least two images have at at least three GCP each,
// but each GCP is allowed to show in one image only, use the GCP
// to transform cameras to ground coordinates.
void align_cameras_to_ground(std::vector<std::vector<Vector3>> const& xyz,
                             std::vector<std::vector<Vector2>> const& pix,
                             std::vector<PinholeModel> & sfm_cams,
                             Matrix3x3 & rotation, 
                             Vector3 & translation,
                             double & scale) {
  
  std::string camera_type = "pinhole";
  bool refine_camera = true;
  bool verbose = false; 

  // Cameras individually aligned to ground using GCP. They may not be
  // self-consistent, and are only used to give an idea of the
  // transform to apply to the unaligned cameras.
  std::vector<PinholeModel> aux_cams;

  int num_cams = sfm_cams.size();
  for (int it = 0; it < num_cams; it++) {
    // Export to the format used by the API
    std::vector<double> pixel_values;
    for (size_t c = 0; c < pix[it].size(); c++) {
      pixel_values.push_back(pix[it][c][0]);
      pixel_values.push_back(pix[it][c][1]);
    }

    vw::CamPtr out_cam(new PinholeModel(sfm_cams[it]));

    bool is_good = (xyz[it].size() >= 3);
    if (is_good) 
      fit_camera_to_xyz(camera_type, refine_camera,  
			xyz[it], pixel_values, verbose, out_cam);
    
    aux_cams.push_back(*((PinholeModel*)out_cam.get()));
  }

  double world_scale = find_median_scale_change(sfm_cams, aux_cams, xyz);
  vw_out() << "Initial guess scale to apply when converting to world coordinates using GCP: "
	   << world_scale << ".\n";

  // So far we aligned both cameras individually to GCP and we got an
  // idea of scale.  Yet we would like to align them without changing
  // the relationship between them, so using a single transform for
  // all not an individual transform for each.  This way we will
  // transform the SfM-computed cameras to the new coordinate system.

  // Start by estimating a such a transform.
  int num_pts = 0;
  for (int it = 0; it < num_cams; it++) {
    bool is_good = (xyz[it].size() >= 3);
    if (is_good) 
      num_pts += pix[it].size();
  }
  
  vw::Matrix<double> in_pts, out_pts;
  in_pts.set_size(3, num_pts);
  out_pts.set_size(3, num_pts);
  
  int col = 0;
  for (int it = 0; it < num_cams; it++) {
    
    bool is_good = (xyz[it].size() >= 3);
    if (is_good) {
      // For each camera, find xyz values in the input cameras
      // that map to GCP. Use the scale for that.
      for (int c = 0; c < xyz[it].size(); c++) {
        
        // Distance from camera center to xyz for the individually aligned cameras
        double len = norm_2(aux_cams[it].camera_center() - xyz[it][c]);
        len = len / world_scale;
        Vector3 trans_xyz = sfm_cams[it].camera_center()
          + len * sfm_cams[it].pixel_to_vector(pix[it][c]);
        for (int row = 0; row < in_pts.rows(); row++) {
          in_pts(row, col)  = trans_xyz[row];
          out_pts(row, col) = xyz[it][c][row];
        }
        col++;
      }
    }
  }
  
  if (col != num_pts) 
    vw_throw( LogicErr() << "Book-keeping failure in aligning cameras to ground.\n");

  // The initial transform to world coordinates
  Vector<double> C;
  vw::math::find_3D_transform(in_pts, out_pts, rotation, translation, scale);

  // Copy into C
  transform_to_vector(C, rotation, translation, scale);

  // Form the pixel vector
  int pixel_vec_len = 0;
  for (size_t it = 0; it < pix.size(); it++) {
    bool is_good = (xyz[it].size() >= 3);
    if (is_good)
      pixel_vec_len += pix[it].size() * 2;
  }
  Vector<double> pixel_vec;
  pixel_vec.set_size(pixel_vec_len);
  int count = 0;
  for (size_t it = 0; it < pix.size(); it++) {
    bool is_good = (xyz[it].size() >= 3);
    if (is_good) {
      for (size_t c = 0; c < pix[it].size(); c++) {
        Vector2 pixel = pix[it][c];
        pixel_vec[2*count  ] = pixel[0];
        pixel_vec[2*count+1] = pixel[1];
        count++;
      }
    }
  }
  if (2*count != pixel_vec_len)
    vw_throw( LogicErr() << "Book-keeping failure in cam_gen.\n");
  
  // Optimize the transform
  double abs_tolerance  = 1e-24;
  double rel_tolerance  = 1e-24;
  int    max_iterations = 2000;
  int status = 0;
  CameraSolveRotTransScale lma_model(xyz, pixel_vec, sfm_cams);
  Vector<double> final_params
    = vw::math::levenberg_marquardt(lma_model, C, pixel_vec,
				    status, abs_tolerance, rel_tolerance,
				    max_iterations);

  Vector<double> final_residual = lma_model(final_params, verbose);
  
  // Bring the cameras to world coordinates
  for (int it = 0; it < num_cams; it++) 
    apply_rot_trans_scale(sfm_cams[it], final_params);

  // Unpack the final vector into a rotation + translation + scale
  vector_to_transform(final_params, rotation, translation, scale);

}

PinholeModel fitPinholeModel(CameraModel const* in_model, 
                             vw::Vector2 const& image_size,
                             std::string const& out_distortion_type,
                             bool force_conversion,
                             int sample_spacing,
                             int rpc_degree,
                             double camera_to_ground_dist) {

  PinholeModel out_model;

  if (out_distortion_type == "TsaiLensDistortion") {
    create_approx_pinhole_model<TsaiLensDistortion>
      (in_model, out_model, image_size, sample_spacing, force_conversion,
        rpc_degree, camera_to_ground_dist);
  }else if (out_distortion_type == "BrownConradyDistortion") {
    create_approx_pinhole_model<BrownConradyDistortion>
      (in_model, out_model, image_size, sample_spacing, force_conversion,
        rpc_degree, camera_to_ground_dist);
  } else if (out_distortion_type == RPCLensDistortion::class_name()) {
    create_approx_pinhole_model<RPCLensDistortion>
      (in_model, out_model, image_size, sample_spacing, force_conversion,
        rpc_degree, camera_to_ground_dist);
  }else{
    vw::vw_throw(vw::ArgumentErr() 
                   << "Unsupported output model type: " << out_distortion_type << "\n");
  }
    
  return out_model;
}
  
}} // end namespace vw::camera

