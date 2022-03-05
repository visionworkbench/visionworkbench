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

// TODO: Move this to RPCLensDistortion() as only that class needs it.
void update_rpc_undistortion(PinholeModel const& model){

  const vw::camera::LensDistortion* distortion = model.lens_distortion();
  std::string lens_name = distortion->name();
  if (lens_name != RPCLensDistortion::class_name())
    return;
  
  // Have to cast away the const-ness. Not nice. Only happens for RPC
  // distortion.
  PinholeModel * pin_ptr = const_cast<PinholeModel*>(&model);

  if (lens_name == RPCLensDistortion::class_name()) {
    RPCLensDistortion * rpc_dist = dynamic_cast<RPCLensDistortion*>
      (const_cast<LensDistortion*>(distortion));
    if (rpc_dist == NULL) 
      vw_throw( ArgumentErr() << "PinholeModel::expecting an " + RPCLensDistortion::class_name() +
                " model." );
    
    // Only update this if we have to
    int sample_spacing = auto_compute_sample_spacing(rpc_dist->image_size());
    int rpc_degree = rpc_dist->rpc_degree();
    if (!rpc_dist->can_undistort()) 
      compute_undistortion<RPCLensDistortion>(*pin_ptr, rpc_dist->image_size(), sample_spacing,
                                              rpc_degree);
  }

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
  
}} // end namespace vw::camera

