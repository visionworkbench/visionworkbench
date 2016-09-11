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


/// \file PinholeModel.h
///
/// This file contains the pinhole camera model.
///
#ifndef __VW_CAMERAMODEL_CAMERAUTILITIES_H__
#define __VW_CAMERAMODEL_CAMERAUTILITIES_H__

#include <vw/Camera/PinholeModel.h>
#include <vw/Camera/CAHVModel.h>
#include <vw/Camera/CAHVORModel.h>
#include <vw/Camera/CAHVOREModel.h>
#include <vw/Camera/CameraTransform.h>
#include <vw/FileIO/DiskImageView.h>
#include <boost/shared_ptr.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/algorithm/string/case_conv.hpp>
#include <boost/algorithm/string/predicate.hpp>

namespace vw {
namespace camera {


/// Load a pinhole camera model of any supported type
inline boost::shared_ptr<vw::camera::CameraModel>
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
inline boost::shared_ptr<vw::camera::CAHVModel>
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
    vw::camera::PinholeModel left_pin(camera_path);
    *(cahv.get()) = vw::camera::strip_lens_distortion(left_pin);

  } else {
    vw_throw(vw::ArgumentErr() << "load_cahv_pinhole_camera_model - unsupported camera file type.\n");
  }

  return cahv;
}

// Given two epipolar rectified CAHV camera models and input images,
// get the aligned version of the images suitable for writing to disk and processing.
template <class ImageInT, class ImageOutT>
void get_epipolar_transformed_images(std::string const& left_camera_file,
                                     std::string const& right_camera_file,
                                     boost::shared_ptr<camera::CameraModel> const left_cahv_camera,
                                     boost::shared_ptr<camera::CameraModel> const right_cahv_camera,
                                     ImageInT  const& left_image_in,
                                     ImageInT  const& right_image_in,
                                     ImageOutT      & left_image_out,
                                     ImageOutT      & right_image_out) {

    // In the epipolar alignment case, the "camera_models" function returns the CAHVModel type!
    CAHVModel* left_epipolar_cahv  = dynamic_cast<CAHVModel*>(vw::camera::unadjusted_model(&(*left_cahv_camera )));
    CAHVModel* right_epipolar_cahv = dynamic_cast<CAHVModel*>(vw::camera::unadjusted_model(&(*right_cahv_camera)));
    if (!left_epipolar_cahv || !right_epipolar_cahv) {
      vw_throw(ArgumentErr() << "load_cahv_pinhole_camera_model: CAHVModel cast failed!\n");
    }

    std::string lcase_file = boost::to_lower_copy(left_camera_file);

    // Remove lens distortion and create epipolar rectified images.
    if (boost::ends_with(lcase_file, ".cahvore")) {
      CAHVOREModel left_cahvore (left_camera_file );
      CAHVOREModel right_cahvore(right_camera_file);
      left_image_out  = camera_transform(left_image_in,  left_cahvore,  *left_epipolar_cahv );
      right_image_out = camera_transform(right_image_in, right_cahvore, *right_epipolar_cahv);
      
    } else if (boost::ends_with(lcase_file, ".cahvor") ||
               boost::ends_with(lcase_file, ".cmod") ) {
      CAHVORModel left_cahvor (left_camera_file );
      CAHVORModel right_cahvor(right_camera_file);
      left_image_out  = camera_transform(left_image_in,  left_cahvor,  *left_epipolar_cahv );
      right_image_out = camera_transform(right_image_in, right_cahvor, *right_epipolar_cahv);

    } else if ( boost::ends_with(lcase_file, ".cahv") ||
                boost::ends_with(lcase_file, ".pin" )) {
      CAHVModel left_cahv (left_camera_file );
      CAHVModel right_cahv(right_camera_file);
      left_image_out  = camera_transform(left_image_in,  left_cahv,  *left_epipolar_cahv );
      right_image_out = camera_transform(right_image_in, right_cahv, *right_epipolar_cahv);

    } else if ( boost::ends_with(lcase_file, ".pinhole") ||
                boost::ends_with(lcase_file, ".tsai") ) {
      PinholeModel left_pin (left_camera_file );
      PinholeModel right_pin(right_camera_file);
      
      left_image_out  = camera_transform(left_image_in,  left_pin,  *left_epipolar_cahv );
      right_image_out = camera_transform(right_image_in, right_pin, *right_epipolar_cahv);

    } else {
      vw_throw(ArgumentErr() << "load_cahv_pinhole_camera_model: unsupported camera file type.\n");
    }

}



}}      // namespace vw::camera

#endif  //__CAMERAMODEL_CAMERAUTILITIES_H__
