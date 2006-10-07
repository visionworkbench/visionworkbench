// __BEGIN_LICENSE__
//
// Copyright (C) 2006 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration
// (NASA).  All Rights Reserved.
// 
// This software is distributed under the NASA Open Source Agreement
// (NOSA), version 1.3.  The NOSA has been approved by the Open Source
// Initiative.  See the file COPYING at the top of the distribution
// directory tree for the complete NOSA document.
// 
// THE SUBJECT SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY WARRANTY OF ANY
// KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT
// LIMITED TO, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL CONFORM TO
// SPECIFICATIONS, ANY IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR
// A PARTICULAR PURPOSE, OR FREEDOM FROM INFRINGEMENT, ANY WARRANTY THAT
// THE SUBJECT SOFTWARE WILL BE ERROR FREE, OR ANY WARRANTY THAT
// DOCUMENTATION, IF PROVIDED, WILL CONFORM TO THE SUBJECT SOFTWARE.
//
// __END_LICENSE__
#ifndef __VW_CAMERA_EXIF_H__
#define __VW_CAMERA_EXIF_H__

#include <vw/Core/Exception.h>
#include <vw/Camera/ExifData.h>

namespace vw { 
namespace camera {

  VW_DEFINE_EXCEPTION(ExifErr, vw::Exception);

  class ExifView {
  private:
    ExifData data;
    
  public:
    ExifView() {}
    
    bool load_exif(std::string const& filename);
    
    // Camera info
    std::string get_make();
    std::string get_model();
    
    // Camera settings
    double get_f_number();
    double get_exposure_time();
    int get_iso();
    
    // APEX equivalents
    // TODO: report in some logical way when value doesn't exist
    double get_aperture_value();
    double get_shutter_speed_value();
    double get_film_sensitivity();
    double get_brightness_value();
  };

}} // namespace vw::camera
  
#endif  // __VW_CAMERA_EXIF_H__
