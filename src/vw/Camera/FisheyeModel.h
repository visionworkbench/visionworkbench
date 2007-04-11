// __BEGIN_LICENSE__
// 
// Copyright (C) 2006 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration
// (NASA).  All Rights Reserved.
// 
// Copyright 2006 Carnegie Mellon University. All rights reserved.
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

/// \file FisheyeModel.h
/// 
/// This file contains the fisheye camera model.
///
/// This camera model supports several different types of fisheye lens
/// models.  In the definitions that follow, r is the distance to the
/// image center, f is the focal length, and theta is the angle
/// between the center of the image and the point r.
///
/// Gnomonical: r = f * tan(theta) (no distortion -- same as the pinhole camera)
/// Equidistant: r = f * theta 
/// Orthographic: r = f * sin(theta)
/// Equal Area: r = 2fsin(theta / 2)
/// Stereographic: r = 2ftan(theta / 2)
///
#ifndef __VW_CAMERAMODEL_FISHEYE_H__
#define __VW_CAMERAMODEL_FISHEYE_H__

#include <vw/Camera/CameraModel.h>
#include <vw/Image/ImageView.h>
#include <vw/Image/PixelTypes.h>
#include <vw/Math/Vector.h>
#include <vw/Math/Matrix.h>

namespace vw { 
namespace camera {

  namespace fisheye {
    struct GnomonicalProjection {
      double operator() (double f, double theta) const { return f * tan(theta); }
    };
    
    /// Preserves distance relationships
    struct EquidistantProjection {
      double operator() (double f, double theta) const { return f * theta; }
    };
    
    struct OrthographicProjection {
      double operator() (double f, double theta) const { return f * sin(theta); }      
    };

    /// Preserves area relationships
    struct EqualAreaProjection {
      double operator() (double f, double theta) const { return 2*f * sin(theta/2); }      
    };

    struct StereographicProjection {
      double operator() (double f, double theta) const { return 2*f * tan(theta/2); }      
    };
  } // namespace fisheye    

  template <class ProjectionT = fisheye::EquidistantProjection>
  class FisheyeModel : public CameraModel {
  public:

    //------------------------------------------------------------------
    // Constructors / Destructors
    //------------------------------------------------------------------
    
    /// Create a new fisheye model
    /// 
    /// Here, field of view is measured from the top left corner of
    /// the image to the bottom right corner.
    FisheyeModel(double focal_length, double field_of_view,
                 int32 image_width, int32 image_height,
                 Vector2 principal_point, ProjectionT projection = ProjectionT()) :
      m_focal_length(focal_length), m_field_of_view(field_of_view), 
      m_width(image_width), m_height(image_height),
      m_principal_point(principal_point), m_projection(projection) {
      
      // Set a default set of extrinsics
      set_extrinsics(Vector3(0,0,0),
                     Quaternion<double>(1,0,0,0));
      
      // If the field of view is measured from the left center of the
      // image to the right center of the image, the pixel size is:
      //      m_pixel_size = field_of_view / (2 * width);

      // If the field of view is measured from the upper left corner
      // of the image to the lower right corner of the image, the
      // pixel size is:
      m_pixel_size = sqrt((field_of_view * field_of_view)/2);
    }

    virtual ~FisheyeModel() {}
    
    //------------------------------------------------------------------
    // Methods
    //------------------------------------------------------------------
    virtual Vector2 point_to_pixel(Vector3 const& point) const {
      vw_throw( vw::NoImplErr() << "FisheyeModel: PointToPixel not yet implemented." );
    }

    virtual Vector3 pixel_to_vector (Vector2 const& pix) const {
      double theta_x = (pix.x() - m_camera_center.x()) * m_pixel_size;
      double theta_y = (pix.y() - m_camera_center.y()) * m_pixel_size;
      double r_x = m_projection(m_focal_length, theta_x);
      double r_y = m_projection(m_focal_length, theta_y);
      return norm_2(Vector3(r_x, r_y, m_focal_length));
    }

    virtual Vector3 camera_center(Vector2 const& pix = Vector2() ) const { return m_camera_center; };
   
    void set_extrinsics(Vector3 camera_center, Quaternion<double> pose) {
      m_camera_center = camera_center;
      m_pose = pose;
    }


  private:
    double m_focal_length;
    double m_field_of_view;
    double m_pixel_size;
    Vector3 m_camera_center;
    Quaternion<double> m_pose;
    Vector2 m_principal_point;
    int m_width, m_height;
    ProjectionT m_projection;
  };
  
}}	// namespace vw::camera

#endif	//__CAMERAMODEL_CAHV_H__
