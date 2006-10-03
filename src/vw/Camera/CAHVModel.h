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

/// \file CAHV.h
/// 
/// A CAHV Camera Model.
/// 
#ifndef __VW_CAMERAMODEL_CAHV_H__
#define __VW_CAMERAMODEL_CAHV_H__

#include <vw/Camera/CameraModel.h>
#include <vw/Image/ImageView.h>
#include <vw/Image/PixelTypes.h>
#include <vw/Math/Vector.h>
#include <vw/Math/Matrix.h>

namespace vw { 
namespace camera {

  class CAHVModel : public CameraModel {
  public:
    //------------------------------------------------------------------
    // Constructors / Destructors
    //------------------------------------------------------------------
    
    /// Initialize an empty camera model where the 4 CAHV vectors are
    /// set to zero.
    CAHVModel() {}  

    /// This constructor takes a filename and reads in a camera model
    /// from the file.  The file may contain either CAHV parameters or
    /// pinhole camera parameters.
    CAHVModel(std::string const& filename);
    
    /// Inialize the CAHV vectors directly.
    CAHVModel(Vector3 C_init, Vector3 A_init,
              Vector3 H_init, Vector3 V_init);
    
    /// Initialize the CAHV vectors indirectly using pinhole camera
    /// parameters.  In this variant, the view matrix is supplied
    /// directly.
    CAHVModel(double f, Vector2 pixelSize,
              double xmin, double xmax, double ymin, double ymax,
              Matrix<double, 4, 4> viewMatrix);
    
    /// Initialize the CAHV vectors indirectly using pinhole camera
    /// parameters.
    CAHVModel(double f, Vector2  pixelSize, double Hc, double Vc,
              Vector3 Cinit, Vector3 Ainit,
              Vector3 Hvec, Vector3 Vvec);

    virtual ~CAHVModel() {}
    
    //------------------------------------------------------------------
    // Methods
    //------------------------------------------------------------------
    virtual Vector3 pixel_to_vector(Vector2 const& pix) const;
    virtual Vector2 vector_to_pixel(Vector3 const& vec) const;
    virtual Vector3 camera_center(Vector2 const& pix = Vector2() ) const { return C; };
   
    //------------------------------------------------------------------
    // Exposed Variables
    //------------------------------------------------------------------
    Vector3   C;
    Vector3   A;
    Vector3   H;
    Vector3   V;

  private:
    void read_cahv   (std::string const& filename);
    void read_pinhole(std::string const& filename);
  };

  void epipolar(CAHVModel &src_camera0, CAHVModel &src_camera1, 
                CAHVModel &dst_camera0, CAHVModel &dst_camera1);

  
}}	// namespace vw::camera

#endif	//__CAMERAMODEL_CAHV_H__
