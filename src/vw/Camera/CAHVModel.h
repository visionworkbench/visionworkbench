// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// \file CAHVModel.h
/// 
/// This file contains the CAHV pinhole camera model.
/// 
#ifndef __VW_CAMERAMODEL_CAHV_H__
#define __VW_CAMERAMODEL_CAHV_H__

#include <vw/Camera/CameraModel.h>
#include <vw/Camera/PinholeModel.h>
#include <vw/Math/Vector.h>
#include <vw/Math/Matrix.h>

namespace vw { 
namespace camera {


  /// The CAHV pinhole camera model has been widely used in NASA
  /// planetary mission for rover navigation and scientific camera
  /// systems.  In the CAHV model, camera intrinsics and extrincsics
  /// are jointly parameterized by four 3-Vectors C,A,H, and V.  The
  /// strength of this model is the mathematical ease with which 3D
  /// points can be projected onto the image plane.  The CAHV model is
  /// computationally efficient, however it does not model lens
  /// distortion.  For this, its close relative the CAHVOR model
  /// should be used.  Refer to \ref vw::camera::CAHVORModel for more
  /// information.
  ///
  /// References: 
  ///
  /// Yakimovsky, Y. and Cunningham R., "A System for Extracting
  /// Three-Dimensional Measurements from a Stereo Pair of TV
  /// Cameras". Computer Graphics and Image Processing 7,
  /// pp. 195-210. (1978)
  ///
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
    
    /// Initialize the CAHV vectors directly in the native CAHV format.
    CAHVModel(Vector3 C_vec, Vector3 A_vec, Vector3 H_vec, Vector3 V_vec) : 
      C(C_vec), A(A_vec), H(H_vec), V(V_vec) {}

    /// Initialize a CAHV model from a pinhole model 
    CAHVModel(PinholeModel const& pin_model) {
      //  Intrinsic parametes (in pixel units)
      double fH, fV, Hc, Vc;
      pin_model.intrinsic_parameters(fH, fV, Hc, Vc);    

      Vector3 u,v,w;
      pin_model.coordinate_frame(u,v,w);

      Matrix<double,3,3> R = pin_model.camera_pose().rotation_matrix();

      Vector3 Hvec = R*u;
      Vector3 Vvec = R*v;

      C = pin_model.camera_center();
      A = R*w;
      H = fH*Hvec + Hc*A;
      V = fV*Vvec + Vc*A;	      
    }
    
    
    CAHVModel operator= (PinholeModel const& pin_model);
    
    virtual std::string type() const { return "CAHV"; }
    
    /// Initialize the CAHV vectors indirectly using pinhole camera
    /// parameters.  In this variant, the view matrix is supplied
    /// directly.
    ///
    /// f         - focal length in world units
    /// pixel_size - (width, height) of a pixel on the image plane in world units
    /// xMin, xMax, yMin, yMax - image plane boundaries in units of pixels
    ///
    /// View matrix transformation:
    /// 
    /// |R00 R01 R02 Tx|
    /// |R10 R11 R12 Ty|
    /// |R20 R21 R22 Tz|
    /// | 0   0   0   1|
    /// 
    CAHVModel(double f, Vector2 pixel_size, 
              double xmin, double xmax, double ymin, double ymax,
              Matrix<double, 4, 4> view_matrix) {

      double fH = f/pixel_size.x();
      double fV = f/pixel_size.y();
      double Hc = -xmin;
      double Vc = -ymin;
      
      Vector3 Hvec(view_matrix[0][0], view_matrix[0][1], view_matrix[0][2]);
      Vector3 Vvec(view_matrix[1][0], view_matrix[1][1], view_matrix[1][2]);
      
      C = Vector3(view_matrix[3][0], view_matrix[3][1], view_matrix[3][2]);
      A = Vector3(view_matrix[2][0], view_matrix[2][1], view_matrix[2][2]);
      H = fH*Hvec + Hc*A;
      V = fV*Vvec + Vc*A;	

    }
    
    /// Initialize the CAHV vectors indirectly using pinhole camera
    /// parameters:
    ///
    /// f         - focal length in world units
    /// pixel_size - (width, height) of a pixel on the image plane in world units
    /// Hc, Vc    - location of projection of C on the image plane in units of
    ///             pixels from lower left hand corner
    ///
    /// C         - center of projection in world coordinates
    /// A         - unit vector perpendicular to image plane in world coords
    ///             pointing out the "front" of the camera
    /// Hvec      - unit horizontal image plane vector in world coords
    /// Vvec      - unit vertical image plane vector in world coords
    /// 
    CAHVModel(double f, Vector2 pixel_size, double Hc, double Vc,
              Vector3 Cinit, Vector3 Ainit, Vector3 Hvec, Vector3 Vvec) {
      double fH = f/pixel_size.x();
      double fV = f/pixel_size.y();
      C = Cinit;
      A = Ainit;
      H = fH*Hvec + Hc*A;
      V = fV*Vvec + Vc*A;
    }

    virtual ~CAHVModel() {}
    
    //------------------------------------------------------------------
    // Methods
    //------------------------------------------------------------------
    virtual Vector2 point_to_pixel(Vector3 const& point) const;
    virtual Vector3 pixel_to_vector (Vector2 const& pix) const;
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

  /// Given two CAHV camera models, this method returns two new camera
  /// models that have been epipolar rectified.
  void epipolar(CAHVModel const src_camera0, CAHVModel const src_camera1, 
                CAHVModel &dst_camera0, CAHVModel &dst_camera1);

  inline std::ostream& operator<<(std::ostream& str, CAHVModel const& model) {
    str << "CAHV camera: \n";
    str << "\tC: " << model.C << "\n";
    str << "\tA: " << model.A << "\n";
    str << "\tH: " << model.H << "\n";
    str << "\tV: " << model.V << "\n";
    return str;
  }

  
}}	// namespace vw::camera

#endif	//__CAMERAMODEL_CAHV_H__
