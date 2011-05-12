// __BEGIN_LICENSE__
// Copyright (C) 2006-2011 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// \file CAHVOREModel.h
///
/// This file contains the CAHVORE camera model.  This camera
/// model is a refinement of the \ref CAHVORModel: it adds extra terms
/// to model fisheye lenses
///
#ifndef __VW_CAMERAMODEL_CAHVORE_H__
#define __VW_CAMERAMODEL_CAHVORE_H__

#include <vw/Core/Log.h>
#include <vw/Camera/CAHVModel.h>

namespace vw {
namespace camera {

  class CAHVOREModel : public CameraModel {
  public:
    typedef CAHVModel linearized_type;

    //------------------------------------------------------------------
    // Constructors / Destructors
    //------------------------------------------------------------------

    /// Default constructor creates CAHVORE vectors equal to 0.
    CAHVOREModel() {}

    /// Read a CAHVORE file from disk.
    CAHVOREModel(std::string const& filename);

    /// Initialize the CAHVORE vectors directly in the native CAHVORE format.
    CAHVOREModel(Vector3 const& C_vec, Vector3 const& A_vec, Vector3 const& H_vec, Vector3 const& V_vec,
                 Vector3 const& O_vec, Vector3 const& R_vec, Vector3 const& E_Vec) :
      C(C_vec), A(A_vec), H(H_vec), V(V_vec), O(O_vec), R(R_vec), E(E_Vec), P(1.0) {}
    CAHVOREModel(Vector3 const& C_vec, Vector3 const& A_vec, Vector3 const& H_vec, Vector3 const& V_vec,
                 Vector3 const& O_vec, Vector3 const& R_vec, Vector3 const& E_Vec, int T, double P_v) :
      C(C_vec), A(A_vec), H(H_vec), V(V_vec), O(O_vec), R(R_vec), E(E_Vec), P(P_v) {
      switch( T ) {
      case 1: P = 1.0; break;
      case 2: P = 0.0; break;
      case 3:
        if ( P < 0 || P > 1 ) vw_throw( ArgumentErr() << "Invalid P value: " << P_v << "\n" );
        break;
      default: vw_throw( ArgumentErr() << "Unknown CAHVORE type: " << T << "\n" );
      }
    }
    CAHVOREModel(Vector3 const& C_vec, Vector3 const& A_vec, Vector3 const& H_vec, Vector3 const& V_vec,
                 Vector3 const& O_vec, Vector3 const& R_vec, Vector3 const& E_Vec, double P_v) :
      C(C_vec), A(A_vec), H(H_vec), V(V_vec), O(O_vec), R(R_vec), E(E_Vec), P(P_v) {
      if ( P < 0 || P > 1 ) vw_throw( ArgumentErr() << "Invalid P value: " << P_v << "\n" );
    }

    virtual ~CAHVOREModel() {}
    virtual std::string type() const { return "CAHVORE"; }

    //------------------------------------------------------------------
    // Methods
    //------------------------------------------------------------------
    virtual Vector2 point_to_pixel(Vector3 const& point) const;
    virtual Vector3 pixel_to_vector(Vector2 const& pix) const;
    virtual Vector3 camera_center(Vector2 const& /*pix*/ = Vector2() ) const { return C; };

    /// Write CAHVORE model to file.
    void write(std::string const& filename);

    //------------------------------------------------------------------
    // Exposed Variables
    //------------------------------------------------------------------
    Vector3   C;
    Vector3   A;
    Vector3   H;
    Vector3   V;
    Vector3   O;
    Vector3   R;
    Vector3   E;
    double    P; // We don't have T as it is redundant information
  private:
    bool check_line( std::istream& istream, char letter );
  };

  // Function to "map" the CAHVORE parameters into CAHV:
  CAHVModel linearize_camera( CAHVOREModel const& camera_model,
                              Vector2i const& cahvore_image_size,
                              Vector2i const& cahv_image_size );

}} // namespace vw::camera

#endif//__VW_CAMERAMODEL_CAHVORE_H__
