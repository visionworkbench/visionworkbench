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
// 
#include <vw/Camera/CAHVModel.h>
#include <boost/algorithm/string.hpp>

using namespace std;

namespace vw {
namespace camera {

  // FIXME -- Disabled for now until Pinhole API stabilizes. -mbroxton
  // Pinhole may have stabilized... make sure to incorporate the u,v,w
  // camera frame vectors.
  CAHVModel CAHVModel::operator= (PinholeModel const& pin_model) {

    //  Pinhole model parameters (in pixel units)
    double fH, fV, Hc, Vc;
    pin_model.intrinsic_parameters(fH, fV, Hc, Vc);

    //  Unit vectors defining camera coordinate frame
    Vector3 u,v,w;
    pin_model.get_coordinate_frame(u,v,w);

    //  The true rotation between world and camera coordinate
    //  frames includes the rotation R --AND-- a rotation from
    //  specifying the directions of increasing u,v,w pixels
    Matrix<double,3,3> R = pin_model.camera_pose().rotation_matrix();
    Matrix<double,3,3> R_uvw;
    select_row(R_uvw,0) = u;
    select_row(R_uvw,1) = v;
    select_row(R_uvw,2) = w;
    Matrix<double,3,3> rot_matrix = R_uvw * R;

    //  Now create the components of the CAHV model...
    Vector3 Hvec = select_row(rot_matrix,0);
    Vector3 Vvec = select_row(rot_matrix,1);
      
    C = pin_model.camera_center();
    A = select_col(transpose(rot_matrix), 2);
    H = fH*Hvec + Hc*A;
    V = fV*Vvec + Vc*A;	      

    return *this;
  }


  /// This constructor takes a filename and reads in a camera model
  /// from the file.  The file may contain either CAHV parameters or
  /// pinhole camera parameters.
  CAHVModel::CAHVModel(std::string const& filename) {
    if (filename.size() == 0)
      vw_throw( IOErr() << "CAHVModel: null file name passed to constructor." );
    
    if (boost::ends_with(filename, ".cahv"))
      read_cahv(filename);
    else if (boost::ends_with(filename, ".pin"))
      read_pinhole(filename);
    else 
      vw_throw( IOErr() << "CAHVModel: Unknown camera file suffix." );
  }

  Vector2 CAHVModel::point_to_pixel(Vector3 const& point) const {
    double dDot = dot_prod(point-C, A);
    return Vector2( dot_prod(point-C, H) / dDot,
                    dot_prod(point-C, V) / dDot );
  }

  Vector3 CAHVModel::pixel_to_vector(Vector2 const& pix) const {
    Vector3 va, vb;
    
    // Vertical component in a plane perpendicular to the vector
    va = V + (-(pix.y()) * A);
    
    // Horizontal component in a plane perpendicular to the vector
    vb = H + (-(pix.x()) * A);
    
    // Find vector
    Vector3 vec = cross_prod(va, vb);
    
    // Normalize vector
    vec *= 1.0 / norm_2(vec);
    
    // The vector VxH should be pointing in the same directions as A,
    // if it isn't (because we have a left handed system), flip the
    // vector.
    Vector3 temp = cross_prod(V, H);
    if (dot_prod(temp, A) < 0.0){
      vec *= -1.0;
      //cout << "CAHV PixelToVector changed sign of vec" << endl;
    }
    return vec;
  }

  // --------------------------------------------------
  //                 Private Methods
  // --------------------------------------------------
  void CAHVModel::read_cahv(std::string const& filename) {

    FILE *cahvFP = fopen(filename.c_str(), "r");
    if (cahvFP == 0)
      vw_throw( IOErr() << "CAHVModel::read_cahv: Could not open file\n" );
    
    char line[4096];
    
    // Scan through comments
    fgets(line, sizeof(line), cahvFP);
    while(line[0] == '#') 
      fgets(line, sizeof(line), cahvFP);
    
    if (sscanf(line,"C = %lf %lf %lf", &C(0), &C(1), &C(2)) != 3) {
      vw_throw( IOErr() << "CAHVModel::read_CAHV: Could not read C vector\n" );
      fclose(cahvFP);
    }

    fgets(line, sizeof(line), cahvFP);
    if (sscanf(line,"A = %lf %lf %lf", &A(0),&A(1), &A(2)) != 3) {
      vw_throw( IOErr() << "CAHVModel::read_CAHV: Could not read A vector\n" );
      fclose(cahvFP);
    }

    fgets(line, sizeof(line), cahvFP);
    if (sscanf(line,"H = %lf %lf %lf", &H(0), &H(1), &H(2)) != 3) {
      vw_throw( IOErr() << "CAHVModel::read_CAHV: Could not read H vector\n" );
      fclose(cahvFP);
    }

    fgets(line, sizeof(line), cahvFP);
    if (sscanf(line,"V = %lf %lf %lf", &V(0), &V(1), &V(2)) != 3) {
      vw_throw( IOErr() << "CAHVModel::read_CAHV: Could not read V vector\n" );
      fclose(cahvFP);
    }
    
    fclose(cahvFP);
  }

  void CAHVModel::read_pinhole(std::string const& filename) {
    FILE *camFP = fopen(filename.c_str(), "r");
    
    if (camFP == 0)
      vw_throw( IOErr() << "CAHVModel::read_pinhole: Could not open file\n" );

    char line[2048];
    double f, fH, fV, Hc, Vc;
    Vector2 pixelSize;
    Vector3 Hvec, Vvec;
  
    // Read intrinsic parameters
    fgets(line, sizeof(line), camFP);
    if (sscanf(line,"f = %lf", &f) != 1) {
      vw_throw( IOErr() << "CAHVModel::read_pinhole: Could not read focal length\n" );
      fclose(camFP);
    }

    fgets(line, sizeof(line), camFP);
    if (sscanf(line,"SP = %lf %lf", &pixelSize.x(), &pixelSize.y()) != 2) {
      vw_throw( IOErr() << "CAHVModel::read_pinhole: Could not read pixel size\n" );
      fclose(camFP);
    }

    fgets(line, sizeof(line), camFP);
    if (sscanf(line,"IC = %lf %lf", &Hc, &Vc) != 2) {
      vw_throw( IOErr() << "CAHVModel::ReadPinhole: Could not read image center pos\n" );
      fclose(camFP);
    }

    // Read extrinsic parameters
    fgets(line, sizeof(line), camFP);
    if (sscanf(line,"C = %lf %lf %lf", &C(0), &C(1), &C(2)) != 3) {
      vw_throw( IOErr() << "CAHVModel::read_pinhole: Could not read C vector\n" );
      fclose(camFP);
    }

    fgets(line, sizeof(line), camFP);
    if (sscanf(line,"A = %lf %lf %lf", &A(0), &A(1), &A(2)) != 3) {
      vw_throw( IOErr() << "CAHVModel::read_pinhole: Could not read A vector\n" );
      fclose(camFP);
    }

    fgets(line, sizeof(line), camFP);
    if (sscanf(line,"Hv = %lf %lf %lf", &Hvec(0), &Hvec(1), &Hvec(2)) != 3) {
      vw_throw( IOErr() << "CAHVModel::read_pinhole: Could not read Hvec\n" );
      fclose(camFP);
    }

    fgets(line, sizeof(line), camFP);
    if (sscanf(line,"Vv = %lf %lf %lf", &Vvec(0), &Vvec(1), &Vvec(2)) != 3) {
      vw_throw( IOErr() << "CAHVModel::read_pinhole: Could not read Vvec\n" );
      fclose(camFP);
    }

    // In the future, we should also read in a view matrix -- LJE
    //     double dummy
    //     if (sscanf(line, "VM = %lf %lf %lf %f %lf %lf %lf %f "
    // 	       "%lf %lf %lf %f %lf %lf %lf %f ",
    // 	       &Hvec(0), &Hvec(1), &Hvec(2), &dummy,
    // 	       &Vvec(0), &Vvec(1), &Vvec(2), &dummy,
    // 	       &A(0), &A(1), &A(2), &dummy,
    // 	       &C(0), &C(1), &C(2), &dummy) != 16)
    //     {
    //       vw_throw( IOErr()
    // 	<< "CAHVModel::ReadPinhole: Could not read view matrix\n" );
    //       fclose(camFP);
    //     }

    fH = f/pixelSize.x();
    fV = f/pixelSize.y();
    
    H = fH*Hvec + Hc*A;
    V = fV*Vvec + Vc*A;
    
    fclose(camFP);
  }


  // This is a re-implementation of the Epipolar math that is easier
  // to read, however it disagrees with the other epipolar code below
  // slightly.  It needs another look to clean up this discrepancy,
  // which is due to the direction chosen for Hvec. -mbroxton

//   void epipolar(CAHVModel const src_camera0, CAHVModel const src_camera1, 
//                 CAHVModel &dst_camera0, CAHVModel &dst_camera1) {

//     double fh, fv, hc, vc;
//     Vector3 f, g;
//     Vector3 A, H, V;

//     // Compute a common image center and scale for the two models
//     hc = dot_prod(src_camera0.H, src_camera0.A) / 2.0 + dot_prod(src_camera1.H, src_camera1.A) / 2.0;
//     vc = dot_prod(src_camera0.V, src_camera0.A) / 2.0 + dot_prod(src_camera1.V, src_camera1.A) / 2.0;

//     // Find the magnitude of the new H and V vectors (this will be the
//     // average of the focal lengths of these cameras).
//     f = cross_prod(src_camera0.A, src_camera0.H);
//     g = cross_prod(src_camera1.A, src_camera1.H); 
//     fh = (norm_2(f) + norm_2(g)) / 2.0;
    
//     f = cross_prod(src_camera0.A, src_camera0.V);
//     g = cross_prod(src_camera1.A, src_camera1.V);
//     fv = (norm_2(f) + norm_2(g)) / 2.0;

//     // Use common center and scale to construct an average A vector
//     Vector3 A_avg  = 0.5 * (src_camera0.A + src_camera1.A);

//     // Then adjust A to be perpindicular to the baseline between the
//     // two imagers.  This will move the epipoles to infinity.
//     Vector3 Hvec = normalize(src_camera1.C - src_camera0.C);
//     Vector3 Vvec = normalize(cross_prod(Hvec,A_avg));
//     A = normalize(cross_prod(Vvec,Hvec));       

//     std::cout << "Fh: " << fh << "   Fv: " << fv << "\n";
//     std::cout << "Ch: " << hc << "   Cv: " << vc << "\n";
//     std::cout << "Hvec: " << Hvec << "   Vvec: " << Vvec << "\n";

//     // Use the standard equations for the epipolar camera model to
//     // determine H and V given the camera intristics and the
//     // horizontal and vertical vectors.
//     H = fh*Hvec + hc*A;
//     V = fv*Vvec + vc*A;	      

//     dst_camera0.C = src_camera0.C;
//     dst_camera0.A = A;
//     dst_camera0.H = H;
//     dst_camera0.V = V;
    
//     dst_camera1.C = src_camera1.C;
//     dst_camera1.A = A;
//     dst_camera1.H = H;
//     dst_camera1.V = V;
//   }

  void epipolar(CAHVModel const src_camera0, CAHVModel const src_camera1, 
                CAHVModel &dst_camera0, CAHVModel &dst_camera1) {

    double hs, hc, vs, vc;
    Vector3 f, g, hp, ap, app, vp;
    Vector3 a, h, v;

    // Compute a common image center and scale for the two models 
    hc = dot_prod(src_camera0.H, src_camera0.A) / 2.0 + dot_prod(src_camera1.H, src_camera1.A) / 2.0;
    vc = dot_prod(src_camera0.V, src_camera0.A) / 2.0 + dot_prod(src_camera1.V, src_camera1.A) / 2.0;

    f = cross_prod(src_camera0.A, src_camera0.H);
    g = cross_prod(src_camera1.A, src_camera1.H); 
    hs = (norm_2(f) + norm_2(g)) / 2.0;
    
    f = cross_prod(src_camera0.A, src_camera0.V);
    g = cross_prod(src_camera1.A, src_camera1.V);
    vs = (norm_2(f) + norm_2(g)) / 2.0;

    // Use common center and scale to construct common A, H, V 
    app  = src_camera0.A + src_camera1.A;

    // Note the directionality of f here, for consistency later
    f = src_camera1.C - src_camera0.C;
    g = cross_prod(app, f); // alter f (CxCy) to be         
    f = cross_prod(g, app); // perpendicular to average A 
    
    if (dot_prod(f, src_camera0.H) > 0)
      hp = f * hs / (norm_2(f));
    else
      hp = -f * hs / (norm_2(f));
    
    app = 0.5 * app;
    g = hp * dot_prod(app,hp) / (hs * hs);
    ap = app -g;
    a = ap/norm_2(ap);
    f = cross_prod(a, hp);
    vp = f * vs / hs;
    f = hc * a;
    h = hp + f;
    
    f = vc * a;
    v = vp + f;
    
    dst_camera0.C = src_camera0.C;
    dst_camera0.A = a;
    dst_camera0.H = h;
    dst_camera0.V = v;
    
    dst_camera1.C = src_camera1.C;
    dst_camera1.A = a;
    dst_camera1.H = h;
    dst_camera1.V = v;
  }

}} // namespace vw::camera
