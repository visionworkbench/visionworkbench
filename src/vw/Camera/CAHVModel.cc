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


//
#include <vw/Core/Exception.h>
#include <vw/Core/Log.h>
#include <vw/Core/TypeDeduction.h>
#include <vw/Math/Matrix.h>
#include <vw/Math/Quaternion.h>
#include <vw/Math/Vector.h>
#include <vw/Camera/CAHVModel.h>
#include <vw/Camera/PinholeModel.h>

#include <boost/algorithm/string/predicate.hpp>

namespace vw {
namespace camera {

  /// This constructor takes a filename and reads in a camera model
  /// from the file.  The file may contain either CAHV parameters or
  /// pinhole camera parameters.
  CAHVModel::CAHVModel(std::string const& filename) {
    if (filename.empty())
      vw_throw( IOErr() << "CAHVModel: null file name passed to constructor." );

    if (boost::ends_with(filename, ".cahv"))
      read_cahv(filename);
    else if (boost::ends_with(filename, ".pin"))
      read_pinhole(filename);
    else
      vw_throw( IOErr() << "CAHVModel: Unknown camera file suffix." );
  }

  /// Initialize the CAHV vectors directly in the native CAHV format.
  CAHVModel::CAHVModel(Vector3 const& C_vec,
                       Vector3 const& A_vec,
                       Vector3 const& H_vec,
                       Vector3 const& V_vec) :
    C(C_vec), A(A_vec), H(H_vec), V(V_vec) {}

  /// Initialize a CAHV model from a pinhole model
  CAHVModel::CAHVModel(PinholeModel const& pin_model) {
    *this = pin_model;
  }

  std::string
  CAHVModel::type() const { return "CAHV"; }

  CAHVModel CAHVModel::operator=(PinholeModel const& pin_model) {

    //  Pinhole model parameters (in nominal units)
    double fH, fV, cH, cV;
    pin_model.intrinsic_parameters(fH, fV, cH, cV);
    
    // Convert the intrinsics to pixel units
    double pitch = pin_model.pixel_pitch();
    if (pitch == 0)
      vw_throw(ArgumentErr() << "CAHVModel::operator=() input pinhole model has zero pitch!");
    
    fH /= pitch;
    fV /= pitch;
    cH /= pitch;
    cV /= pitch;

    //  Unit vectors defining camera coordinate frame
    Vector3 u,v,w;
    pin_model.coordinate_frame(u,v,w);

    //  The true rotation between world and camera coordinate
    //  frames includes the rotation R --AND-- a rotation from
    //  specifying the directions of increasing u,v,w pixels
    Matrix<double,3,3> R = pin_model.get_rotation_matrix();

    //  Now create the components of the CAHV model.
    Vector3 Hvec = R*u;
    Vector3 Vvec = R*v;

    C = pin_model.camera_center();
    A = R*w;
    H = fH*Hvec + cH*A;
    V = fV*Vvec + cV*A;

    return *this;
  }

  // Convert CAHV to Pinhole. if the CAHV model was created from
  // Pinhole using operator=(PinholeModel const& pin_model), the
  // original pinhole model should be recovered. Reference:
  // https://agupubs.onlinelibrary.wiley.com/doi/epdf/10.1029/2003JE002199
  vw::camera::PinholeModel CAHVModel::toPinhole() {
    
    double Hc = dot_prod(A, H);           // optical center x
    double fH = norm_2(cross_prod(A, H)); // focal length x
    double Vc = dot_prod(A, V);           // optical center y
    double fV = norm_2(cross_prod(A, V)); // focal length y
    vw::Vector3 H_prime = (H - Hc * A)/fH;
    vw::Vector3 V_prime = (V - Vc * A)/fV;

    // Normalize
    H_prime = H_prime/norm_2(H_prime);
    V_prime = V_prime/norm_2(V_prime);
    A = A/norm_2(A);
    
    // world-to-cam
    vw::Matrix<double> R;
    R.set_size(3, 3);
    for (int col = 0; col < 3; col++) {
      R(0, col) = H_prime[col];
      R(1, col) = V_prime[col]; // plus or minus?
      // Below we do not use a minus sign, unlike the paper
      // above. That because the VW Pinhole camera model does not
      // change the sign of the x and y pixel coordinates after
      // projecting into the camera.
      R(2, col) = A[col];
    }

    // Use no distortion, and set pixel pitch to 1
    LensDistortion const* distortion = NULL;
    double pixel_pitch = 1.0;
    return vw::camera::PinholeModel(C, inverse(R), fH, fV, Hc, Vc, distortion, pixel_pitch);
  }
  
  CAHVModel::CAHVModel(double f, Vector2 const& pixel_size,
                       double xmin, double /*xmax*/, double ymin, double /*ymax*/,
                       Matrix<double, 4, 4> const& view_matrix) {

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

  CAHVModel::CAHVModel(double f, Vector2 const& pixel_size, double Hc, double Vc,
                       Vector3 const& Cinit, Vector3 const& Ainit,
                       Vector3 const& Hvec, Vector3 const& Vvec) {
    double fH = f/pixel_size.x();
    double fV = f/pixel_size.y();
    C = Cinit;
    A = Ainit;
    H = fH*Hvec + Hc*A;
    V = fV*Vvec + Vc*A;
  }

  Vector2 CAHVModel::point_to_pixel(Vector3 const& point) const {
    double dDot = dot_prod(point-C, A);
    return Vector2( dot_prod(point-C, H) / dDot,
                    dot_prod(point-C, V) / dDot );
  }

  Vector3 CAHVModel::pixel_to_vector(Vector2 const& pix) const {

    // Find vector
    Vector3 vec =
      normalize(cross_prod(V - pix.y() * A,
                           H - pix.x() * A));

    // The vector VxH should be pointing in the same directions as A,
    // if it isn't (because we have a left handed system), flip the vector.
    if (dot_prod(cross_prod(V, H), A) < 0.0)
      vec *= -1.0;
    return vec;
  }

  Vector3 CAHVModel::camera_center(Vector2 const& pix) const {
    return C;
  };

  // --------------------------------------------------
  //                 Private Methods
  // --------------------------------------------------
  void CAHVModel::read_cahv(std::string const& filename) {

    try {
      std::ifstream input(filename.c_str(), std::ifstream::in);
      input.exceptions(std::ifstream::failbit | std::ifstream::badbit);

      vw_out(InfoMessage, "camera") << "Reading CAHV file: "
                                    << filename << ".\n";

      char r1, r2;

      while (true) {
        input.ignore(1024, 'C');
        input >> r1;
        if (r1 == '=')
          break;
      }
      input >> C[0] >> C[1] >> C[2];

      input >> r1 >> r2;
      if (r1 != 'A' || r2 != '=')
        vw_throw( IOErr() << "CAHVModel: Could not read A vector\n" );
      input >> A(0) >> A(1) >> A(2);

      input >> r1 >> r2;
      if (r1 != 'H' || r2 != '=')
        vw_throw( IOErr() << "CAHVModel: Could not read H vector\n" );
      input >> H(0) >> H(1) >> H(2);

      input >> r1 >> r2;
      if (r1 != 'V' || r2 != '=')
        vw_throw( IOErr() << "CAHVModel: Could not read V vector\n" );
      input >> V(0) >> V(1) >> V(2);

    } catch ( const std::ifstream::failure& e ) {
      vw_throw( IOErr() << "CAHVModel: Could not read file: " << filename << " (" << e.what() << ")" );
    }
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
    if (!fgets(line, sizeof(line), camFP) ||
        sscanf(line,"f = %lf", &f) != 1) {
      vw_throw( IOErr() << "CAHVModel::read_pinhole: Could not read focal length\n" );
      fclose(camFP);
    }

    if (!fgets(line, sizeof(line), camFP) ||
        sscanf(line,"SP = %lf %lf", &pixelSize.x(), &pixelSize.y()) != 2) {
      vw_throw( IOErr() << "CAHVModel::read_pinhole: Could not read pixel size\n" );
      fclose(camFP);
    }

    if (!fgets(line, sizeof(line), camFP) ||
        sscanf(line,"IC = %lf %lf", &Hc, &Vc) != 2) {
      vw_throw( IOErr() << "CAHVModel::ReadPinhole: Could not read image center pos\n" );
      fclose(camFP);
    }

    // Read extrinsic parameters
    if (!fgets(line, sizeof(line), camFP) ||
        sscanf(line,"C = %lf %lf %lf", &C(0), &C(1), &C(2)) != 3) {
      vw_throw( IOErr() << "CAHVModel::read_pinhole: Could not read C vector\n" );
      fclose(camFP);
    }

    if (!fgets(line, sizeof(line), camFP) ||
        sscanf(line,"A = %lf %lf %lf", &A(0), &A(1), &A(2)) != 3) {
      vw_throw( IOErr() << "CAHVModel::read_pinhole: Could not read A vector\n" );
      fclose(camFP);
    }

    if (!fgets(line, sizeof(line), camFP) ||
        sscanf(line,"Hv = %lf %lf %lf", &Hvec(0), &Hvec(1), &Hvec(2)) != 3) {
      vw_throw( IOErr() << "CAHVModel::read_pinhole: Could not read Hvec\n" );
      fclose(camFP);
    }

    if (!fgets(line, sizeof(line), camFP) ||
        sscanf(line,"Vv = %lf %lf %lf", &Vvec(0), &Vvec(1), &Vvec(2)) != 3) {
      vw_throw( IOErr() << "CAHVModel::read_pinhole: Could not read Vvec\n" );
      fclose(camFP);
    }

    fH = f/pixelSize.x();
    fV = f/pixelSize.y();

    H = fH*Hvec + Hc*A;
    V = fV*Vvec + Vc*A;

    fclose(camFP);
  }

  void epipolar(CAHVModel const &src_camera0, CAHVModel const &src_camera1,
                CAHVModel       &dst_camera0, CAHVModel       &dst_camera1) {

    // Compute a common image center and scale for the two models
    double hc = dot_prod(src_camera0.H, src_camera0.A) / 2.0 +
                dot_prod(src_camera1.H, src_camera1.A) / 2.0;
    double vc = dot_prod(src_camera0.V, src_camera0.A) / 2.0 +
                dot_prod(src_camera1.V, src_camera1.A) / 2.0;

    double hs = norm_2(cross_prod(src_camera0.A, src_camera0.H))/2.0 +
                norm_2(cross_prod(src_camera1.A, src_camera1.H))/2.0;

    double vs = norm_2(cross_prod(src_camera0.A, src_camera0.V))/2.0 +
                norm_2(cross_prod(src_camera1.A, src_camera1.V))/2.0;

    // Use common center and scale to construct common A, H, V
    Vector3 app  = src_camera0.A + src_camera1.A;

    // Note the directionality of 1 to 0, for consistency later
    Vector3 f = cross_prod(cross_prod(app, src_camera1.C - src_camera0.C), app);

    Vector3 hp;
    if (dot_prod(f, src_camera0.H) > 0)
      hp = f * hs / (norm_2(f));
    else
      hp = -f * hs / (norm_2(f));

    app *= 0.5;
    Vector3 g  = hp * dot_prod(app,hp) / (hs * hs);
    Vector3 a  = normalize(app - g);
    Vector3 vp = cross_prod(a, hp) * vs / hs;

    // The two output cameras are identical except that they are each
    //  colocated with their corresponding input camera.
    dst_camera0.C = src_camera0.C;
    dst_camera1.C = src_camera1.C;

    dst_camera0.A = dst_camera1.A = a;
    dst_camera0.H = dst_camera1.H = hp + hc * a;
    dst_camera0.V = dst_camera1.V = vp + vc * a;
  }

  void CAHVModel::write(std::string const& filename ) {
    try {
      std::ofstream output(filename.c_str(), std::ofstream::out);
      output.exceptions(std::ofstream::failbit | std::ofstream::badbit);
      output.precision(20);

      vw_out(InfoMessage, "camera") << "Writing CAHV file: " << filename << "\n";

      output << "C = " << C[0] << " " << C[1] << " " << C[2] << "\n"
             << "A = " << A[0] << " " << A[1] << " " << A[2] << "\n"
             << "H = " << H[0] << " " << H[1] << " " << H[2] << "\n"
             << "V = " << V[0] << " " << V[1] << " " << V[2] << "\n";
    } catch ( const std::ofstream::failure& e ) {
      vw_throw( IOErr() << "CAHVModel: Could not write file: " << filename << "(" << e.what() << ")" );
    }
  }

  std::ostream& operator<<(std::ostream& str, CAHVModel const& model) {
    str << "CAHV camera: \n";
    str << "\tC: " << model.C << "\n";
    str << "\tA: " << model.A << "\n";
    str << "\tH: " << model.H << "\n";
    str << "\tV: " << model.V << "\n";
    return str;
  }

}} // namespace vw::camera
