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


#include <vw/Math/Functors.h>
#include <vw/Math/Quaternion.h>
#include <vw/Camera/CameraModel.h>

#include <sstream>
#include <iostream>


std::ostream& vw::camera::operator<<(std::ostream& ostr,
                                     vw::camera::AdjustedCameraModel const& cam ) {
  ostr.precision(18);
  ostr << "AdjustedCameraModel(Trans: " << cam.m_translation
       << " Rot:          " << cam.m_rotation
       << " Pixel offset: " << cam.m_pixel_offset
       << " Scale:        " << cam.m_scale
       << " Cam:          " << cam.m_camera->type() << ")\n";
  return ostr;
}

namespace vw {
namespace camera {

Quaternion<double>
CameraModel::camera_pose(Vector2 const& pix) const {
  vw_throw( NoImplErr() << "CameraModel: this camera model has not implemented camera_pose()" );
  return Quaternion<double>();
}

AdjustedCameraModel::AdjustedCameraModel(boost::shared_ptr<CameraModel> camera_model,
                                         Vector3 const& translation, Quat const& rotation,
                                         Vector2 const& pixel_offset, double scale) :
  m_camera(camera_model),
  m_translation(translation),
  m_rotation(rotation),
  m_rotation_inverse(inverse(rotation)),
  m_pixel_offset(pixel_offset),
  m_scale(scale) {

  // Set as the rotation center the old camera center for pixel (0,0).
  // It is important to note that for linescan cameras, each line has
  // its own camera center. It is not a problem that we use a single
  // fixed rotation center. We just need to pick some point to rotate
  // about. Ideally it would be the camera center for the pixel in the
  // middle of the image, if we know the image size. But the current
  // choice should be just as good.

  // This can throw an exception, so be careful. If (0, 0) does not
  // succeed, keep on trying on the diagonal.
  m_rotation_center = Vector3(); // Initialize to center of planet (for lack of anything better)
  for (int i = 0; i < 100000; i+= 10) {
    try {
      m_rotation_center = m_camera->camera_center(Vector2(i, i));
      break; // success
    }catch(...){}
  }
}

AdjustedCameraModel::~AdjustedCameraModel() {}

std::string AdjustedCameraModel::type() const { return "Adjusted"; }

Vector3 AdjustedCameraModel::translation() const { return m_translation; }

Quat AdjustedCameraModel::rotation() const { return m_rotation; }

Matrix<double,3,3> AdjustedCameraModel::rotation_matrix() const { return m_rotation.rotation_matrix(); }

Vector3 AdjustedCameraModel::axis_angle_rotation() const {
  Quat quat = this->rotation();
  return quat.axis_angle();
}

void AdjustedCameraModel::set_rotation(Quat const& rotation) {
  m_rotation = rotation;
  m_rotation_inverse = inverse(m_rotation);
}

Vector2 AdjustedCameraModel::pixel_offset() const { return m_pixel_offset; }

double AdjustedCameraModel::scale() const { return m_scale; }

Vector3 AdjustedCameraModel::adjusted_point(Vector3 const& point) const {
  Vector3 offset_pt = point-m_rotation_center-m_translation;
  Vector3 new_pt = m_rotation_inverse.rotate(offset_pt) + m_rotation_center;
  return new_pt;
}

Vector2 AdjustedCameraModel::point_to_pixel (Vector3 const& point) const {
  Vector3 new_pt = this->adjusted_point(point);
  return (m_camera->point_to_pixel(new_pt) - m_pixel_offset)/m_scale;
}

Vector3 AdjustedCameraModel::pixel_to_vector (Vector2 const& pix) const {
  return m_rotation.rotate(m_camera->pixel_to_vector(m_scale*pix + m_pixel_offset));
}

Vector3 AdjustedCameraModel::camera_center(Vector2 const& pix) const {
  Vector3 old_ct = m_camera->camera_center(m_scale*pix + m_pixel_offset);
  return m_rotation.rotate(old_ct - m_rotation_center) + m_rotation_center + m_translation;
  // Old and incorrect formula below
  //return m_camera->camera_center(m_scale*pix + m_pixel_offset) + m_translation;
}

Quat AdjustedCameraModel::camera_pose(Vector2 const& pix) const {
  return m_rotation*m_camera->camera_pose(m_scale*pix + m_pixel_offset);
}

// Modify the adjustments by applying on top of them a scale*rotation + translation
// transform with the origin at the center of the planet (such as output
// by pc_align's forward or inverse computed alignment transform). 
void AdjustedCameraModel::apply_transform(vw::Matrix4x4 const& M){

  Matrix3x3 R = submatrix(M, 0, 0, 3, 3);
  Vector3   T;
  for (int r = 0; r < 3; r++) 
    T[r] = M(r, 3);

  double scale = pow(det(R), 1.0/3.0);
  for (size_t r = 0; r < R.rows(); r++)
    for (size_t c = 0; c < R.cols(); c++)
      R(r, c) /= scale;
  
  // The logic is the following. A point on the camera body gets transformed
  // by the given adjustments as:
  // x -> m_rotation.rotate(x - m_rotation_center) + m_rotation_center + m_translation.
  // If on top of that we apply the current R and T,
  // that would be:
  // R*(m_rotation.rotate(x - m_rotation_center) + m_rotation_center + m_translation) + T
  // which is same as
  // (R*m_rotation).rotate(x-m_rotation_center) + m_rotation_center +
  // R*(m_rotation_center + m_translation) + T - m_rotation_center.
  // Note that below we also incorporated the scale.
  m_rotation    = Quat(R)*m_rotation;
  m_translation = R*scale*(m_rotation_center + m_translation) + T - m_rotation_center;
}
  
// Return the camera adjustment as a 4x4 rotation + translation
// transform (upper-left 3x3 block is the rotation, the first 3
// elements in last column is the translation, so the style of
// pc_align is used) as when applied around the planet center.
// The adjusted camera center and pixel to vector are obtained
// from the unadjusted ones by applying the translation
// and rotation returned here.
// No scale factor is applied here.
vw::Matrix4x4 AdjustedCameraModel::ecef_transform() const {

  // Use the fact that the adjustment transform applied to a point x
  // on the camera acts like m_rotation.rotate(x - m_rotation_center)
  // + m_rotation_center + m_translation.

  // Throw error if there is a scale factor.
  if (std::abs(m_scale - 1.0) > 1e-10)
    vw::vw_throw(vw::NoImplErr() 
      << "ecef_transform() does not support a scale factor.\n");
    
  Matrix3x3 R = rotation_matrix();
  Vector3   S = -R * m_rotation_center + m_rotation_center + m_translation;

  Matrix4x4 T;
  T.set_identity();
    
  submatrix(T, 0, 0, 3, 3) = R;
  for (int row = 0; row < 3; row++) 
    T(row, 3) = S[row];

  return T;
}
  
void AdjustedCameraModel::write(std::string const& filename) {
  std::ofstream ostr(filename.c_str());
  ostr.precision(18);
  ostr << m_translation[0] << " " << m_translation[1]
       << " " << m_translation[2] << "\n";
  ostr << m_rotation.w() << " " << m_rotation.x() << " "
       << m_rotation.y() << " " << m_rotation.z() << "\n";
  ostr << m_pixel_offset.x() << " " << m_pixel_offset.y() << "\n";
  ostr << m_scale << "\n";

}

// TODO(oalexan1): Support reading and writing m_rotation_center,
// rather than always computing it on the fly. 
  
void AdjustedCameraModel::read(std::string const& filename) {
  Vector4 c;
  Vector3 pos;
  Vector2 pixel_offset;
  double scale = 1.0;

  // Start with a blank model, and try to read data into it
  *this = AdjustedCameraModel(m_camera);

  std::ifstream istr(filename.c_str());
  if (istr >> pos[0] >> pos[1] >> pos[2])
    this->set_translation(pos);

  if (istr >> c[0] >> c[1] >> c[2] >> c[3])
    this->set_rotation(Quat(c));

  if (istr >> pixel_offset[0] >> pixel_offset[1])
    this->set_pixel_offset(pixel_offset);

  if (istr >> scale)
    this->set_scale(scale);
}

CameraModel* unadjusted_model(CameraModel * cam){

  AdjustedCameraModel *acam
    = dynamic_cast<AdjustedCameraModel*>(cam);

  if (acam != NULL)
    return acam->unadjusted_model().get();

  return cam;
}

const CameraModel* unadjusted_model(const CameraModel * cam){

  const AdjustedCameraModel *acam
    = dynamic_cast<const AdjustedCameraModel*>(cam);

  if (acam != NULL)
    return acam->unadjusted_model().get();

  return cam;
}

boost::shared_ptr<CameraModel> unadjusted_model(boost::shared_ptr<CameraModel> cam) {

  AdjustedCameraModel *acam
    = dynamic_cast<AdjustedCameraModel*>(cam.get());

  if (acam != NULL)
    return acam->unadjusted_model();

  return cam;
}

} // End namespace camera
} // End namespace vw
