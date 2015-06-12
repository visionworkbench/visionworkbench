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

using namespace vw;
using namespace vw::camera;

Quaternion<double>
CameraModel::camera_pose(Vector2 const& pix) const {
  vw_throw( NoImplErr() << "CameraModel: this camera model has not implemented camera_pose()" );
  return Quaternion<double>();
}

AdjustedCameraModel::AdjustedCameraModel(boost::shared_ptr<CameraModel> camera_model) : m_camera(camera_model) {
  m_rotation = Quat(math::identity_matrix<3>());
  m_rotation_inverse = Quat(math::identity_matrix<3>());
}

AdjustedCameraModel::AdjustedCameraModel(boost::shared_ptr<CameraModel> camera_model,
                                         Vector3 const& translation, Quat const& rotation, Vector2 pixel_offset) :
  m_camera(camera_model), m_translation(translation), m_rotation(rotation), m_rotation_inverse(inverse(rotation)), m_pixel_offset(pixel_offset) {}

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

Vector2 AdjustedCameraModel::point_to_pixel (Vector3 const& point) const {
  Vector3 offset_pt = point-m_camera->camera_center(Vector2(0,0))-m_translation;
  Vector3 new_pt = m_rotation_inverse.rotate(offset_pt) + m_camera->camera_center(Vector2(0,0));
  return m_camera->point_to_pixel(new_pt) - m_pixel_offset;
}

Vector3 AdjustedCameraModel::pixel_to_vector (Vector2 const& pix) const {
  return m_rotation.rotate(m_camera->pixel_to_vector(pix + m_pixel_offset));
}

Vector3 AdjustedCameraModel::camera_center(Vector2 const& pix) const {
  return m_camera->camera_center(pix + m_pixel_offset) + m_translation;
}

Quat AdjustedCameraModel::camera_pose(Vector2 const& pix) const {
  return m_rotation*m_camera->camera_pose(pix + m_pixel_offset);
}

void AdjustedCameraModel::write(std::string const& filename) {
  std::ofstream ostr(filename.c_str());
  ostr << m_translation[0] << " " << m_translation[1]
       << " " << m_translation[2] << "\n";
  ostr << m_rotation.w() << " " << m_rotation.x() << " "
       << m_rotation.y() << " " << m_rotation.z() << "\n";
  ostr << m_pixel_offset.x() << " " << m_pixel_offset.y() << "\n";

}

void AdjustedCameraModel::read(std::string const& filename) {
  Vector4 c;
  Vector3 pos;
  Vector2 pixel_offset;
  std::ifstream istr(filename.c_str());
  istr >> pos[0] >> pos[1] >> pos[2];
  istr >> c[0] >> c[1] >> c[2] >> c[3];

  if (istr >> pixel_offset[0] >> pixel_offset[1])
    m_pixel_offset = pixel_offset;

  this->set_translation(pos);
  this->set_rotation(Quat(c));
}

std::ostream& camera::operator<<(std::ostream& ostr,
                                 AdjustedCameraModel const& cam ) {
  ostr << "AdjustedCameraModel(Trans: " << cam.m_translation
       << " Rot: "          << cam.m_rotation
       << " Pixel offset: " << cam.m_pixel_offset
       << " Cam: " << cam.m_camera->type() << ")\n";
  return ostr;
}

vw::camera::CameraModel* vw::camera::unadjusted_model(vw::camera::CameraModel * cam){

  vw::camera::AdjustedCameraModel *acam
    = dynamic_cast<vw::camera::AdjustedCameraModel*>(cam);

  if (acam != NULL)
    return acam->unadjusted_model().get();

  return cam;
}

const vw::camera::CameraModel* vw::camera::unadjusted_model(const vw::camera::CameraModel * cam){

  const vw::camera::AdjustedCameraModel *acam
    = dynamic_cast<const vw::camera::AdjustedCameraModel*>(cam);

  if (acam != NULL)
    return acam->unadjusted_model().get();

  return cam;
}
