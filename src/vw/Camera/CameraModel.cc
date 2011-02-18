// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#include <vw/Camera/CameraModel.h>

using namespace vw;
using namespace vw::camera;

Vector3 AdjustedCameraModel::axis_angle_rotation() const {
  Quat quat = this->rotation();
  return quat.axis_angle();
}

void AdjustedCameraModel::set_rotation(Quat const& rotation) {
  m_rotation = rotation;
  m_rotation_inverse = inverse(m_rotation);
}

Vector2 AdjustedCameraModel::point_to_pixel (Vector3 const& point) const {
  Vector3 offset_pt = point-m_camera->camera_center(Vector2(0,0))-m_translation;
  Vector3 new_pt = m_rotation_inverse.rotate(offset_pt) + m_camera->camera_center(Vector2(0,0));
  return m_camera->point_to_pixel(new_pt);
}

Vector3 AdjustedCameraModel::pixel_to_vector (Vector2 const& pix) const {
  return m_rotation.rotate(m_camera->pixel_to_vector(pix));
}

Vector3 AdjustedCameraModel::camera_center(Vector2 const& pix) const {
  return m_camera->camera_center(pix) + m_translation;
}

Quat AdjustedCameraModel::camera_pose(Vector2 const& pix) const {
  return m_rotation*m_camera->camera_pose(pix);
}

void AdjustedCameraModel::write(std::string const& filename) {
  std::ofstream ostr(filename.c_str());
  ostr << m_translation[0] << " " << m_translation[1]
       << " " << m_translation[2] << "\n";
  ostr << m_rotation.w() << " " << m_rotation.x() << " "
       << m_rotation.y() << " " << m_rotation.z() << "\n";
}

void AdjustedCameraModel::read(std::string const& filename) {
  Vector4 c;
  Vector3 pos;
  std::ifstream istr(filename.c_str());
  istr >> pos[0] >> pos[1] >> pos[2];
  istr >> c[0] >> c[1] >> c[2] >> c[3];
  this->set_translation(pos);
  this->set_rotation(Quat(c));
}

std::ostream& camera::operator<<(std::ostream& ostr,
           AdjustedCameraModel const& cam ) {
  ostr << "AdjustedCameraModel(Trans: " << cam.m_translation << " Rot: "
       << cam.m_rotation << " Cam: " << cam.m_camera->type() << ")\n";
  return ostr;
}
