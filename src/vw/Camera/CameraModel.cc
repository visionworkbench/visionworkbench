// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#include <vw/Camera/CameraModel.h>

using namespace vw;

Vector3
vw::camera::AdjustedCameraModel::axis_angle_rotation() const {
  Quaternion<double> quat = this->rotation();
  return quat.axis_angle();
}

void
vw::camera::AdjustedCameraModel::set_translation(Vector3 const& translation) {
  m_translation = translation;
}

void
vw::camera::AdjustedCameraModel::set_rotation(Quaternion<double> const& rotation) {
  m_rotation = rotation;
  m_rotation_inverse = inverse(m_rotation);
}

void
vw::camera::AdjustedCameraModel::set_rotation(Matrix<double,3,3> const& rotation) {
  m_rotation = Quaternion<double>(rotation);
  m_rotation_inverse = inverse(m_rotation);
}

void
vw::camera::AdjustedCameraModel::set_axis_angle_rotation(Vector3 const& axis_angle) {
  Quaternion<double> rot = axis_angle_to_quaternion(axis_angle);
  this->set_rotation(rot);
}

Vector2
vw::camera::AdjustedCameraModel::point_to_pixel (Vector3 const& point) const {
  Vector3 offset_pt = point-m_camera->camera_center(Vector2(0,0))-m_translation;
  Vector3 new_pt = m_rotation_inverse.rotate(offset_pt) + m_camera->camera_center(Vector2(0,0));
  return m_camera->point_to_pixel(new_pt);
}

Vector3
vw::camera::AdjustedCameraModel::pixel_to_vector (Vector2 const& pix) const {
  return m_rotation.rotate(m_camera->pixel_to_vector(pix));
}

Vector3
vw::camera::AdjustedCameraModel::camera_center(Vector2 const& pix) const {
  return m_camera->camera_center(pix) + m_translation;
}

Quaternion<double>
vw::camera::AdjustedCameraModel::camera_pose(Vector2 const& pix) const {
  return m_rotation*m_camera->camera_pose(pix);
}

void
vw::camera::AdjustedCameraModel::write(std::string const& filename) {
  std::ofstream ostr(filename.c_str());
  ostr << m_translation[0] << " " << m_translation[1]
       << " " << m_translation[2] << "\n";
  ostr << m_rotation.w() << " " << m_rotation.x() << " "
       << m_rotation.y() << " " << m_rotation.z() << "\n";
}

void
vw::camera::AdjustedCameraModel::read(std::string const& filename) {
  Quaternion<double> rot;
  Vector3 pos;
  std::ifstream istr(filename.c_str());
  istr >> pos[0] >> pos[1] >> pos[2];
  istr >> rot.w() >> rot.x() >> rot.y() >> rot.z();
  this->set_translation(pos);
  this->set_rotation(rot);
}

std::ostream&
vw::camera::operator<<(std::ostream& ostr,
                       camera::AdjustedCameraModel const& cam ) {
  ostr << "AdjustedCameraModel(Trans: " << cam.m_translation << " Rot: "
       << cam.m_rotation << " Cam: " << cam.m_camera->type() << ")\n";
  return ostr;
}
