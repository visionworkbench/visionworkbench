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

// This file is included in LinescanModel.h

template <class PositionFuncT, class PoseFuncT>
LinescanModel<PositionFuncT, PoseFuncT>::LinescanModel( int     number_of_lines,
                                                        int     samples_per_line,
                                                        int     sample_offset,
                                                        double  focal_length,
                                                        double  along_scan_pixel_size,
                                                        double  across_scan_pixel_size,
                                                        std::vector<double> const& line_times,
                                                        Vector3 pointing_vec,
                                                        Vector3 u_vec,
                                                        PositionFuncT const& position_func,
                                                        PoseFuncT const& pose_func) : m_position_func(position_func),
                                                                                      m_pose_func(pose_func) {

  VW_ASSERT(int(line_times.size()) == number_of_lines,
            ArgumentErr() << "LinescanModel: number of line integration times does not match the number of scanlines.\n");

  // Intrinsics
  m_number_of_lines  = number_of_lines;
  m_samples_per_line = samples_per_line;
  m_sample_offset    = sample_offset;
  m_focal_length     = focal_length;
  m_along_scan_pixel_size  = along_scan_pixel_size;
  m_across_scan_pixel_size = across_scan_pixel_size;

  m_line_times = line_times;

  m_pointing_vec = normalize(pointing_vec);
  m_u_vec        = normalize(u_vec);
}

template <class PositionFuncT, class PoseFuncT>
LinescanModel<PositionFuncT, PoseFuncT>::LinescanModel( int     number_of_lines,
                                                        int     samples_per_line,
                                                        int     sample_offset,
                                                        double  focal_length,
                                                        double  along_scan_pixel_size,
                                                        double  across_scan_pixel_size,
                                                        double  line_integration_time,
                                                        Vector3 pointing_vec,
                                                        Vector3 u_vec,
                                                        PositionFuncT const& position_func,
                                                        PoseFuncT const& pose_func) : m_position_func(position_func),
                                                                                      m_pose_func(pose_func) {

  // Intrinsics
  m_number_of_lines  = number_of_lines;
  m_samples_per_line = samples_per_line;
  m_sample_offset    = sample_offset;
  m_focal_length     = focal_length;
  m_along_scan_pixel_size = along_scan_pixel_size;
  m_across_scan_pixel_size = across_scan_pixel_size;

  m_line_times.resize(number_of_lines);
  double sum = 0;
  for (int i = 0; i < number_of_lines; ++i) {
    m_line_times[i] = sum;
    sum += line_integration_time;
  }

  m_pointing_vec = normalize(pointing_vec);
  m_u_vec        = normalize(u_vec);
}

template <class PositionFuncT, class PoseFuncT>
Vector3 LinescanModel<PositionFuncT, PoseFuncT>::pixel_to_vector(Vector2 const& pix) const {
  double u = pix[0]; // Column
  double v = pix[1]; // Line

  // Check to make sure that this is a valid pixel
  if (int(round(v)) < 0 || int(round(v)) >= int(m_line_times.size()))
    vw_throw( PixelToRayErr() << "LinescanModel: requested pixel "
              << pix << " is not on a valid scanline." );

  // The view_matrix takes vectors from the camera (extrinsic)
  // coordinate system to the world frame
  //
  // The position and veloctiy are not actually needed, since we are
  // purely interested in returning the direction of the ray at this
  // point and not its origin.
  //
  double approx_line_time = interp_line_time(pix[1]);

  // The viewplane is the [pointing_vec cross u_vec] plane of the
  // camera coordinate system.  Assuming the origin of the
  // coordinate system is at the center of projection, the image
  // plane is at pointing_vec = +f, and the pixel position in
  // camera coordinates is:
  double pixel_pos_u = (u + m_sample_offset) * m_across_scan_pixel_size;
  Vector3 pixel_direction =
    pixel_pos_u * m_u_vec + m_focal_length * m_pointing_vec;

  // Transform to world coordinates using the rigid rotation
  Quat pose = m_pose_func(approx_line_time);
  return normalize(pose.rotate(pixel_direction));
}

template <class PositionFuncT, class PoseFuncT>
Vector3 LinescanModel<PositionFuncT, PoseFuncT>::camera_center(Vector2 const& pix ) const {
  // Check to make sure that this is a valid pixel
  if (int(round(pix[1])) < 0 ||
      int(round(pix[1])) >= int(m_line_times.size()))
    vw_throw( PixelToRayErr() << "LinescanModel: requested pixel " << pix << " is not on a valid scanline." );

  double approx_line_time = interp_line_time(pix[1]);

  return m_position_func(approx_line_time);
}

template <class PositionFuncT, class PoseFuncT>
Quaternion<double> LinescanModel<PositionFuncT, PoseFuncT>::camera_pose(Vector2 const& pix) const {
  // Check to make sure that this is a valid pixel
  if (int(round(pix[1])) < 0 || int(round(pix[1])) >= int(m_line_times.size()))
    vw_throw( PixelToRayErr() << "LinescanModel::camera_pose(): requested pixel " << pix << " is not on a valid scanline." );

  double approx_line_time = interp_line_time(pix[1]);

  return m_pose_func(approx_line_time);
}

template <class PositionFuncT, class PoseFuncT>
double LinescanModel<PositionFuncT, PoseFuncT>::interp_line_time(double row) const {
  int    y     = int(floor(row));
  double normy = row - y;
  double approx_line_time = double( m_line_times[y] + (m_line_times[y+1] - m_line_times[y]) * normy );
  return approx_line_time;
}


template <class PositionFuncT, class PoseFuncT>
std::ostream& operator<<( std::ostream& os, LinescanModel<PositionFuncT, PoseFuncT> const& camera_model) {
  os << "\n-------------------- Linescan Camera Model -------------------\n\n";
  os << " Camera center @ origin :   " << camera_model.camera_center(Vector2(0,0)) << "\n";
  os << " Number of Lines        :   " << camera_model.number_of_lines()        << "\n";
  os << " Samples per Line       :   " << camera_model.samples_per_line()       << "\n";
  os << " Sample Offset          :   " << camera_model.sample_offset()          << "\n";
  os << " Focal Length           :   " << camera_model.focal_length()           << "\n";
  os << " Across Scan Pixel Size :   " << camera_model.across_scan_pixel_size() << "\n";
  os << " Along Scan Pixel Size  :   " << camera_model.along_scan_pixel_size()  << "\n";
  os << "\n------------------------------------------------------------------------\n\n";
  return os;
}

