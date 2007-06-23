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


#include <vw/Camera/PinholeModel.h>
#include <vw/Math/EulerAngles.h>

void vw::camera::PinholeModel::read_file(std::string const& filename) {

  char line[2048];
  double fu, fv, cu, cv;
  Vector3 u_direction, v_direction, w_direction;
  Vector3 C;
  Matrix3x3 R;

  FILE *cam_file = fopen(filename.c_str(), "r");
  if (cam_file == 0) vw_throw( IOErr() << "CAHVModel::read_pinhole: Could not open file\n" );
    
  // Read intrinsic parameters
  fgets(line, sizeof(line), cam_file);
  if (sscanf(line,"fu = %lf", &fu) != 1) {
    fclose(cam_file);
    vw_throw( IOErr() << "read_pinhole(): Could not read x focal length\n" );
  }

  fgets(line, sizeof(line), cam_file);
  if (sscanf(line,"fv = %lf", &fv) != 1) {
    fclose(cam_file);
    vw_throw( IOErr() << "read_pinhole(): Could not read y focal length\n" );
  }

  fgets(line, sizeof(line), cam_file);
  if (sscanf(line,"cu = %lf", &cu) != 1) {
    fclose(cam_file);
    vw_throw( IOErr() << "read_pinhole(): Could not read x principal point\n" );
  }

  fgets(line, sizeof(line), cam_file);
  if (sscanf(line,"cv = %lf", &cv) != 1) {
    fclose(cam_file);
    vw_throw( IOErr() << "read_pinhole(): Could not read y principal point\n" );
  }

  fgets(line, sizeof(line), cam_file);
  if (sscanf(line,"u_direction = %lf %lf %lf", &u_direction(0), &u_direction(1), &u_direction(2)) != 3) {
    fclose(cam_file);
    vw_throw( IOErr() << "read_pinhole(): Could not read u direction vector\n" );
  }

  fgets(line, sizeof(line), cam_file);
  if (sscanf(line,"v_direction = %lf %lf %lf", &v_direction(0), &v_direction(1), &v_direction(2)) != 3) {
    fclose(cam_file);
    vw_throw( IOErr() << "read_pinhole(): Could not read v direction vector\n" );
  }

  fgets(line, sizeof(line), cam_file);
  if (sscanf(line,"w_direction = %lf %lf %lf", &w_direction(0), &w_direction(1), &w_direction(2)) != 3) {
    fclose(cam_file);
    vw_throw( IOErr() << "read_pinhole(): Could not read w direction vector\n" );
  }

  // Read extrinsic parameters
  fgets(line, sizeof(line), cam_file);
  if (sscanf(line,"C = %lf %lf %lf", &C(0), &C(1), &C(2)) != 3) {
    fclose(cam_file);
    vw_throw( IOErr() << "read_pinhole: Could not read C (camera center) vector\n" );
  }
  
  fgets(line, sizeof(line), cam_file);
  if ( sscanf(line, "R = %lf %lf %lf %lf %lf %lf %lf %lf %lf",
              &R(0,0), &R(0,1), &R(0,2),
              &R(1,0), &R(1,1), &R(1,2),
              &R(2,0), &R(2,1), &R(2,2)) != 9 ) {
      fclose(cam_file);
      vw_throw( IOErr() << "read_pinhole(): Could not read rotation matrix\n" );
  }
  fclose(cam_file);
  
  m_u_direction = u_direction;
  m_v_direction = v_direction;
  m_w_direction = w_direction;

  m_fu = fu;
  m_fv = fv;
  m_cu = cu;
  m_cv = cv;
  m_camera_center = C;
  
  m_rotation = R;
  this->rebuild_camera_matrix();
}
