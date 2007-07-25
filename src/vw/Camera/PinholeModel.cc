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


// Reads in a file containing parameters of a pinhole model with
// a tsai lens distortion model. An example is provided at the end of this file.
void vw::camera::PinholeModel::read_file(std::string const& filename) {

  char line[2048];
  double fu, fv, cu, cv;
  Vector3 u_direction, v_direction, w_direction;
  Vector3 C;
  Matrix3x3 R;
  Vector4 distortion_params(0,0,0,0);


  FILE *cam_file = fopen(filename.c_str(), "r");
  if (cam_file == 0) vw_throw( IOErr() << "PinholeModel::read_file: Could not open file\n" );
    
  // Read intrinsic parameters
  fgets(line, sizeof(line), cam_file);
  if (sscanf(line,"fu = %lf", &fu) != 1) {
    fclose(cam_file);
    vw_throw( IOErr() << "PinholeModel::read_file(): Could not read x focal length\n" );
  }

  fgets(line, sizeof(line), cam_file);
  if (sscanf(line,"fv = %lf", &fv) != 1) {
    fclose(cam_file);
    vw_throw( IOErr() << "PinholeModel::read_file(): Could not read y focal length\n" );
  }

  fgets(line, sizeof(line), cam_file);
  if (sscanf(line,"cu = %lf", &cu) != 1) {
    fclose(cam_file);
    vw_throw( IOErr() << "PinholeModel::read_file(): Could not read x principal point\n" );
  }

  fgets(line, sizeof(line), cam_file);
  if (sscanf(line,"cv = %lf", &cv) != 1) {
    fclose(cam_file);
    vw_throw( IOErr() << "PinholeModel::read_file(): Could not read y principal point\n" );
  }

  fgets(line, sizeof(line), cam_file);
  if (sscanf(line,"u_direction = %lf %lf %lf", &u_direction(0), &u_direction(1), &u_direction(2)) != 3) {
    fclose(cam_file);
    vw_throw( IOErr() << "PinholeModel::read_file(): Could not read u direction vector\n" );
  }

  fgets(line, sizeof(line), cam_file);
  if (sscanf(line,"v_direction = %lf %lf %lf", &v_direction(0), &v_direction(1), &v_direction(2)) != 3) {
    fclose(cam_file);
    vw_throw( IOErr() << "PinholeModel::read_file(): Could not read v direction vector\n" );
  }

  fgets(line, sizeof(line), cam_file);
  if (sscanf(line,"w_direction = %lf %lf %lf", &w_direction(0), &w_direction(1), &w_direction(2)) != 3) {
    fclose(cam_file);
    vw_throw( IOErr() << "PinholeModel::read_file(): Could not read w direction vector\n" );
  }

  // Read extrinsic parameters
  fgets(line, sizeof(line), cam_file);
  if (sscanf(line,"C = %lf %lf %lf", &C(0), &C(1), &C(2)) != 3) {
    fclose(cam_file);
    vw_throw( IOErr() << "PinholeModel::read_file: Could not read C (camera center) vector\n" );
  }
  
  fgets(line, sizeof(line), cam_file);
  if ( sscanf(line, "R = %lf %lf %lf %lf %lf %lf %lf %lf %lf",
              &R(0,0), &R(0,1), &R(0,2),
              &R(1,0), &R(1,1), &R(1,2),
              &R(2,0), &R(2,1), &R(2,2)) != 9 ) {
      fclose(cam_file);
      vw_throw( IOErr() << "PinholeModel::read_file(): Could not read rotation matrix\n" );
  }

  // Read distortion parameters.
   fgets(line, sizeof(line), cam_file);
  if (sscanf(line,"k1 = %lf", &distortion_params[0] ) != 1) {
    fclose(cam_file);
    vw_throw( IOErr() << "PinholeModel::read_file(): Could not read tsai distortion parameter k1\n" );
  }
  
  fgets(line, sizeof(line), cam_file);
  if (sscanf(line,"k2 = %lf", &distortion_params[1] ) != 1) {
    fclose(cam_file);
    vw_throw( IOErr() << "PinholeModel::read_file(): Could not read tsai distortion parameter k2\n" );
  }

  fgets(line, sizeof(line), cam_file);
  if (sscanf(line,"p1 = %lf", &distortion_params[2] ) != 1) {
    fclose(cam_file);
    vw_throw( IOErr() << "PinholeModel::read_file(): Could not read tsai distortion parameter p1\n" );
  }

  fgets(line, sizeof(line), cam_file);
  if (sscanf(line,"p2 = %lf", &distortion_params[3] ) != 1) {
    fclose(cam_file);
    vw_throw( IOErr() << "PinholeModel::read_file(): Could not read tsai distortion parameter p2\n" );
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

  if( (distortion_params(0) == 0)
      && (distortion_params(1) == 0)
      && (distortion_params(2) == 0)
      && (distortion_params(3) == 0))
    {
      // Distortion model is null
      m_distortion_model_ptr = boost::shared_ptr<LensDistortion>(new NullLensDistortion());
      m_distortion_model_ptr->set_parent_camera_model(this);
    } else{
      // Create a Tsai distortion model 
      m_distortion_model_ptr = boost::shared_ptr<LensDistortion>(new TsaiLensDistortion(distortion_params));
      m_distortion_model_ptr->set_parent_camera_model(this);
    }
}


//   Write parameters of an exiting PinholeModel into a .tsai file for later use.
//   FIXME: does not output distortion parameters

void vw::camera::PinholeModel::write_file(std::string const& filename) const {
  std::ofstream cam_file(filename.c_str());
  if( !cam_file.is_open() ) vw_throw( IOErr() << "PinholeModel::write_file: Could not open file\n" );
  
  cam_file << "fu = " << m_fu << "\n";
  cam_file << "fv = " << m_fv << "\n";
  cam_file << "cu = " << m_cu << "\n";
  cam_file << "cv = " << m_cv << "\n";
  cam_file << "u_direction = " << m_u_direction[0] << " " << m_u_direction[1] << " " << m_u_direction[2] << "\n";
  cam_file << "v_direction = " << m_v_direction[0] << " " << m_v_direction[1] << " " << m_v_direction[2] << "\n";
  cam_file << "w_direction = " << m_w_direction[0] << " " << m_w_direction[1] << " " << m_w_direction[2] << "\n";
  cam_file << "C = " << m_camera_center[0] << " " << m_camera_center[1] << " " << m_camera_center[2] << "\n";
  cam_file << "R = " << m_rotation(0,0) << " " << m_rotation(0,1) << " " << m_rotation(0,2) << " " << m_rotation(1,0) << " " << m_rotation(1,1) << " " << m_rotation(1,2) << " " << m_rotation(2,0) << " " << m_rotation(2,1) << " " << m_rotation(2,2) << "\n";

  //  Distortion Parameters. This should be implemented by overloading the
  //  << operator for distortion models

  cam_file << m_distortion_model_ptr << "\n";
  cam_file << "\n" << "\n" << " Parameters for a Pinhole camera model with tsai lens distortion model." << "\n";
  cam_file.close();
}



/* Contents of a sample .tsai file:

fu = 611.651
fv = 610.216
cu = 500.829
cv = 396.22
u_direction = 1 0 0
v_direction = 0 1 0
w_direction = 0 0 1
C = -0.328711 -0.0637059 -0.828905
R = 0.000412095 -0.99998 0.00624732 0.409245 0.00586886 0.912405 -0.912424 0.00218069 0.40924
k1 = 0
k2 = 0
p1 = 0
p2 = 0


 Parameters for a Pinhole camera model with tsai lens distortion model.


*/
