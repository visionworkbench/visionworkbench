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
#include <vw/Cartography/GeoTransform.h>

// libProj.4
#include <projects.h>

// Boost
#include <boost/algorithm/string.hpp>

// Vision Workbench
#include <vw/Image/ImageView.h>

static char** split_proj4_string(std::string const& proj4_str, int &num_strings) {
  std::vector<std::string> arg_strings;
  std::string trimmed_proj4_str = boost::trim_copy(proj4_str);
  boost::split( arg_strings, trimmed_proj4_str, boost::is_any_of(" ") ); 
      
  char** strings = new char*[arg_strings.size()]; 
  for ( unsigned i = 0; i < arg_strings.size(); ++i ) {
    strings[i] = new char[2048];
    strncpy(strings[i], arg_strings[i].c_str(), 2048);
  }
  num_strings = arg_strings.size();
  return strings;
}

namespace vw {
namespace cartography {

  void GeoTransform::init_from_georefs(GeoReference const& src_georef, GeoReference const& dst_georef) {
    // Set up two proj.4 projections
    src_proj_str = src_georef.proj4_str();
    dst_proj_str = dst_georef.proj4_str();
    
    // proj.4 is expecting the parameters to be split up into seperate
    // c-style strings.  
    int src_n, dst_n;
    char** src_strings = split_proj4_string(src_proj_str, src_n);
    char** dst_strings = split_proj4_string(dst_proj_str, dst_n);
    
    src_proj = pj_init(src_n, src_strings);
    dst_proj = pj_init(dst_n, dst_strings);
    
    for (int i = 0; i < src_n; i++) 
      delete [] src_strings[i];
    for (int i = 0; i < dst_n; i++) 
      delete [] dst_strings[i];
    delete [] src_strings;
    delete [] dst_strings;
      
    m_inv_src_transform = inverse(src_georef.transform());
    m_inv_dst_transform = inverse(dst_georef.transform());
  }

  GeoTransform::~GeoTransform() {
    pj_free((PJ*)src_proj);
    pj_free((PJ*)dst_proj);

  }
  
  Vector2 GeoTransform::reverse(Vector2 const& v) const {
    XY projected;  
    LP unprojected;
    
    // Apply the affine transformation to the projected point from
    // the destination projection.  This transforms a coordinate in
    // pixel space to a coordinate in the geographically projected
    // space.
    Matrix<double,3,3> dst_transform = m_dst_georef.transform();
    projected.u = v[0] * dst_transform(0,0) + v[1] * dst_transform(0,1) + dst_transform(0,2) / 
      (dst_transform(2,0) + dst_transform(2,1) + dst_transform(2,2));
    projected.v = v[0] * dst_transform(1,0) + v[1] * dst_transform(1,1) + dst_transform(1,2) /
      (dst_transform(2,0) + dst_transform(2,1) + dst_transform(2,2));

    // Proj.4 expects the (lon,lat) pair to be in radians, so we
    // must make a conversion if the CS in geographic (lat/lon).
    if ( !m_dst_georef.is_projected() ) {
      projected.u *= DEG_TO_RAD;
      projected.v *= DEG_TO_RAD;
    }
    
    // Apply the chain of transformations to move the point from the
    // destination projection to the source projection.
    unprojected = pj_inv(projected, (PJ*)dst_proj);
    projected = pj_fwd(unprojected, (PJ*)src_proj);
    
    // Proj.4 expects the (lon,lat) pair to be in radians, so we
    // must make a conversion if the CS in geographic (lat/lon).
    if ( !m_src_georef.is_projected() ) {
      projected.u *= RAD_TO_DEG;
      projected.v *= RAD_TO_DEG;
    }
    
    // Return the point multiplied by the inverse of the source
    // affine transform.  This places the projected coordinate back
    // into pixel space.
    return Vector2(projected.u * m_inv_src_transform(0,0) + projected.v * m_inv_src_transform(0,1) + m_inv_src_transform(0,2) / 
                     (m_inv_src_transform(2,0) + m_inv_src_transform(2,1) + m_inv_src_transform(2,2)), 
                   projected.u * m_inv_src_transform(1,0) + projected.v * m_inv_src_transform(1,1) + m_inv_src_transform(1,2) / 
                     (m_inv_src_transform(2,0) + m_inv_src_transform(2,1) + m_inv_src_transform(2,2)));
  }

  Vector2 GeoTransform::forward(Vector2 const& v) const {
    XY projected;  
    LP unprojected;
    Matrix<double,3,3> src_transform = m_src_georef.transform();
    projected.u = v[0] * src_transform(0,0) + v[1] * src_transform(0,1) + src_transform(0,2) / 
      (src_transform(2,0) + src_transform(2,1) + src_transform(2,2));
    projected.v = v[0] * src_transform(1,0) + v[1] * src_transform(1,1) + src_transform(1,2) /
      (src_transform(2,0) + src_transform(2,1) + src_transform(2,2));

    if ( !m_src_georef.is_projected() ) {
      projected.u *= DEG_TO_RAD;
      projected.v *= DEG_TO_RAD;
    }
    
    unprojected = pj_inv(projected, (PJ*)src_proj);
    projected = pj_fwd(unprojected, (PJ*)dst_proj);
    
    if ( !m_dst_georef.is_projected() ) {
      projected.u *= RAD_TO_DEG;
      projected.v *= RAD_TO_DEG;
    }

    return Vector2(projected.u * m_inv_dst_transform(0,0) + projected.v * m_inv_dst_transform(0,1) + m_inv_dst_transform(0,2) / (m_inv_dst_transform(2,0) + m_inv_dst_transform(2,1) + m_inv_dst_transform(2,2)), 
                   projected.u * m_inv_dst_transform(1,0) + projected.v * m_inv_dst_transform(1,1) + m_inv_dst_transform(1,2) / (m_inv_dst_transform(2,0) + m_inv_dst_transform(2,1) + m_inv_dst_transform(2,2)));
  }
    

  void reproject_point_image(ImageView<Vector3> const& point_image,
                             GeoReference const& src_georef,
                             GeoReference const& dst_georef) {

    PJ *src_proj, *dst_proj;

    std::string src_proj_str = src_georef.proj4_str();
    std::string dst_proj_str = dst_georef.proj4_str();

    // proj.4 is expecting the parameters to be split up into seperate
    // c-style strings.  
    int src_n, dst_n;
    char** src_strings = split_proj4_string(src_proj_str, src_n);
    char** dst_strings = split_proj4_string(dst_proj_str, dst_n);
    
    src_proj = pj_init(src_n, src_strings);
    dst_proj = pj_init(dst_n, dst_strings);

    // Delete the split up strings allocated in split_proj4_string()
    for (int i = 0; i < src_n; i++) 
      delete [] src_strings[i];
    for (int i = 0; i < dst_n; i++) 
      delete [] dst_strings[i];
    delete [] src_strings;
    delete [] dst_strings;

    XY projected;  
    LP unprojected;

    // Iterate over the image, transforming the first two coordinates
    // in the Vector one at a time.  The third coordinate is taken to
    // be the altitude value, and this value is not touched.
    for (unsigned j=0; j < point_image.rows(); ++j) {
      for (unsigned i=0; i < point_image.cols(); ++i) {
        if (point_image(i,j) != Vector3()) {
          
          if ( !src_georef.is_projected() ) {
            projected.u = point_image(i,j).x() * DEG_TO_RAD;
            projected.v = point_image(i,j).y() * DEG_TO_RAD;
          } else {
            projected.u = point_image(i,j).x();
            projected.v = point_image(i,j).y();
          }

          unprojected = pj_inv(projected, src_proj);
          projected = pj_fwd(unprojected, dst_proj);
          
          point_image(i,j).x() = projected.u;
          point_image(i,j).y() = projected.v;
        }
      }
    }

    pj_free(src_proj);
    pj_free(dst_proj);
  }

  
}} // namespace vw::cartography

