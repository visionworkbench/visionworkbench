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

// Vision Workbench
#include <vw/Image/ImageView.h>

// Proj.4
#include <projects.h>


namespace vw {
namespace cartography {

  // Constructor
    GeoTransform::GeoTransform(GeoReference const& src_georef, GeoReference const& dst_georef) :
    m_src_georef(src_georef), m_dst_georef(dst_georef) {
    const std::string src_datum = m_src_georef.datum().proj4_str();
    const std::string dst_datum = m_dst_georef.datum().proj4_str();

    // This optimizes in the common case where the two images are
    // already in the same map projection, and we need only apply
    // the affine transform.  This will break, of course, as soon as
    // we have mare than one type of GeoReference object, but it
    // makes life faster for now. -mbroxton
    if (m_src_georef.proj4_str() == m_dst_georef.proj4_str())
      m_skip_map_projection = true;
    else
      m_skip_map_projection = false;

    // This optimizes the case where the two datums are the same, 
    // and thus we don't need to call proj to convert between them 
    // as we transform.
    if(src_datum == dst_datum) {
      m_skip_datum_conversion = true;
    } else {
      // Set up the various variables for proj.
      m_skip_datum_conversion = false;

      // The source proj4 context.
      std::stringstream ss_src;
      // We convert lat/long to lat/long regardless of what the 
      // source or destination georef uses.
      ss_src << "+proj=latlong " << src_datum;
      m_src_datum = boost::shared_ptr<ProjContext>(new ProjContext(ss_src.str()));
      CHECK_PROJ_INIT_ERROR(ss_src.str().c_str());

      // The destination proj4 context.
      std::stringstream ss_dst;
      ss_dst << "+proj=latlong " << dst_datum;
      m_dst_datum = boost::shared_ptr<ProjContext>(new ProjContext(ss_dst.str()));
      CHECK_PROJ_INIT_ERROR(ss_dst.str().c_str());
    }
  }
  // Performs a forward or reverse datum conversion.
  Vector2 GeoTransform::datum_convert(Vector2 const& v, bool forward) const {
    double x = v[0];
    double y = v[1];
    double z = 0;

    if(forward)
      pj_transform(m_src_datum->proj_ptr(), m_dst_datum->proj_ptr(), 1, 0, &x, &y, &z);
    else
      pj_transform(m_dst_datum->proj_ptr(), m_src_datum->proj_ptr(), 1, 0, &x, &y, &z);
    CHECK_PROJ_ERROR;

    return Vector2(x, y);
  }

  void reproject_point_image(ImageView<Vector3> const& point_image,
                             GeoReference const& src_georef,
                             GeoReference const& dst_georef) {

    GeoTransform gtx(src_georef, dst_georef);

    // Iterate over the image, transforming the first two coordinates
    // in the Vector one at a time.  The third coordinate is taken to
    // be the altitude value, and this value is not touched.
    for (int32 j=0; j < point_image.rows(); ++j) {
      for (int32 i=0; i < point_image.cols(); ++i) {
        if (point_image(i,j) != Vector3()) {
          Vector2 in(point_image(i,j)[0], point_image(i,j)[1]);
          Vector2 out = gtx.forward(in);
          point_image(i,j).x() = out[0];
          point_image(i,j).y() = out[1];
        }
      }
    }
  }

  
}} // namespace vw::cartography

