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
#ifndef __VW_CARTOGRAPHY_GEOTRANSFORM_H__
#define __VW_CARTOGRAPHY_GEOTRANSFORM_H__

#include <vw/Math/Vector.h>
#include <vw/Image/Transform.h>
#include <vw/Cartography/GeoReference.h>

namespace vw {
namespace cartography {

  class GeoTransform : public TransformBase<GeoTransform> {
    
    GeoReference m_src_georef, m_dst_georef;
    Matrix<double,3,3> m_inv_src_transform, m_inv_dst_transform;; 
    std::string src_proj_str, dst_proj_str;
    void *src_proj, *dst_proj;
    
    void init_from_georefs(GeoReference const& src_georef, GeoReference const& dst_georef);

  public:
    /// Normal constructor
    GeoTransform(GeoReference const& src_georef, GeoReference const& dst_georef) :
      m_src_georef(src_georef), m_dst_georef(dst_georef) {
      init_from_georefs(src_georef, dst_georef);
      // For debugging:
      //       std::cout << src_proj_str << "\n";
      //       std::cout << dst_proj_str << "\n";
      VW_ASSERT(src_proj_str.size() != 0, ArgumentErr() << "GeoTransform: source georeference not sufficiently well-defined.  Proj.4 string was empty.");
      VW_ASSERT(dst_proj_str.size() != 0, ArgumentErr() << "GeoTransform: destination georeference not sufficiently well-defined.  Proj.4 string was empty.");
    }

    /// Copy Constructor
    GeoTransform(GeoTransform const& copy) {
      m_src_georef = copy.m_src_georef;
      m_dst_georef = copy.m_dst_georef;
      init_from_georefs(copy.m_src_georef, copy.m_dst_georef);
    }

    /// Copy Assignment
    GeoTransform& operator=(GeoTransform const& copy) {
      m_src_georef = copy.m_src_georef;
      m_dst_georef = copy.m_dst_georef;
      init_from_georefs(copy.m_src_georef, copy.m_dst_georef);
      return *this;
    }

    /// Destructor
    ~GeoTransform();
    
    /// Given a pixel coordinate of an image in a destination
    /// georeference frame, this routine computes the corresponding
    /// pixel from an image in the source georeference frame.
    Vector2 reverse(Vector2 const& v) const;

    /// Given a pixel coordinate of an image in a source
    /// georeference frame, this routine computes the corresponding
    /// pixel the destination (transformed) image.
    Vector2 forward(Vector2 const& v) const;
  };


}} // namespace vw::cartography

#endif // __GEO_TRANSFORM_H__
