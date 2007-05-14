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

  class GeoTransform : public TransformHelper<GeoTransform,ContinuousFunction,ContinuousFunction> {
    
    GeoReferenceBase const& m_src_georef;
    GeoReferenceBase const& m_dst_georef;

  public:
    /// Normal constructor
    GeoTransform(GeoReferenceBase const& src_georef, GeoReferenceBase const& dst_georef) :
      m_src_georef(src_georef), m_dst_georef(dst_georef) {}

    /// Given a pixel coordinate of an image in a destination
    /// georeference frame, this routine computes the corresponding
    /// pixel from an image in the source georeference frame.
    Vector2 reverse(Vector2 const& v) const {
      return m_src_georef.lonlat_to_pixel(m_dst_georef.pixel_to_lonlat(v));
    }

    /// Given a pixel coordinate of an image in a source
    /// georeference frame, this routine computes the corresponding
    /// pixel the destination (transformed) image.
    Vector2 forward(Vector2 const& v) const {
      return m_dst_georef.lonlat_to_pixel(m_src_georef.pixel_to_lonlat(v));
    }

  };


  /// Reproject an image whose pixels contain 3D points (usually in
  /// some spherical coordinate system).  Important note: we assume
  /// that the 3D points already have the affine transform applied to
  /// them (they correspond to real 3D coordinates and not pixel
  /// coordinates in an image), therefor the affine transform portion
  /// of the georeference is completely ignored by the function.  It
  /// does not matter what affine transform you are using in the
  /// src_georef or dst_georef.
  ///
  /// Important Note: The convention here is that the Vector3 contains
  /// the ordered triple: (longitude, latitude, altitude). 
  void reproject_point_image(ImageView<Vector3> const& point_image,
                             GeoReferenceBase const& src_georef,
                             GeoReferenceBase const& dst_georef); 

}} // namespace vw::cartography

#endif // __GEO_TRANSFORM_H__
