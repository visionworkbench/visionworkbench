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

/// \file KMLQuadTreeConfig.h
/// 
/// A configuration class that provides callbacks for
/// QuadTreeGenerator that generate multi-resolution
/// KML overlays.
/// 
#ifndef __VW_MOSAIC_KMLQUADTREECONFIG_H__
#define __VW_MOSAIC_KMLQUADTREECONFIG_H__

#include <vw/Mosaic/QuadTreeGenerator.h>

namespace vw {
namespace mosaic {

  struct KMLQuadTreeConfigData;

  class KMLQuadTreeConfig {
  public:
    KMLQuadTreeConfig();

    // Set the extents (in degrees) of the quadtree.
    void set_longlat_bbox( BBox2 const& bbox );

    // Set the title of the root KML file.
    void set_title( std::string const& title );

    // Set the maximum level of detail (in pixels) at which each resolution  
    // of the quadtree is displayed.
    void set_max_lod_pixels( int32 pixels );

    // Set the drawOrder offset.  Overlays with positive offets are drawn on top.
    void set_draw_order_offset( int32 offset );

    // Set an option string of additional metadata to be included in the root KML file.
    void set_metadata( std::string const& data );
    
    // Configure the given quadtree to generate this KML.  This enables image 
    // cropping and sets the image path function, branch function, and metadata 
    // function.  If you intend to override any of these, be sure to do so 
    // *after* calling configure() or your changes will be overwritten.
    void configure( QuadTreeGenerator& qtree ) const;

  private:
    // The implementation is stored in a shared pointer so that it can 
    // be safely bound to the quadtree callbacks in colsures even if 
    // this config object goes out of scope.
    boost::shared_ptr<KMLQuadTreeConfigData> m_data;

    // This object is non-copyable for now.
    KMLQuadTreeConfig(KMLQuadTreeConfig const&);
    KMLQuadTreeConfig& operator=(KMLQuadTreeConfig const&);
  };

} // namespace mosaic
} // namespace vw

#endif // __VW_MOSAIC_KMLQUADTREECONFIG_H__
