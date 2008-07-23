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

/// \file TMSQuadTreeGenerator.h
/// 
/// A subclass of ImageQuadTreeGenerator that generates multi-resolution 
/// TMS overlays.
/// 
#ifndef __VW_MOSAIC_TMSQUADTREEGENERATOR_H__
#define __VW_MOSAIC_TMSQUADTREEGENERATOR_H__

#include <iomanip>
#include <vw/Mosaic/QuadTreeGenerator.h>

namespace vw {
namespace mosaic {

  template <class PixelT>
  class TMSQuadTreeGenerator : public vw::mosaic::ImageQuadTreeGenerator<PixelT>
  {
    typedef vw::mosaic::ImageQuadTreeGenerator<PixelT> base_type;

  protected:
    int max_lod_pixels;
    int draw_order_offset;
    std::string kml_title;
    mutable std::ostringstream root_node_tags;

  public:
    // Constructor, templatized on the source image type.  The supplied 
    // bounding box sets the size in degrees of the TMS overlay, and has 
    // nothing in particular to do with the source image or pixels.
    template <class ImageT>
    TMSQuadTreeGenerator( std::string const& tree_name,
                          ImageViewBase<ImageT> const& source )
      : vw::mosaic::ImageQuadTreeGenerator<PixelT>( tree_name, source )
    {
      base_type::set_crop_images( false );
    }
    
    virtual ~TMSQuadTreeGenerator() {}
    
    virtual std::string compute_image_path( std::string const& name ) const {
      fs::path path( base_type::m_base_dir, fs::native );

      Vector2i pos(0,0);
      for ( int i=1; i<(int)name.length(); ++i ) {
        pos *= 2;
        if( name[i]=='0' ) pos += Vector2i(0,1);
        else if( name[i]=='1' ) pos += Vector2i(1,1);
        else if( name[i]=='3' ) pos += Vector2i(1,0);
      }
      std::ostringstream oss;
      oss << name.length()-1 << "/" << pos.x() << "/" << pos.y();
      path /= oss.str();

      return path.native_file_string();
    }

    virtual void write_meta_file( typename base_type::PatchInfo const& info ) const
    {
    }

  };

} // namespace mosaic
} // namespace vw

#endif // __VW_MOSAIC_TMSQUADTREEGENERATOR_H__
