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
#ifndef __VW_MOSAIC_UNIVIEWQUADTREEGENERATOR_H__
#define __VW_MOSAIC_UNIVIEWQUADTREEGENERATOR_H__

#include <iomanip>
#include <vw/Mosaic/QuadTreeGenerator.h>

namespace vw {
namespace mosaic {

  template <class PixelT>
  class UniviewQuadTreeGenerator : public vw::mosaic::ImageQuadTreeGenerator<PixelT>
  {
    typedef vw::mosaic::ImageQuadTreeGenerator<PixelT> base_type;

  protected:
    int max_lod_pixels;
    int draw_order_offset;
    std::string kml_title;
    mutable std::ostringstream root_node_tags;
    bool m_terrain;

    // Override m_source to be an ImageView instead of an ImageViewRef, so
    // we can perform the pixel math without changing the original 
    // image.
    ImageView<PixelT> m_source;

  public:
    // Constructor, templatized on the source image type.  The supplied 
    // bounding box sets the size in degrees of the Uniview overlay, and has 
    // nothing in particular to do with the source image or pixels.
    template <class ImageT>
    UniviewQuadTreeGenerator( std::string const& tree_name,
                              ImageViewBase<ImageT> const& source,
                              bool terrain = false )
      : vw::mosaic::ImageQuadTreeGenerator<PixelT>( tree_name, source )
    {
      base_type::set_crop_images( false );
      m_terrain = terrain;
    }
    
    virtual ~UniviewQuadTreeGenerator() {}
    
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
      if (name.length() == 1) 
        oss << "global";
      else
        oss << name.length()-2 << "/" << pos.y() << "/" << pos.x();
      path /= oss.str();

      return path.native_file_string();
    }

    virtual void write_meta_file( typename base_type::PatchInfo const& info ) const
    {
    }

    // We have to override write_image because it is templateized on 
    // PixelT, and we need that to not be the case. The code is practically
    // identical anyway.
    template <class OutPixelT>
    void write_image( typename base_type::PatchInfo const& info, ImageView<OutPixelT> const& image ) const {
      ImageBuffer buf = image.buffer();
      DiskImageResource *r = DiskImageResource::create( info.filepath, buf.format );
      r->write( buf, BBox2i(0,0,buf.format.cols,buf.format.rows) );
      delete r;
    }

    // We override write_patch so it performs the math on the image to 
    // translate it into a signed integer (if terrain mode is turned on), 
    // before writing it out to disk.
    virtual void write_patch( ImageView<PixelT> const& image, std::string const& name, int32 level, int32 x, int32 y ) const {
      typename base_type::PatchInfo info;
      info.name = name;
      info.filepath = compute_image_path( name );
      info.level = level;
      int32 scale = 1 << level;
      int32 interior_size = this->m_patch_size - this->m_patch_overlap;
      Vector2i position( x*interior_size - this->m_patch_overlap/2, y*interior_size - this->m_patch_overlap/2 );
      info.image_bbox = BBox2i( Vector2i(0,0), Vector2i(image.cols(), image.rows()) );
      info.visible_bbox = info.region_bbox = info.image_bbox;
      info.visible_bbox.contract( this->m_patch_overlap/2 );
      ImageView<PixelT> patch_image = image;
      if( this->m_crop_images ) {
        info.image_bbox = nonzero_data_bounding_box( image );
        if( info.image_bbox.empty() ) {
          vw_out(DebugMessage, "mosaic") << "\tIgnoring empty image: " << name << std::endl;
          return;
        }
        if( info.image_bbox.width() != int(this->m_patch_size) || info.image_bbox.height() != int(this->m_patch_size) )
          patch_image = crop( image, info.image_bbox );
        info.visible_bbox.crop( info.image_bbox );
      }
      info.image_bbox = scale * (info.image_bbox + position);
      info.visible_bbox = scale * (info.visible_bbox + position);
      info.region_bbox = scale * (info.region_bbox + position);
      std::string output_image_type = this->m_output_image_type;
      if( m_terrain )
        output_image_type = "png";
      if( output_image_type == "auto" ) {
        if( is_opaque( patch_image ) ) output_image_type = "jpg";
        else output_image_type = "png";
      }
      info.filepath += "." + output_image_type;

      create_directories( fs::path( info.filepath, fs::native ).branch_path() );

      vw_out(InfoMessage, "mosaic") << "\tSaving image: " << info.filepath << "\t" << patch_image.cols() << "x" << patch_image.rows() << std::endl;
      // Transform image into 16-bit unsigned integer, if terrain.
      if( m_terrain ) {
        // This could probably be accomplished with a functor.
        // Note that this is only ever done when the pixel type of 
        // patch_image is PixelGray or PixelGrayA, so it's perfectly safe 
        // to take the lowest channel of the pixel. We can't just call the
        // v() function because PixelRGB(A) does not have that, which will
        // cause a compiler error (even though this code never gets called
        // for those pixel types, the compiler can't tell that). You can 
        // ignore the "comparison is always false" compiler warnings for
        // this reason.
        ImageView<PixelGray<uint16> > out_image;
        out_image.set_size( patch_image.cols(), patch_image.rows() );
        for(int i=0; i < patch_image.cols(); i++) {
          for(int j=0; j < patch_image.rows(); j++) {
            if(patch_image(i,j)[0] < 0)
              out_image(i,j) = PixelGray<uint16>((uint16)(patch_image(i,j)[0] + 65536));
            else
              out_image(i,j) = PixelGray<uint16>((uint16)(patch_image(i,j)[0]));
          }
        }
        write_image( info, out_image );
      } else {
        write_image( info, patch_image );
      }

      if( this->m_write_meta_file ) write_meta_file(info);
    }

  };

} // namespace mosaic
} // namespace vw

#endif // __VW_MOSAIC_UNIVIEWQUADTREEGENERATOR_H__
