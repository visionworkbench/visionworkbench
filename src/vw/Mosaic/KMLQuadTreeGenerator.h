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

/// \file KMLQuadTreeGenerator.h
/// 
/// A subclass of ImageQuadTreeGenerator that generates multi-resolution 
/// KML overlays.
/// 
#ifndef __VW_MOSAIC_KMLQUADTREEGENERATOR_H__
#define __VW_MOSAIC_KMLQUADTREEGENERATOR_H__

#include <iomanip>
#include <vw/Mosaic/QuadTreeGenerator.h>

namespace vw {
namespace mosaic {

  template <class PixelT>
  class KMLQuadTreeGenerator : public vw::mosaic::ImageQuadTreeGenerator<PixelT>
  {
    typedef vw::mosaic::ImageQuadTreeGenerator<PixelT> base_type;

  protected:
    BBox2 total_longlat_bbox;
    int max_lod_pixels;
    int draw_order_offset;
    std::string kml_title;
    std::string kml_metadata;
    mutable std::ostringstream root_node_tags;

    std::string kml_latlonaltbox( BBox2 const& longlat_bbox ) const {
      BBox2 bbox = longlat_bbox;
      bbox.crop( total_longlat_bbox );
      std::ostringstream tag;
      tag << std::setprecision(10);
      tag << "<LatLonAltBox>"
          << "<north>" << bbox.min().y() << "</north>"
          << "<south>" << bbox.max().y() << "</south>"
          << "<east>" << bbox.max().x() << "</east>"
          << "<west>" << bbox.min().x() << "</west>"
          << "</LatLonAltBox>";
      return tag.str();
    }

    std::string kml_network_link( std::string const& name, std::string const& href, BBox2 const& longlat_bbox ) const {
      std::ostringstream tag;
      tag << std::setprecision(10);
      tag << "  <NetworkLink>\n"
          << "    <name>" + name + "</name>\n"
          << "    <Region>" + kml_latlonaltbox(longlat_bbox) + "<Lod><minLodPixels>128</minLodPixels><maxLodPixels>-1</maxLodPixels></Lod></Region>\n"
          << "    <Link><href>" + href + "</href><viewRefreshMode>onRegion</viewRefreshMode></Link>\n"
          << "  </NetworkLink>\n";
      return tag.str();
    }

    std::string kml_ground_overlay( std::string const& href, BBox2 const& bbox, int draw_order, int max_lod_pixels ) const {
      std::ostringstream tag;
      tag << std::setprecision(10);
      tag << "  <GroundOverlay>\n"
          << "    <Region>" << kml_latlonaltbox(bbox) << "<Lod><minLodPixels>128</minLodPixels><maxLodPixels>" << max_lod_pixels << "</maxLodPixels></Lod></Region>\n"
          << "    <name>" << href << "</name>\n"
          << "    <Icon><href>" << href << "</href></Icon>\n"
          << "    " << kml_latlonaltbox(bbox) << "\n"
          << "    <drawOrder>" << draw_order+draw_order_offset << "</drawOrder>\n"
          << "  </GroundOverlay>\n";
      return tag.str();
    }

    BBox2 pixels_to_longlat( BBox2i const& image_bbox ) const {
      double width = this->m_source.cols(), height = this->m_source.rows();
      // Fractional bounding-box
      BBox2 fbb( image_bbox.min().x()/width, image_bbox.min().y()/height,
                 image_bbox.width()/width, image_bbox.height()/height );
      BBox2 bb( total_longlat_bbox.min().x()+fbb.min().x()*total_longlat_bbox.width(),
                total_longlat_bbox.max().y()-fbb.min().y()*total_longlat_bbox.height(),
                fbb.width()*total_longlat_bbox.width(), -fbb.height()*total_longlat_bbox.height() );
      return bb;
    }

    /* This function calculates whether we merge the top two parts of the 
     * quadtree or not, based on the bounding box we're generating a 
     * quadtree for, the current level, and the dimensions of the image.
     *
     * What the function does is calculate whether the the new aspect ratio
     * of the current center lattiude point is larger than twice the 
     * current aspect ratio.
    */
    bool calculate_merge(int32 level, BBox2i const& bbox, double center_lat) const {
      bool r;

      r = (level != this->m_tree_levels - 1);
      r = r && (1.0 / cos(M_PI / 180 * center_lat) >= 2 * bbox.width() / bbox.height());

      return r;
    }

    /* We override generate_branch in order to implement the behavior where
     * we stop using 4 branches at each node, and instead use 3 (merging 
     * the top or bottom two) as |latitude| increases.
    */
    virtual ImageView<PixelT> generate_branch( std::string name, int32 level, BBox2i const& bbox, int32 width_patch_size, const ProgressCallback &progress_callback )
    {
      progress_callback.report_progress(0);
      if(progress_callback.abort_requested())
        vw_throw(Aborted() << "Aborted by ProgressCallback" );
      ImageView<PixelT> image;
      bool wide_image;

      if(bbox.width() != bbox.height()) wide_image = true;
      else wide_image = false;
      int height = this->m_patch_size;
      int width = width_patch_size;

      // Reject patches that fall otuside the crop region
      if(!bbox.intersects( this->m_crop_bbox ) ) {
        vw_out(DebugMessage, "mosaic") << "\tIgnoring empty image: " << name << std::endl;
        image.set_size(width, height);
        return image;
      }

      // Reject patches that fail the interior intersection check.
      // This effectively prunes branches of the tree with no source data.
      if( ! (*(this->m_sparse_tile_check))(bbox) ) {
        vw_out(DebugMessage, "mosaic") << "\tIgnoring empty branch: " << name << std::endl;
        image.set_size(width, height);
        return image;
      }

      // Base case: rasterize the highest resolution tile.
      if(level == 0) {
        BBox2 data_bbox = bbox;
        data_bbox.crop(this->m_crop_bbox);
        image = crop(this->m_source, data_bbox);

        if( data_bbox != bbox )
            image = edge_extend( image, bbox - data_bbox.min(), ZeroEdgeExtension() );
      } else {
        ImageView<PixelT> big_image(2*width, 2*height);
        SubProgressCallback spcb0(progress_callback, 0., .25);
        SubProgressCallback spcb1(progress_callback, .25, .5);
        SubProgressCallback spcb2(progress_callback, .5, .75);
        SubProgressCallback spcb3(progress_callback, .75, 1.);
        SubProgressCallback spcb4(progress_callback, 0., .5);
        SubProgressCallback spcb5(progress_callback, .5, 1.);

        double center_lat = total_longlat_bbox.max().y() - bbox.center().y()*total_longlat_bbox.height() / this->m_source.rows();
        bool merge = calculate_merge(level, bbox, center_lat);

        // Top half of the quadtree.
        if(center_lat > 0 && merge) {
          // Merge the top half of this part of the quadtree.
          crop(big_image, 0, 0, 2*width_patch_size, this->m_patch_size) =
              resize( generate_branch( name+"4", level-1, BBox2i(bbox.min().x(), bbox.min().y(), bbox.width(), bbox.height()/2), 2*width_patch_size, spcb4), 2*width_patch_size, this->m_patch_size );
        } else {
          // Do not merge the top half.
          crop(big_image, 0, 0, width_patch_size, this->m_patch_size) =
              generate_branch( name+"0", level-1, BBox2i(bbox.min().x(), bbox.min().y(), bbox.width()/2, bbox.height()/2), width_patch_size, spcb0 );
          crop(big_image, width_patch_size, 0, width_patch_size, this->m_patch_size) =
              generate_branch( name+"1", level-1, BBox2i(bbox.min().x() + bbox.width()/2, bbox.min().y(), bbox.width()/2, bbox.height()/2), width_patch_size, spcb1 );
        }
        
        // Bottom half of the quadtree.
        if(center_lat < 0 && merge) {
          crop(big_image, 0, this->m_patch_size, 2*width_patch_size, this->m_patch_size) =
            resize( generate_branch( name+"5", level-1, BBox2i(bbox.min().x(), bbox.min().y()+bbox.height()/2, bbox.width(), bbox.height()/2), 2*width_patch_size, spcb4 ), 2*width_patch_size, this->m_patch_size );
        } else {
          crop(big_image, 0, this->m_patch_size, width_patch_size, this->m_patch_size) =
              generate_branch( name+"2", level-1, BBox2i(bbox.min().x(), bbox.min().y()+bbox.height()/2, bbox.width()/2, bbox.height()/2), width_patch_size, spcb2 );
          crop(big_image, width_patch_size, this->m_patch_size, width_patch_size, this->m_patch_size) =
              generate_branch( name+"3", level-1, BBox2i(bbox.min().x()+bbox.width()/2, bbox.min().y()+bbox.height()/2, bbox.width()/2, bbox.height()/2), width_patch_size, spcb3 );
        }

        std::vector<float> kernel(2);
        kernel[0] = 0.5; kernel[1] = 0.5;
        image.set_size(width, height);
        rasterize( subsample( separable_convolution_filter( big_image, kernel, kernel, 1, 1), 2), image );
      }
      write_patch(image, name, level, bbox);

      progress_callback.report_progress(1);
      return image;
    }

    // We must override write_patch too.
    void write_patch( ImageView<PixelT> const& image, std::string const& name, int32 level, BBox2i const& bbox) const {
      std::string ext;
      typename base_type::PatchInfo info;

      info.name = name;
      if( is_opaque( image ) ) ext = ".jpg";
      else ext = ".png";
      info.filepath = compute_image_path( name ) + ext;
      info.level = level;
      info.image_bbox = bbox;

      create_directories( fs::path( info.filepath, fs::native ).branch_path() );

      vw_out(InfoMessage, "mosaic") << "\tSaving image: " << info.filepath << "\t" << image.cols() << "x" << image.rows() << std::endl;
      this->write_image( info, image );
      write_meta_file( info );
    }

  public:
    // Constructor, templatized on the source image type.  The supplied 
    // bounding box sets the size in degrees of the KML overlay, and has 
    // nothing in particular to do with the source image or pixels.
    template <class ImageT>
    KMLQuadTreeGenerator( std::string const& tree_name,
                          ImageViewBase<ImageT> const& source,
                          BBox2 total_longlat_bbox )
      : vw::mosaic::ImageQuadTreeGenerator<PixelT>( tree_name, source ),
        total_longlat_bbox( total_longlat_bbox ),
        max_lod_pixels( -1 ),
        draw_order_offset( 0 )
    { }
    
    void set_max_lod_pixels( int pixels ) {
      max_lod_pixels = pixels;
    }

    void set_draw_order_offset( int offset ) {
      draw_order_offset = offset;
    }

    void set_kml_title( std::string const& name ) {
      kml_title = name;
    }

    void set_metadata( std::string const& data ) {
      kml_metadata = data;
    }

    virtual ~KMLQuadTreeGenerator() {}
    
    virtual std::string compute_image_path( std::string const& name ) const {
      fs::path path( base_type::m_base_dir, fs::native );

      if( name.length() == 1 ) {
        path /= change_extension( path, "" ).leaf();
      }
      else {
        for ( int i=1; i<(int)name.length() - (int)base_type::m_levels_per_directory; i+=base_type::m_levels_per_directory ) {
          path /= name.substr( i, base_type::m_levels_per_directory );
        }
        path /= name.substr( 1, name.length()-1 );
      }

      return path.native_file_string();
    }

    virtual void write_meta_file( typename base_type::PatchInfo const& info ) const
    {
      BBox2 bbox = pixels_to_longlat( this->m_crop_bbox );
      fs::path file_path( info.filepath, fs::native );
      int base_len = file_path.branch_path().native_file_string().size() + 1;

      fs::path kml_path = change_extension( file_path, ".kml" );
      fs::ofstream kml( kml_path );
      kml << std::setprecision(10);

      kml << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
      kml << "<kml><Folder>\n";
      int children = 0;

      bool root_node = ( info.level == this->m_tree_levels-1 );

      if( root_node ) {
        if( ! kml_title.empty() )
          root_node_tags << "  <name>" << kml_title << "</name>\n";
        if( ! kml_metadata.empty() )
          root_node_tags << "  <Metadata>" << kml_metadata << "</Metadata>\n";
        double lon = (bbox.min().x()+bbox.max().x())/2;
        double lat = (bbox.min().y()+bbox.max().y())/2;
        double range = 1e5 * (bbox.width()*cos(M_PI/180*lat)-bbox.height());
        if( range > 1.2e7 ) range = 1.2e7;
        root_node_tags << "  <LookAt><longitude>" << lon << "</longitude><latitude>" << lat << "</latitude><range>" << range << "</range></LookAt>\n";
        root_node_tags << "  <Style><ListStyle><listItemType>checkHideChildren</listItemType></ListStyle></Style>\n";
      }

      double center_lat = total_longlat_bbox.max().y() - info.image_bbox.center().y() * total_longlat_bbox.height()/this->m_source.rows();
      bool merge = calculate_merge(info.level, info.image_bbox, center_lat);

      BBox2 bb = pixels_to_longlat( info.image_bbox );

      // Compute the top half.
      if(center_lat > 0 && merge) {
        // We call our merged top half "kml4".
        std::string kml4 = compute_image_path( info.name+"4" ) + ".kml";
        if( exists( fs::path( kml4, fs::native ) ) ) {
          children++;
          kml << kml_network_link( "4", kml4.substr(base_len), BBox2(bb.min().x(), bb.min().y(), bb.width(), bb.height()/2) );
        }
      } else {
        // Normal behavior for the top half.
        std::string kml0 = compute_image_path( info.name+"0" ) + ".kml";
        if( exists( fs::path( kml0, fs::native ) ) ) {
          children++;
          kml << kml_network_link( "0", kml0.substr(base_len), BBox2(bb.min().x(), bb.min().y(), bb.width()/2, bb.height()/2) );
        }
        std::string kml1 = compute_image_path( info.name+"1" ) + ".kml";
        if( exists( fs::path( kml1, fs::native ) ) ) {
          children++;
          kml << kml_network_link( "1", kml1.substr(base_len), BBox2(bb.min().x()+bb.width()/2, bb.min().y(), bb.width()/2, bb.height()/2) );
        }
      }

      // Now compute the bottom half.
      if(center_lat < 0 && merge) {
        // We use 'kml5' to represent a merged bottom half of the quadtree.
        std::string kml5 = compute_image_path( info.name+"5" ) + ".kml";
        if( exists( fs::path( kml5, fs::native ) ) ) {
          children++;
          kml << kml_network_link( "5", kml5.substr(base_len), BBox2(bb.min().x(), bb.min().y()+bb.height()/2, bb.width(), bb.height()/2) );
        } 
      } else {
        std::string kml2 = compute_image_path( info.name+"2" ) + ".kml";
        if( exists( fs::path( kml2, fs::native ) ) ) {
          children++;
          kml << kml_network_link( "2", kml2.substr(base_len), BBox2(bb.min().x(), bb.min().y()+bb.height()/2, bb.width()/2, bb.height()/2) );
        }
        std::string kml3 = compute_image_path( info.name+"3" ) + ".kml";
        if( exists( fs::path( kml3, fs::native ) ) ) {
          children++;
          kml << kml_network_link( "3", kml3.substr(base_len), BBox2(bb.min().x()+bb.width()/2, bb.min().y()+bb.height()/2, bb.width()/2, bb.height()/2 ) );
        }
      }

      int max_lod = max_lod_pixels;
      if(!children) max_lod = -1;
      kml << kml_ground_overlay( file_path.leaf(), bb, this->m_tree_levels - info.level, max_lod );

      if( root_node ) kml << root_node_tags.str();
      kml << "</Folder></kml>\n";
    }

    // Overriding generate as well.
    void generate( const ProgressCallback &progress_callback = ProgressCallback::dummy_instance() ) {
      int32 maxdim = std::max(this->m_source.cols(), this->m_source.rows());
      this->m_tree_levels = 1 + int32( ceil( log( maxdim/(double)(this->m_patch_size) ) / log(2.0) ) );

      vw_out(DebugMessage, "mosaic") << "Using patch size: " << this->m_patch_size << " pixels" << std::endl;
      vw_out(DebugMessage, "mosaic") << "Generating " << this->m_base_dir << " quadtree with " << this->m_tree_levels << " levels." << std::endl;

      int tree_size = this->m_patch_size * (1 << (this->m_tree_levels-1));
      generate_branch("r", this->m_tree_levels-1, BBox2i(0,0,tree_size,tree_size), this->m_patch_size, progress_callback);
      progress_callback.report_finished();
    }

  };


} // namespace mosaic
} // namespace vw

#endif // __VW_MOSAIC_KMLQUADTREEGENERATOR_H__
