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

/// \file QuadTreeGenerator.h
/// 
/// A class that generates filesystem-based quadtrees of large images.
/// 
#ifndef __VW_MOSAIC_QUADTREEGENERATOR_H__
#define __VW_MOSAIC_QUADTREEGENERATOR_H__

#include <vector>
#include <map>
#include <string>
#include <fstream>

#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/filesystem/convenience.hpp>

#include <vw/Math/BBox.h>
#include <vw/Image/ImageView.h>
#include <vw/Image/ImageViewRef.h>
#include <vw/Image/EdgeExtension.h>
#include <vw/Image/Manipulation.h>
#include <vw/Image/Algorithms.h>
#include <vw/Image/Filter.h>
#include <vw/FileIO/DiskImageResource.h>
#include <vw/Core/ProgressCallback.h>
#include <vw/Mosaic/SparseTileCheck.h>

namespace vw {
namespace mosaic {

  namespace fs = boost::filesystem;

  template <class PixelT>
  class ImageQuadTreeGenerator {
  public:

    template <class ImageT>
    ImageQuadTreeGenerator( std::string const& tree_name, ImageViewBase<ImageT> const& source )
      : m_source( source.impl() ),
        m_crop_bbox( 0, 0, m_source.cols(), m_source.rows() ),
        m_output_image_type( "png" ),
        m_patch_size( 256 ),
        m_patch_overlap( 0 ),
        m_levels_per_directory( 3 ),
        m_crop_images( true )
    {
      fs::path tree_path( tree_name, fs::native );
      fs::path base_path = tree_path.branch_path() / tree_path.leaf();
      m_base_dir = base_path.native_directory_string();
      m_sparse_tile_check = boost::shared_ptr<SparseTileCheckBase>(new SparseTileCheck<ImageT>(source.impl()));
    }

    virtual ~ImageQuadTreeGenerator() {}

    struct PatchInfo {
      std::string name, filepath;
      int32 level;
      BBox2i image_bbox;
      BBox2i visible_bbox;
      BBox2i region_bbox;
    };

    virtual std::string compute_image_path( std::string const& name ) const {
      fs::path path( m_base_dir, fs::native );
      
      for (int i= 0; i< (int)name.length() - (int)m_levels_per_directory; i += m_levels_per_directory) {
        path /= name.substr( i, m_levels_per_directory );
      }
      path /= name;

      return path.native_file_string();
    }

    virtual void write_image( PatchInfo const& info, ImageView<PixelT> const& image ) const {
      ImageBuffer buf = image.buffer();
      DiskImageResource *r = DiskImageResource::create( info.filepath, buf.format );
      r->write( buf, BBox2i(0,0,buf.format.cols,buf.format.rows) );
      delete r;
    }

    virtual void write_meta_file( PatchInfo const& info ) const {
      int32 scale = 1 << info.level;
      fs::path filepath( info.filepath, fs::native );
      fs::ofstream outfile( change_extension( filepath, ".bbx" ) );
      outfile << scale << "\n";
      outfile << info.image_bbox.min().x() << "\n";
      outfile << info.image_bbox.min().y() << "\n";
      outfile << info.image_bbox.width() << "\n";
      outfile << info.image_bbox.height() << "\n";
      outfile << info.visible_bbox.min().x() << "\n";
      outfile << info.visible_bbox.min().y() << "\n";
      outfile << info.visible_bbox.width() << "\n";
      outfile << info.visible_bbox.height() << "\n";
      outfile << m_source.cols() << "\n";
      outfile << m_source.rows() << "\n";
    }

    void generate( const ProgressCallback &progress_callback = ProgressCallback::dummy_instance() ) {
      int32 maxdim = std::max( m_source.cols(), m_source.rows() );
      m_tree_levels = 1 + int32( ceilf( log( maxdim/(float)(m_patch_size-m_patch_overlap) ) / log(2.0) ) );
      m_patch_cache.resize( m_tree_levels );
      m_filename_cache.resize( m_tree_levels );

      vw_out(DebugMessage) << "Using patch size: " << m_patch_size << " pixels" << std::endl;
      vw_out(DebugMessage) << "Using patch overlap: " << m_patch_overlap << " pixels" << std::endl;
      vw_out(DebugMessage) << "Generating patch files of type: " << m_output_image_type << std::endl;
      vw_out(DebugMessage) << "Generating " << m_base_dir << " quadtree with " << m_tree_levels << " levels." << std::endl;

      try {
        generate_branch( "r", m_tree_levels-1, 0, 0, progress_callback );
        progress_callback.report_finished();
      } catch (Aborted) {
        progress_callback.report_aborted();
      }
    }

    BBox2i const& get_crop_bbox() const {
      return m_crop_bbox;
    }

    void set_crop_bbox( BBox2i const& bbox ) {
      if( bbox.min().x() < 0 || bbox.min().y() < 0 ||
          bbox.max().x() > int(m_source.cols()) || bbox.max().y() > int(m_source.rows()) )
        vw_throw( ArgumentErr() << "Requested QuadTree bounding box exceeds source dimensions!" );
      m_crop_bbox = bbox;
    }

    std::string const& get_output_image_file_type() const {
      return m_output_image_type;
    }

    void set_output_image_file_type( std::string const& extension ) {
      m_output_image_type = extension;
    }

    int32 get_patch_size() const {
      return m_patch_size;
    }

    void set_patch_size( int32 size ) {
      m_patch_size = size;
    }

    int32 get_patch_overlap() const {
      return m_patch_overlap;
    }

    void set_patch_overlap( int32 overlap ) {
      m_patch_overlap = overlap;
    }

    int32 get_levels_per_directory() const {
      return m_levels_per_directory;
    }

    void set_levels_per_directory( int32 levels_per_directory ) {
      m_levels_per_directory = levels_per_directory;
    }

    bool get_crop_images() const {
      return m_crop_images;
    }

    void set_crop_images( bool crop ) {
      m_crop_images = crop;
    }

    std::string const& get_base_dir() const {
      return m_base_dir;
    }

    void set_base_dir( std::string const& base_dir ) {
      m_base_dir = base_dir;
    }
    
  protected:
    std::string m_base_dir;
    ImageViewRef<PixelT> m_source;
    BBox2i m_crop_bbox;
    std::string m_output_image_type;
    int32 m_patch_size;
    int32 m_patch_overlap;
    int32 m_levels_per_directory;
    bool m_crop_images;
    int32 m_tree_levels;
    std::vector<std::map<std::pair<int32,int32>,ImageView<PixelT> > > m_patch_cache;
    std::vector<std::map<std::pair<int32,int32>,std::string> > m_filename_cache;
    boost::shared_ptr<SparseTileCheckBase> m_sparse_tile_check;
    
    void write_patch( ImageView<PixelT> const& image, std::string const& name, int32 level, int32 x, int32 y ) const {
      PatchInfo info;
      info.name = name;
      info.filepath = compute_image_path( name );
      info.level = level;
      int32 scale = 1 << level;
      int32 interior_size = m_patch_size - m_patch_overlap;
      Vector2i position( x*interior_size-m_patch_overlap/2, y*interior_size-m_patch_overlap/2 );
      info.image_bbox = BBox2i( Vector2i(0,0), Vector2i(image.cols(), image.rows()) );
      info.visible_bbox = info.region_bbox = info.image_bbox;
      info.visible_bbox.contract( m_patch_overlap/2 );
      ImageView<PixelT> patch_image = image;
      if( m_crop_images ) {
        info.image_bbox = nonzero_data_bounding_box( image );
        if( info.image_bbox.empty() ) {
          vw_out(DebugMessage) << "\tIgnoring empty image: " << name << std::endl;
          return;
        }
        if( info.image_bbox.width() != int(m_patch_size) || info.image_bbox.height() != int(m_patch_size) )
          patch_image = crop( image, info.image_bbox );
        info.visible_bbox.crop( info.image_bbox );
      }
      info.image_bbox = scale * (info.image_bbox + position);
      info.visible_bbox = scale * (info.visible_bbox + position);
      info.region_bbox = scale * (info.region_bbox + position);
      std::string output_image_type = m_output_image_type;
      if( output_image_type == "auto" ) {
        if( is_opaque( patch_image ) ) output_image_type = "jpg";
        else output_image_type = "png";
      }
      info.filepath += "." + output_image_type;

      create_directories( fs::path( info.filepath, fs::native ).branch_path() );

      vw_out(InfoMessage+1) << "\tSaving image: " << info.filepath << "\t" << patch_image.cols() << "x" << patch_image.rows() << std::endl;
      write_image( info, patch_image );
      write_meta_file( info );
    }

    ImageView<PixelT> generate_branch( std::string name, int32 level, int32 x, int32 y, 
                                       const ProgressCallback &progress_callback ) 
    {
      progress_callback.report_progress(0);
      if (progress_callback.abort_requested()) {
        vw_throw( Aborted() << "Aborted by ProgressCallback" );
      }
      ImageView<PixelT> image;
      int32 scale = 1 << level;
      int32 interior_size = m_patch_size - m_patch_overlap;
      BBox2i patch_bbox = scale * BBox2i( Vector2i(x, y), Vector2i(x+1, y+1) ) * interior_size;

      // Reject patches that fall outside the crop region
      if( ! patch_bbox.intersects( m_crop_bbox ) ) {
        vw_out(DebugMessage) << "\tIgnoring empty image: " << name << std::endl;
        image.set_size( interior_size, interior_size );
        return image;
      }

      // Reject patches that fail the interior intersection check.
      // This effectively prunes branches of the tree with no source
      // data.
      if( ! (*m_sparse_tile_check)(patch_bbox) ) {
        vw_out(DebugMessage) << "\tIgnoring empty branch: " << name << std::endl;
        image.set_size( interior_size, interior_size );
        return image;
      }

      // Base case: rasterize the highest resolution tile.
      if( level == 0 ) {
        vw_out(DebugMessage) << "ImageQuadTreeGenerator rasterizing region " << patch_bbox << std::endl;
        BBox2i data_bbox = patch_bbox;
        data_bbox.crop( m_crop_bbox );
        image = crop( m_source, data_bbox );
        if( data_bbox != patch_bbox ) {
          image = edge_extend( image, patch_bbox-data_bbox.min(), ZeroEdgeExtension() );
        }
      }
      else {
        ImageView<PixelT> big_image( 2*interior_size, 2*interior_size );
        SubProgressCallback spcb0( progress_callback,  0., .25 );
        SubProgressCallback spcb1( progress_callback, .25, .5 );
        SubProgressCallback spcb2( progress_callback, .5 , .75 );
        SubProgressCallback spcb3( progress_callback, .75, 1.);

        crop( big_image, 0, 0, interior_size, interior_size ) = generate_branch( name + "0", level-1, 2*x, 2*y, spcb0 );
        crop( big_image, interior_size, 0, interior_size, interior_size ) = generate_branch( name + "1", level-1, 2*x+1, 2*y, spcb1 );
        crop( big_image, 0, interior_size, interior_size, interior_size ) = generate_branch( name + "2", level-1, 2*x, 2*y+1, spcb2 );
        crop( big_image, interior_size, interior_size, interior_size, interior_size ) = generate_branch( name + "3", level-1, 2*x+1, 2*y+1, spcb3 );
        std::vector<float> kernel(2); kernel[0]=0.5; kernel[1]=0.5;
        image.set_size( interior_size, interior_size );
        rasterize( subsample( separable_convolution_filter( big_image, kernel, kernel, 1, 1 ), 2 ), image );
      }
      // If there's no patch overlap, we're done and we can just write it out
      if( m_patch_overlap == 0 ) {
        write_patch( image, name, level, x, y );
      }
      // Otherwise this interior affects up to nine patches, each of which 
      // requires some special consideration.  There may be a cleaner way 
      // to structure this, but this works for now.
      else {
        int overlap = m_patch_overlap / 2;
        if( x > 0 ) {
          if( y > 0 ) {
            // Top left
            ImageView<PixelT> &top_left = m_patch_cache[level][std::make_pair(x-1,y-1)];
            crop( top_left, interior_size+overlap, interior_size+overlap, overlap, overlap ) = crop( image, 0, 0, overlap, overlap );
            write_patch( top_left, m_filename_cache[level][std::make_pair(x-1,y-1)], level, x-1, y-1 );
            m_patch_cache[level].erase(std::make_pair(x-1,y-1));
            m_filename_cache[level].erase(std::make_pair(x-1,y-1));
          }
          // Left
          ImageView<PixelT> &left = m_patch_cache[level][std::make_pair(x-1,y)];
          crop( left, interior_size+overlap, overlap, overlap, interior_size ) = crop( image, 0, 0, overlap, interior_size );
          if( scale*interior_size*(y+1) < m_source.rows() ) {
            // Bottom left
            ImageView<PixelT> &bot_left = m_patch_cache[level][std::make_pair(x-1,y+1)];
            bot_left.set_size( m_patch_size, m_patch_size );
            crop( bot_left, interior_size+overlap, 0, overlap, overlap ) = crop( image, 0, interior_size-overlap, overlap, overlap );
          }
          else {
            write_patch( left, m_filename_cache[level][std::make_pair(x-1,y)], level, x-1, y );
            m_patch_cache[level].erase(std::make_pair(x-1,y));
            m_filename_cache[level].erase(std::make_pair(x-1,y));
          }
        }
        if( y > 0 ) {
          // Top
          ImageView<PixelT> &top = m_patch_cache[level][std::make_pair(x,y-1)];
          crop( top, overlap, interior_size+overlap, interior_size, overlap ) = crop( image, 0, 0, interior_size, overlap );
          if( !( scale*interior_size*(x+1) < m_source.cols() ) ) {
            write_patch( top, m_filename_cache[level][std::make_pair(x,y-1)], level, x, y-1 );
            m_patch_cache[level].erase(std::make_pair(x,y-1));
            m_filename_cache[level].erase(std::make_pair(x,y-1));
          }
        }
        // Center
        ImageView<PixelT> &center = m_patch_cache[level][std::make_pair(x,y)];
        center.set_size( m_patch_size, m_patch_size );
        crop( center, overlap, overlap, interior_size, interior_size ) = image;
        if( scale*interior_size*(x+1) < m_source.cols() || scale*interior_size*(y+1) < m_source.rows() ) {
          m_filename_cache[level][std::make_pair(x,y)] = name;
        }
        else {
          write_patch( center, name, level, x, y );
          m_patch_cache[level].erase(std::make_pair(x,y-1));
        }
        if( scale*interior_size*(y+1) < m_source.rows() ) {
          // Bottom
          ImageView<PixelT> &bot = m_patch_cache[level][std::make_pair(x,y+1)];
          bot.set_size( m_patch_size, m_patch_size );
          crop( bot, overlap, 0, interior_size, overlap ) = crop( image, 0, interior_size-overlap, interior_size, overlap );
        }
        if( scale*interior_size*(x+1) < m_source.cols() ) {
          if( y > 0 ) {
            // Top right
            ImageView<PixelT> &top_right = m_patch_cache[level][std::make_pair(x+1,y-1)];
            top_right.set_size( m_patch_size, m_patch_size );
            crop( top_right, 0, interior_size+overlap, overlap, overlap ) = crop( image, interior_size-overlap, 0, overlap, overlap );
          }
          // Right
          ImageView<PixelT> &right = m_patch_cache[level][std::make_pair(x+1,y)];
          right.set_size( m_patch_size, m_patch_size );
          crop( right, 0, overlap, overlap, interior_size ) = crop( image, interior_size-overlap, 0, overlap, interior_size );
          if( scale*interior_size*(y+1) < m_source.rows() ) {
            // Bottom right
            ImageView<PixelT> &bot_right = m_patch_cache[level][std::make_pair(x+1,y+1)];
            bot_right.set_size( m_patch_size, m_patch_size );
            crop( bot_right, 0, 0, overlap, overlap ) = crop( image, interior_size-overlap, interior_size-overlap, overlap, overlap );
          }
        }
      }

      progress_callback.report_progress(1);
      return image;
    }
  };

} // namespace mosaic
} // namespace vw

#endif // __VW_MOSAIC_QUADTREEGENERATOR_H__
