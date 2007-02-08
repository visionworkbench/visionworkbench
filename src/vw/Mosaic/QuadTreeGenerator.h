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

#include <vw/Math/BBox.h>
#include <vw/Image/ImageView.h>
#include <vw/Image/ImageViewRef.h>
#include <vw/Image/EdgeExtension.h>
#include <vw/Image/Manipulation.h>
#include <vw/Image/Algorithms.h>
#include <vw/Image/Filter.h>
#include <vw/FileIO/DiskImageResource.h>

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
        m_crop_images( true )
    {
      fs::path tree_path( tree_name, fs::native );
      m_tree_name = tree_path.leaf();
      fs::path base_path = tree_path.branch_path() / ( tree_path.leaf() + ".qtree" );
      m_base_dir = base_path.native_directory_string();
    }

    virtual ~ImageQuadTreeGenerator() {}

    virtual void write_meta_file( std::string const& name, unsigned level, BBox2i const& image_bbox, 
                                  BBox2i const& visible_bbox, BBox2i const& region_bbox ) const {
      unsigned scale = 1 << level;
      fs::path base_path( m_base_dir, fs::native );
      fs::ofstream outfile( base_path/(name+".bbx") );
      outfile << scale << "\n";
      outfile << image_bbox.min().x() << "\n";
      outfile << image_bbox.min().y() << "\n";
      outfile << image_bbox.max().x()-image_bbox.min().x() << "\n";
      outfile << image_bbox.max().y()-image_bbox.min().y() << "\n";
      outfile << visible_bbox.min().x() << "\n";
      outfile << visible_bbox.min().y() << "\n";
      outfile << visible_bbox.max().x()-visible_bbox.min().x() << "\n";
      outfile << visible_bbox.max().y()-visible_bbox.min().y() << "\n";
      outfile << m_source.cols() << "\n";
      outfile << m_source.rows() << "\n";
    }


    void generate() {
      unsigned maxdim = std::max( m_source.cols(), m_source.rows() );
      m_tree_levels = 1 + unsigned( ceilf( log( maxdim/(float)(m_patch_size-m_patch_overlap) ) / log(2.0) ) );
      m_patch_cache.resize( m_tree_levels );
      m_filename_cache.resize( m_tree_levels );

      fs::path dir_path( m_base_dir );
      fs::path top_branch_path = dir_path / m_tree_name;

      vw_out(DebugMessage) << "Using patch size: " << m_patch_size << " pixels" << std::endl;
      vw_out(DebugMessage) << "Using patch overlap: " << m_patch_overlap << " pixels" << std::endl;
      vw_out(DebugMessage) << "Generating patch files of type: " << m_output_image_type << std::endl;
      vw_out(DebugMessage) << "Generating " << top_branch_path.native_directory_string() << " quadtree with " << m_tree_levels << " levels." << std::endl;

      if( fs::exists( dir_path ) ) {
        if( ! fs::is_directory( dir_path ) )
          vw_throw( IOErr() << "Path " << dir_path.native_directory_string() << " is not a directory!  Remove it first." );
      }
      else {
        fs::create_directory( dir_path );
      }
      generate_branch( m_tree_name, m_tree_levels-1, 0, 0 );
    }

    void set_crop_bbox( BBox2i const& bbox ) {
      if( bbox.min().x() < 0 || bbox.min().y() < 0 ||
          bbox.max().x() > int(m_source.cols()) || bbox.max().y() > int(m_source.rows()) )
        vw_throw( ArgumentErr() << "Requested QuadTree bounding box exceeds source dimensions!" );
      m_crop_bbox = bbox;
    }

    void set_output_image_file_type( std::string const& extension ) {
      m_output_image_type = extension;
    }

    void set_patch_size( unsigned size ) {
      m_patch_size = size;
    }

    void set_patch_overlap( unsigned overlap ) {
      m_patch_overlap = overlap;
    }

    void set_crop_images( bool crop ) {
      m_crop_images = crop;
    }

    void set_base_dir( std::string const& base_dir ) {
      m_base_dir = base_dir;
    }
    
  protected:
    std::string m_tree_name;
    std::string m_base_dir;
    ImageViewRef<PixelT> m_source;
    BBox2i m_crop_bbox;
    std::string m_output_image_type;
    unsigned m_patch_size;
    unsigned m_patch_overlap;
    bool m_crop_images;
    unsigned m_tree_levels;
    std::vector<std::map<std::pair<unsigned,unsigned>,ImageView<PixelT> > > m_patch_cache;
    std::vector<std::map<std::pair<unsigned,unsigned>,std::string> > m_filename_cache;
    
    void write_patch( ImageView<PixelT> const& image, std::string const& name, unsigned level, unsigned x, unsigned y ) const {
      unsigned scale = 1 << level;
      unsigned interior_size = m_patch_size - m_patch_overlap;
      Vector2i position( x*interior_size-m_patch_overlap/2, y*interior_size-m_patch_overlap/2 );
      BBox2i image_bbox( Vector2i(0,0), Vector2i(image.cols(), image.rows()) );
      BBox2i visible_bbox = image_bbox, region_bbox = image_bbox;
      visible_bbox.contract( m_patch_overlap/2 );
      fs::path base_path( m_base_dir, fs::native );
      if( m_crop_images ) {
        image_bbox = nonzero_data_bounding_box( image );
        if( image_bbox.empty() ) {
          vw_out(DebugMessage) << "\tIgnoring empty image: " << name << std::endl;
          fs::remove( base_path/name );
          return;
        }
        if( image_bbox.width() == int(m_patch_size) && image_bbox.height() == int(m_patch_size) )
          write_image( (base_path/(name + "." + m_output_image_type)).native_file_string(), image );
        else
          write_image( (base_path/(name + "." + m_output_image_type)).native_file_string(), 
                       ImageView<PixelT>( crop( image, image_bbox.min().x(), image_bbox.min().y(), image_bbox.width(), image_bbox.height() ) ) );
        visible_bbox.crop( image_bbox );
      }
      else {
        write_image( (base_path/(name + "." + m_output_image_type)).native_file_string(), image );
      }
      image_bbox += position;
      visible_bbox += position;
      region_bbox += position;
      write_meta_file( name, level, image_bbox*scale, visible_bbox*scale, region_bbox*scale );
    }

    ImageView<PixelT> generate_branch( std::string name, unsigned level, unsigned x, unsigned y ) {
      ImageView<PixelT> image;
      unsigned scale = 1 << level;
      unsigned interior_size = m_patch_size - m_patch_overlap;
      BBox2i patch_bbox = scale * BBox2i( Vector2i(x, y), Vector2i(x+1, y+1) ) * interior_size;
      fs::path base_path( m_base_dir, fs::native );
      fs::path path = base_path / name;
      if( ! patch_bbox.intersects( m_crop_bbox ) ) {
        vw_out(DebugMessage) << "\tIgnoring empty image: " << path.native_file_string() << std::endl;
        image.set_size( interior_size, interior_size );
        return image;
      }
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
        fs::create_directory( path );
        ImageView<PixelT> big_image( 2*interior_size, 2*interior_size );
        crop( big_image, 0, 0, interior_size, interior_size ) = generate_branch( name + "/0", level-1, 2*x, 2*y );
        crop( big_image, interior_size, 0, interior_size, interior_size ) = generate_branch( name + "/1", level-1, 2*x+1, 2*y );
        crop( big_image, 0, interior_size, interior_size, interior_size ) = generate_branch( name + "/2", level-1, 2*x, 2*y+1 );
        crop( big_image, interior_size, interior_size, interior_size, interior_size ) = generate_branch( name + "/3", level-1, 2*x+1, 2*y+1 );
        std::vector<float> kernel(2); kernel[0]=0.5; kernel[1]=0.5;
        image.set_size( m_patch_size, m_patch_size );
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
      return image;
    }
  };

} // namespace mosaic
} // namespace vw

#endif // __VW_MOSAIC_QUADTREEGENERATOR_H__
