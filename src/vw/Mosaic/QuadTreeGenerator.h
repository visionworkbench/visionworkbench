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

#include <vw/Math/BBox.h>
#include <vw/Image/ImageView.h>
#include <vw/Image/EdgeExtend.h>
#include <vw/Image/Manipulation.h>
#include <vw/Image/Algorithms.h>
#include <vw/Image/Filter.h>
#include <vw/FileIO/DiskImageResource.h>

namespace vw {
namespace mosaic {

  namespace fs = boost::filesystem;

  void write_bbox_file( std::string const& name, int scale, BBox2i const& image_bbox, 
                        BBox2i const& visible_bbox, int full_x_dim, int full_y_dim ) {
    std::ofstream outfile( name.c_str() );
    outfile << scale << "\n";
    outfile << image_bbox.min().x() << "\n";
    outfile << image_bbox.min().y() << "\n";
    outfile << image_bbox.max().x()-image_bbox.min().x() << "\n";
    outfile << image_bbox.max().y()-image_bbox.min().y() << "\n";
    outfile << visible_bbox.min().x() << "\n";
    outfile << visible_bbox.min().y() << "\n";
    outfile << visible_bbox.max().x()-visible_bbox.min().x() << "\n";
    outfile << visible_bbox.max().y()-visible_bbox.min().y() << "\n";
    outfile << full_x_dim << "\n";
    outfile << full_y_dim << "\n";
  }

  template <class PixelT>
  class ImageQuadTreeGenerator {
  public:

    struct Settings {
      std::string tree_name;
      std::string output_image_type;
      int patch_size;
      int patch_overlap;
      bool crop_images;
    };

    class PatchSource {
    public:
      virtual ~PatchSource() {}
      virtual ImageView<PixelT> generate_patch( BBox2i const& bbox ) const = 0;
      virtual int cols() const = 0;
      virtual int rows() const = 0;
      virtual bool intersects( BBox2i const& bbox ) const {
        if( bbox.max().x()<=0 || bbox.max().y()<=0 || bbox.min().x() >= cols() || bbox.min().y() >= rows() ) return false;
        else return true;
      }
    };

    ImageQuadTreeGenerator( Settings const& settings, PatchSource const& source )
      : settings( settings ), source( source ) {}
    
    void generate() {
      int maxdim = std::max( source.cols(), source.rows() );
      int tree_levels = 1 + int( ceilf( log( maxdim/(float)(settings.patch_size-settings.patch_overlap) ) / log(2.0) ) );
      patch_cache.resize( tree_levels );
      filename_cache.resize( tree_levels );

      fs::path tree_path( settings.tree_name, fs::native );
      fs::path dir_path = tree_path.branch_path() / ( tree_path.leaf() + ".qtree" );
      fs::path top_branch_path = dir_path / tree_path.leaf();

      std::cout << "Using patch size: " << settings.patch_size << " pixels" << std::endl;
      std::cout << "Using patch overlap: " << settings.patch_overlap << " pixels" << std::endl;
      std::cout << "Generating patch files of type: " << settings.output_image_type << std::endl;
      std::cout << "Generating " << dir_path.native_file_string() << " quadtree with " << tree_levels << " levels." << std::endl;

      if( fs::exists( dir_path ) )
        throw IOErr() << "Path " << dir_path.native_file_string() << " already exists!  Remove it first.";
      fs::create_directory( dir_path );
      generate_branch( top_branch_path, tree_levels-1, 0, 0 );
    }

    
  private:
    Settings settings;
    PatchSource const& source;
    std::vector<std::map<std::pair<int,int>,ImageView<PixelT> > > patch_cache;
    std::vector<std::map<std::pair<int,int>,std::string> > filename_cache;
    
    void write_patch( ImageView<PixelT> const& image, std::string const& name, int level, int x, int y ) const {
      int scale = 1 << level;
      int interior_size = settings.patch_size - settings.patch_overlap;
      Vector<int,2> position( x*interior_size-settings.patch_overlap/2, y*interior_size-settings.patch_overlap/2 );
      BBox2i image_bbox( Vector<int,2>(0,0), Vector<int,2>(image.cols(), image.rows()) );
      BBox2i visible_bbox = image_bbox;
      visible_bbox.contract( settings.patch_overlap/2 );
      if( settings.crop_images ) {
        image_bbox = bounding_box( image );
        if( image_bbox.empty() ) {
          std::cout << "\tIgnoring empty image: " + name + ".png\n";
          fs::remove( fs::path( name, fs::native ) );
          return;
        }
        if( image_bbox.width() == int(settings.patch_size) && image_bbox.height() == int(settings.patch_size) )
          write_image( name + "." + settings.output_image_type, image );
        else
          write_image( name + "." + settings.output_image_type, ImageView<PixelT>( crop( image, image_bbox.min().x(), image_bbox.min().y(), image_bbox.width(), image_bbox.height() ) ) );
        visible_bbox.crop( image_bbox );
      }
      else {
        write_image( name + "." + settings.output_image_type, image );
      }
      image_bbox += position;
      visible_bbox += position;
      write_bbox_file( name + ".bbx", scale, image_bbox*scale, visible_bbox*scale, source.cols(), source.rows() );
    }

    ImageView<PixelT> generate_branch( fs::path const& path, int level, int x, int y ) {
      ImageView<PixelT> image;
      int scale = 1 << level;
      int interior_size = settings.patch_size - settings.patch_overlap;
      BBox2i patch_bbox = scale * BBox2i( Vector<int,2>(x, y), Vector<int,2>(x+1, y+1) ) * interior_size;
      if( ! source.intersects( patch_bbox ) ) {
        std::cout << "\tIgnoring empty image: " + path.native_file_string() + ".png\n";
        image.set_size( interior_size, interior_size );
        return image;
      }
      if( level == 0 ) {
        image = source.generate_patch( patch_bbox );
      }
      else {
        fs::create_directory( path );
        ImageView<PixelT> big_image( 2*interior_size, 2*interior_size );
        crop( big_image, 0, 0, interior_size, interior_size ) = generate_branch( path / "0", level-1, 2*x, 2*y );
        crop( big_image, interior_size, 0, interior_size, interior_size ) = generate_branch( path / "1", level-1, 2*x+1, 2*y );
        crop( big_image, 0, interior_size, interior_size, interior_size ) = generate_branch( path / "2", level-1, 2*x, 2*y+1 );
        crop( big_image, interior_size, interior_size, interior_size, interior_size ) = generate_branch( path / "3", level-1, 2*x+1, 2*y+1 );
        std::vector<float> kernel(2); kernel[0]=0.5; kernel[1]=0.5;
        image.set_size( settings.patch_size, settings.patch_size );
        rasterize( subsample( separable_convolution_filter( big_image, kernel, kernel, 1, 1 ), 2 ), image );
      }
      // If there's no patch overlap, we're done and we can just write it out
      if( settings.patch_overlap == 0 ) {
        write_patch( image, path.native_file_string(), level, x, y );
      }
      // Otherwise this interior affects up to nine patches, each of which 
      // requires some special consideration.  There may be a cleaner way 
      // to structure this, but this works for now.
      else {
        int overlap = settings.patch_overlap / 2;
        if( x > 0 ) {
          if( y > 0 ) {
            // Top left
            ImageView<PixelT> &top_left = patch_cache[level][std::make_pair(x-1,y-1)];
            crop( top_left, interior_size+overlap, interior_size+overlap, overlap, overlap ) = crop( image, 0, 0, overlap, overlap );
            write_patch( top_left, filename_cache[level][std::make_pair(x-1,y-1)], level, x-1, y-1 );
            patch_cache[level].erase(std::make_pair(x-1,y-1));
            filename_cache[level].erase(std::make_pair(x-1,y-1));
          }
          // Left
          ImageView<PixelT> &left = patch_cache[level][std::make_pair(x-1,y)];
          crop( left, interior_size+overlap, overlap, overlap, interior_size ) = crop( image, 0, 0, overlap, interior_size );
          if( scale*interior_size*(y+1) < source.rows() ) {
            // Bottom left
            ImageView<PixelT> &bot_left = patch_cache[level][std::make_pair(x-1,y+1)];
            bot_left.set_size( settings.patch_size, settings.patch_size );
            crop( bot_left, interior_size+overlap, 0, overlap, overlap ) = crop( image, 0, interior_size-overlap, overlap, overlap );
          }
          else {
            write_patch( left, filename_cache[level][std::make_pair(x-1,y)], level, x-1, y );
            patch_cache[level].erase(std::make_pair(x-1,y));
            filename_cache[level].erase(std::make_pair(x-1,y));
          }
        }
        if( y > 0 ) {
          // Top
          ImageView<PixelT> &top = patch_cache[level][std::make_pair(x,y-1)];
          crop( top, overlap, interior_size+overlap, interior_size, overlap ) = crop( image, 0, 0, interior_size, overlap );
          if( !( scale*interior_size*(x+1) < source.cols() ) ) {
            write_patch( top, filename_cache[level][std::make_pair(x,y-1)], level, x, y-1 );
            patch_cache[level].erase(std::make_pair(x,y-1));
            filename_cache[level].erase(std::make_pair(x,y-1));
          }
        }
        // Center
        ImageView<PixelT> &center = patch_cache[level][std::make_pair(x,y)];
        center.set_size( settings.patch_size, settings.patch_size );
        crop( center, overlap, overlap, interior_size, interior_size ) = image;
        if( scale*interior_size*(x+1) < source.cols() || scale*interior_size*(y+1) < source.rows() ) {
          filename_cache[level][std::make_pair(x,y)] = path.native_file_string();
        }
        else {
          write_patch( center, path.native_file_string(), level, x, y );
          patch_cache[level].erase(std::make_pair(x,y-1));
        }
        if( scale*interior_size*(y+1) < source.rows() ) {
          // Bottom
          ImageView<PixelT> &bot = patch_cache[level][std::make_pair(x,y+1)];
          bot.set_size( settings.patch_size, settings.patch_size );
          crop( bot, overlap, 0, interior_size, overlap ) = crop( image, 0, interior_size-overlap, interior_size, overlap );
        }
        if( scale*interior_size*(x+1) < source.cols() ) {
          if( y > 0 ) {
            // Top right
            ImageView<PixelT> &top_right = patch_cache[level][std::make_pair(x+1,y-1)];
            top_right.set_size( settings.patch_size, settings.patch_size );
            crop( top_right, 0, interior_size+overlap, overlap, overlap ) = crop( image, interior_size-overlap, 0, overlap, overlap );
          }
          // Right
          ImageView<PixelT> &right = patch_cache[level][std::make_pair(x+1,y)];
          right.set_size( settings.patch_size, settings.patch_size );
          crop( right, 0, overlap, overlap, interior_size ) = crop( image, interior_size-overlap, 0, overlap, interior_size );
          if( scale*interior_size*(y+1) < source.rows() ) {
            // Bottom right
            ImageView<PixelT> &bot_right = patch_cache[level][std::make_pair(x+1,y+1)];
            bot_right.set_size( settings.patch_size, settings.patch_size );
            crop( bot_right, 0, 0, overlap, overlap ) = crop( image, interior_size-overlap, interior_size-overlap, overlap, overlap );
          }
        }
      }
      return image;
    }
  };


  template <class ImageT>
  class ImageViewPatchSource : public ImageQuadTreeGenerator<typename ImageT::pixel_type>::PatchSource {
    ImageT image;
  public:
    typedef typename ImageT::pixel_type pixel_type;
    
    ImageViewPatchSource( ImageT const& image ) : image( image ) {}
    
    virtual ~ImageViewPatchSource() {}

    ImageView<pixel_type> generate_patch( BBox2i const& bbox ) const {
      BBox2i image_bbox(0,0,image.cols(),image.rows());
      image_bbox.crop( bbox );
      ImageView<pixel_type> result( bbox.width(), bbox.height(), image.planes() );
      crop( result, image_bbox-bbox.min() ) = crop( image, image_bbox );
      return result;
    }
    
    virtual int cols() const {
      return image.cols();
    }

    virtual int rows() const {
      return image.rows();
    }
  };

} // namespace mosaic
} // namespace vw

#endif // __VW_MOSAIC_QUADTREEGENERATOR_H__
