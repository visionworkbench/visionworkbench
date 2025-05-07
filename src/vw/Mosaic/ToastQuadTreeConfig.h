// __BEGIN_LICENSE__
//  Copyright (c) 2006-2013, United States Government as represented by the
//  Administrator of the National Aeronautics and Space Administration. All
//  rights reserved.
//
//  The NASA Vision Workbench is licensed under the Apache License,
//  Version 2.0 (the "License"); you may not use this file except in
//  compliance with the License. You may obtain a copy of the License at
//  http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
// __END_LICENSE__


/// \file ToastQuadTreeConfig.h
///
/// Configures a QuadTreeGenerator to generates TOAST tiles for use
/// with e.g. WorldWide Telescope.  This completely replaces the
/// QuadTreeGenerator's Processor engine, so it breaks many of the
/// optional features of QuadTreeGenerator.
///
/// This assumes that its input image is already a global image in
/// TOAST projection, and requires that it be a square image with
/// dimensions 255*2^n+1.  See Cartography/ToastTransform.h.
///
#ifndef __VW_MOSAIC_TOASTQUADTREECONFIG_H__
#define __VW_MOSAIC_TOASTQUADTREECONFIG_H__

#include <boost/filesystem/path.hpp>

#include <vw/Mosaic/QuadTreeGenerator.h>

namespace vw {
namespace mosaic {

  namespace fs = boost::filesystem;

  template <class PixelT>
  class ToastProcessor : public QuadTreeGenerator::ProcessorBase {
    ImageViewRef<PixelT> m_source;

    struct CacheEntry {
      int32 level, x, y;
      ImageView<PixelT> tile;
    };
    typedef std::list<CacheEntry> cache_t;
    cache_t m_cache;

  public:
    template <class ImageT>
    ToastProcessor( QuadTreeGenerator *qtree, ImageT const& source )
      : ProcessorBase( qtree ), m_source( source )
    {}

    void generate( BBox2i const& /*region_bbox*/, const ProgressCallback &progress_callback ) {
      // Note: We ignore the region_bbox parameter!
      for( int32 level = qtree->get_tree_levels()-1; level>=0; --level ) {
        double progress = progress_callback.progress();
        SubProgressCallback spc(progress_callback, progress, 1-(1-progress)/4);
        generate_branch( level, 0, 0, 0, spc );
      }
      progress_callback.report_progress(1);
    }

    // Determine whether or not the given bounding box contains valid
    // data in the source image.  This is just an optimization
    // heuristic: this function is allowed to return false positives,
    // though never false negatives.
    bool sparse_check( BBox2i const& bbox ) const {
      if( ! qtree->sparse_image_check() ) return true;
      QuadTreeGenerator::sparse_image_check_type sparse = qtree->sparse_image_check();

      // Skip a bunch of pathological cases by just never
      // declaring sparsity for very large bounding boxes.
      if( bbox.width() > m_source.cols()/2 || bbox.height() > m_source.rows()/2 ) return true;
      BBox2i center(0,0,m_source.cols(),m_source.rows());
      // Early-out for the most common case: no edge effects
      if( center.contains(bbox) ) {
        return sparse(bbox);
      }
      // Otherwise, we consider all nine possible regions, taking the
      // crazy TOAST edge extension into account.  We return true if
      // any of the nine blocks may overlap valid data.
      bool valid = false;
      // First, the center section
      center.crop(bbox);
      if( ! center.empty() ) {
        valid |= sparse(center);
      }
      // Left (flipped 180)
      BBox2i left( -(m_source.cols()-1), 0, m_source.cols()-1, m_source.rows() );
      left.crop(bbox);
      if( ! left.empty() ) {
        left = BBox2i( 1-left.max().x(), m_source.rows()-left.max().y(), left.width(), left.height() );
        valid |= sparse(left);
      }
      // Right (flipped 180)
      BBox2i right( m_source.cols(), 0, m_source.cols()-1, m_source.rows() );
      right.crop(bbox);
      if( ! right.empty() ) {
        right = BBox2i( 2*m_source.cols()-1-right.max().x(), m_source.rows()-right.max().y(), right.width(), right.height() );
        valid |= sparse(right);
      }
      // Top (flipped 180)
      BBox2i top( 0, -(m_source.rows()-1), m_source.cols(), m_source.rows()-1 );
      top.crop(bbox);
      if( ! top.empty() ) {
        top = BBox2i( m_source.cols()-top.max().x(), 1-top.max().y(), top.width(), top.height() );
        valid |= sparse(top);
      }
      // Bottom (flipped 180)
      BBox2i bottom( 0, m_source.rows(), m_source.cols(), m_source.rows()-1 );
      bottom.crop(bbox);
      if( ! bottom.empty() ) {
        bottom = BBox2i( m_source.cols()-bottom.max().x(), 2*m_source.rows()-1-bottom.max().y(), bottom.width(), bottom.height() );
        valid |= sparse(bottom);
      }
      // Top left
      BBox2i topleft( -(m_source.cols()-1), -(m_source.rows()-1), m_source.cols()-1, m_source.rows()-1 );
      topleft.crop(bbox);
      if( ! topleft.empty() ) {
        topleft += Vector2( m_source.cols()-1, m_source.rows()-1 );
        valid |= sparse(topleft);
      }
      // Top right
      BBox2i topright( m_source.cols(), -(m_source.rows()-1), m_source.cols()-1, m_source.rows()-1 );
      topright.crop(bbox);
      if( ! topright.empty() ) {
        topright += Vector2( -m_source.cols(), m_source.rows()-1 );
        valid |= sparse(topright);
      }
      // Bottom left
      BBox2i botleft( -(m_source.cols()-1), m_source.rows(), m_source.cols()-1, m_source.rows()-1 );
      botleft.crop(bbox);
      if( ! botleft.empty() ) {
        botleft += Vector2( m_source.cols()-1, -m_source.rows() );
        valid |= sparse(botleft);
      }
      // Bottom right
      BBox2i botright( m_source.cols(), m_source.rows(), m_source.cols()-1, m_source.rows()-1 );
      botright.crop(bbox);
      if( ! botright.empty() ) {
        botright += Vector2( -m_source.cols(), -m_source.rows() );
        valid |= sparse(botright);
      }
      return valid;
    }

    // Read a previously-written tile in from disk.  Cache the most
    // recently accessed tiles, since each will be used roughly four
    // times.
    ImageView<PixelT> load_tile( int32 level, int32 x, int32 y ) {
      int32 num_tiles = 1 << level;
      if( x==-1 ) {
        if( y==-1 ) {
          return load_tile(level, num_tiles-1, num_tiles-1);
        }
        if( y==num_tiles ) {
          return load_tile(level, num_tiles-1, 0);
        }
        ImageView<PixelT> tile = load_tile(level, 0, num_tiles-1-y);
        if( tile.is_valid_image() ) return rotate_180(tile);
        else return tile;
      }
      if( x==num_tiles ) {
        if( y==-1 ) {
          return load_tile(level, 0, num_tiles-1);
        }
        if( y==num_tiles ) {
          return load_tile(level, 0, 0);
        }
        ImageView<PixelT> tile = load_tile(level, num_tiles-1, num_tiles-1-y);
        if( tile.is_valid_image() ) return rotate_180(tile);
        else return tile;
      }
      if( y==-1 ) {
        ImageView<PixelT> tile = load_tile(level, num_tiles-1-x, 0);
        if( tile.is_valid_image() ) return rotate_180(tile);
        else return tile;
      }
      if( y==num_tiles ) {
        ImageView<PixelT> tile = load_tile(level, num_tiles-1-x, num_tiles-1);
        if( tile.is_valid_image() ) return rotate_180(tile);
        else return tile;
      }

      // Check the cache
      for( typename cache_t::iterator i=m_cache.begin(); i!=m_cache.end(); ++i ) {
        if( i->level==level && i->x==x && i->y==y ) {
          CacheEntry e = *i;
          m_cache.erase(i);
          m_cache.push_front(e);
          return e.tile;
        }
      }

      // Read it in from disk
      fs::path path( qtree->get_name() );
      std::ostringstream filename;
      filename << level << "/" << x << "/" << y << "." << qtree->get_file_type();
      path /= filename.str();
      ImageView<PixelT> tile;
      if( exists(path) ) {
        read_image( tile, path.string() );
      }

      // Save it in the cache.  The cache size of 1024 tiles was chosen
      // somewhat arbitrarily.
      if( m_cache.size() >= 1024 )
        m_cache.pop_back();
      CacheEntry e;
      e.level = level;
      e.x = x;
      e.y = y;
      e.tile = tile;
      m_cache.push_front(e);

      return tile;
    }

    void generate_branch( int32 branch_level, int32 level, int32 x, int32 y, ProgressCallback const& progress_callback ) {
      progress_callback.report_progress(0);
      progress_callback.abort_if_requested();

      int32 tile_size = qtree->get_tile_size();
      int32 scale = 1<<(qtree->get_tree_levels()-level-1);
      int32 minx = x*scale*(tile_size-1);
      int32 miny = y*scale*(tile_size-1);

      // Compute the region of influence for this branch.
      int32 region_minx = minx - (scale-1);
      int32 region_miny = miny - (scale-1);
      int32 region_size = scale*tile_size + (scale-1);
      BBox2i region_bbox(region_minx, region_miny, region_size, region_size);

      // Early-out for sparse images
      if( ! sparse_check(region_bbox) ) {
        progress_callback.report_progress(1);
        return;
      }

      if( branch_level > level ) {
        generate_branch( branch_level, level+1, 2*x,   2*y,   SubProgressCallback(progress_callback, 0.00, 0.25) );
        generate_branch( branch_level, level+1, 2*x+1, 2*y,   SubProgressCallback(progress_callback, 0.25, 0.50) );
        generate_branch( branch_level, level+1, 2*x,   2*y+1, SubProgressCallback(progress_callback, 0.50, 0.75) );
        generate_branch( branch_level, level+1, 2*x+1, 2*y+1, SubProgressCallback(progress_callback, 0.75, 1.00) );
      }
      else {
        fs::path path( qtree->get_name() );
        std::ostringstream filename;
        filename << level << "/" << x << "/" << y << "." << qtree->get_file_type();
        path /= filename.str();
        ImageView<PixelT> tile;
        if( level == qtree->get_tree_levels()-1 ) {
          // Leaf nodes
          tile = crop(m_source, minx, miny, tile_size, tile_size);
        }
        else {
          // Branch nodes
          ImageView<PixelT> super(4*tile_size-3, 4*tile_size-3);

          // In the WWT implementation of TOAST the pixel centers
          // (rather than the than pixel corners) are grid-aligned, so
          // we need to use an odd-sized antialiasing kernel instead of
          // the usual 2x2 box filter.  The following 5-pixel kernel was
          // optimized to avoid the extra blurring associated with using
          // a kernel wider than 2 pixels.  Math was involved.
          std::vector<float> kernel(5);
          kernel[0] = kernel[4] = -0.0344;
          kernel[1] = kernel[3] = 0.2135;
          kernel[2] = 0.6418;

          for( int j=-1; j<3; ++j ) {
            for( int i=-1; i<3; ++i ) {
              ImageView<PixelT> child = load_tile(level+1,2*x+i,2*y+j);
              if( child.is_valid_image() ) crop(super,(tile_size-1)*(i+1),(tile_size-1)*(j+1),tile_size,tile_size) = child;
            }
          }

          tile = subsample( crop( separable_convolution_filter( super, kernel, kernel, NoEdgeExtension() ),
                                  tile_size-1, tile_size-1, 2*tile_size, 2*tile_size ), 2 );

        }

        if( ! is_transparent(tile) ) {
          create_directories( path.parent_path() );
          write_image( path.string(), tile );
        }

        progress_callback.report_progress(1);
      }

    }
  };

  class ToastQuadTreeConfig {
  public:
    template <class ImageT>
    void configure( QuadTreeGenerator& qtree, ImageT const& image ) const {
      int32 tile_size = qtree.get_tile_size();
      int32 levels = qtree.get_tree_levels();
      VW_ASSERT(image.cols()==image.rows(), ArgumentErr() << "ToastQuadTreeConfig requires a square source image.");
      VW_ASSERT(image.cols()==((tile_size-1)*(1<<(levels-1))+1), ArgumentErr() << "ToastQuadTreeConfig requires a source image with dimensions (tilesize-1)*2^numlevels+1.");
      boost::shared_ptr<QuadTreeGenerator::ProcessorBase> processor( new ToastProcessor<typename ImageT::pixel_type>(&qtree, image) );
      qtree.set_processor( processor );
    }
  };

}} // namespace vw::mosaic

#endif // __VW_MOSAIC_TOASTQUADTREECONFIG_H__
