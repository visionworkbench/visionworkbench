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

// Turn off warnings about things we can't control
#define BOOST_ALLOW_DEPRECATED_HEADERS
#define BOOST_BIND_GLOBAL_PLACEHOLDERS  
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#pragma GCC diagnostic ignored "-Wunknown-pragmas"
#include <boost/function.hpp>
#pragma GCC diagnostic pop

#include <vw/Core/ProgressCallback.h>
#include <vw/Core/Stopwatch.h>
#include <vw/Image/EdgeExtension.h>
#include <vw/Image/SparseImageCheck.h>
#include <vw/Image/ImageView.h>
#include <vw/Image/ImageViewRef.h>
#include <vw/Image/Algorithms.h>
#include <vw/Image/Filter.h>
#include <vw/Image/Manipulation.h>


namespace vw {
namespace mosaic {

  // This feature, or something like it, should be implemented better and
  // moved somewhere into the Image module.
  /// ???
  template <class PixelT>
  ImageView<PixelT> box_subsample( ImageView<PixelT> const& image, Vector2i const& scale ) {
    std::vector<double> xkernel(scale.x()), ykernel(scale.y());
    for( int x=0; x<scale.x(); ++x ) 
      xkernel[x] = 1.0 / scale.x();
    for( int y=0; y<scale.y(); ++y ) 
      ykernel[y] = 1.0 / scale.y();
    
    ImageView<PixelT> result( image.cols()/scale.x(), image.rows()/scale.y() );
    rasterize( subsample( separable_convolution_filter( image, xkernel, ykernel, scale.x()-1, scale.y()-1 ), scale.x(), scale.y() ), result );
    return result;
  }
  
  /// A class that generates filesystem-based quadtrees of large images.
  class QuadTreeGenerator {
  public:
    struct TileInfo {
      std::string name, filepath, filetype;
      BBox2i image_bbox, region_bbox;
    };

    typedef boost::function<std::string(QuadTreeGenerator const&, std::string const&)> 
        image_path_func_type;
    typedef boost::function<std::vector<std::pair<std::string,BBox2i> >(QuadTreeGenerator const&, std::string const&, BBox2i const&)> 
        branch_func_type;
    typedef boost::function<boost::shared_ptr<DstImageResource>(QuadTreeGenerator const&, TileInfo const&, ImageFormat const&)> 
        tile_resource_func_type;
    typedef boost::function<void(QuadTreeGenerator const&, TileInfo const&)> 
        metadata_func_type;
    typedef boost::function<bool(BBox2i const&)> 
        sparse_image_check_type;

    class ProcessorBase {
    protected:
      QuadTreeGenerator *qtree;
    public:
      ProcessorBase( QuadTreeGenerator *qtree ) : qtree(qtree) {}
      virtual ~ProcessorBase() {}
      virtual void generate( BBox2i const& bbox, const ProgressCallback &progress_callback ) = 0;
    };

    template <class ImageT>
    QuadTreeGenerator( ImageViewBase<ImageT> const& image, std::string const& tree_name = "output.qtree" )
      : m_tree_name( tree_name ),
        m_tile_size( 256 ),
        m_file_type( "png" ),
        m_crop_bbox(),
        m_crop_images( false ),
        m_cull_images( false ),
        m_dimensions( image.impl().cols(), image.impl().rows() ),
        m_processor( new Processor<typename ImageT::pixel_type>( this, image.impl() ) ),
        m_image_path_func( simple_image_path() ),
        m_branch_func( default_branch_func() ),
        m_tile_resource_func( default_tile_resource_func() ),
        m_metadata_func(),
        m_sparse_image_check( SparseImageCheck<ImageT>(image.impl()) )
    {}

    virtual ~QuadTreeGenerator() {}

    void set_processor( boost::shared_ptr<ProcessorBase> const& processor ) {
      m_processor = processor;
    }

    void generate( const ProgressCallback &progress_callback = ProgressCallback::dummy_instance() );

    void set_crop_bbox( BBox2i const& bbox ) {
      VW_ASSERT( BBox2i(Vector2i(), m_dimensions).contains(bbox),
                 ArgumentErr() << "Requested QuadTree bounding box exceeds source dimensions!" );
      m_crop_bbox = bbox;
    }

    // Compute number of tree levels required with a downsample factor of 2
    int32 get_tree_levels() const {
      int32 maxdim      = (std::max)( m_dimensions.x(), m_dimensions.y() ); // Get largest dimension
      int32 tree_levels = 1 + int32( ceil( log( maxdim/(double)(m_tile_size) ) / log(2.0) ) );
      if (tree_levels < 1) 
        tree_levels = 1; // Can't have less than one level!
      return tree_levels;
    }

    // Simple "get" functions
    std::string const& get_name()        const { return m_tree_name;   }
    BBox2i      const& get_crop_bbox()   const { return m_crop_bbox;   }
    std::string const& get_file_type()   const { return m_file_type;   }
    int32              get_tile_size()   const { return m_tile_size;   }
    Vector2i    const& get_dimensions()  const { return m_dimensions;  }
    bool               get_crop_images() const { return m_crop_images; }
    bool               get_cull_images() const { return m_cull_images; }
    sparse_image_check_type const& sparse_image_check() const { return m_sparse_image_check; }


    // Simple "set" functions
    void set_name              (std::string             const& name              ) {m_tree_name          = name;              }
    void set_file_type         (std::string             const& extension         ) {m_file_type          = extension;         }
    void set_tile_size         (int32                          size              ) {m_tile_size          = size;              }
    void set_crop_images       (bool                           crop              ) {m_crop_images        = crop;              }
    void set_cull_images       (bool                           cull              ) {m_cull_images        = cull;              }
    void set_image_path_func   (image_path_func_type           image_path_func   ) {m_image_path_func    = image_path_func;   }
    void set_branch_func       (branch_func_type        const& branch_func       ) {m_branch_func        = branch_func;       }
    void set_tile_resource_func(tile_resource_func_type const& tile_resource_func) {m_tile_resource_func = tile_resource_func;}
    void set_metadata_func     (metadata_func_type             metadata_func     ) {m_metadata_func      = metadata_func;     }
    void set_sparse_image_check(sparse_image_check_type const& func              ) {m_sparse_image_check = func;              }


    std::string image_path( std::string const& name ) const {
      return m_image_path_func( *this, name );
    }

    std::vector<std::pair<std::string,BBox2i> > branches( std::string const& name, BBox2i const& region ) const {
      return m_branch_func( *this, name, region );
    }

    boost::shared_ptr<DstImageResource> tile_resource(TileInfo const& info, ImageFormat const& format) const {
      return m_tile_resource_func( *this, info, format );
    }

    void make_tile_metadata( TileInfo const& info ) const {
      if( m_metadata_func ) {
        m_metadata_func( *this, info );
      }
    }

    // Here we have a set of functor classes which contain different implementations of a path generation function

    /// Makes paths of the form "path/name/r0132.jpg".  More accurately: [qtree.get_name()]/r[name]
    struct simple_image_path {
      std::string operator()( QuadTreeGenerator const& qtree, std::string const& name );
    };

    /// Makes paths of the form "path/name/r01/r0132.jpg".
    struct tiered_image_path {
      std::string operator()( QuadTreeGenerator const& qtree, std::string const& name, int32 levels_per_directory = 3 );
    };

    /// Makes paths of the form "path/name/013/0132.jpg" (and "path/name/name.jpg" for the top level)
    struct named_tiered_image_path {
      std::string operator()( QuadTreeGenerator const& qtree, std::string const& name, int32 levels_per_directory = 3 );
    };


    //   
  
    /// The default quad-tree branching function
    struct default_branch_func {
      std::vector<std::pair<std::string,BBox2i> > operator()(QuadTreeGenerator const& qtree, std::string const& name, BBox2i const& region);
    };

    /// The default resource function, creates standard disk image resources
    struct default_tile_resource_func {
      boost::shared_ptr<DstImageResource> operator()( QuadTreeGenerator const& qtree, TileInfo const& info, ImageFormat const& format );
    };

  protected:
  
    /// Secret class that contains all the high level tree generation logic
    template <class PixelT>
    class Processor : public ProcessorBase {
      ImageViewRef<PixelT> m_source;

    public:
      /// Construct the image with the qtree object and the full resolution source image
      template <class ImageT>
      Processor( QuadTreeGenerator *qtree, ImageT const& source )
        : ProcessorBase( qtree ), m_source( source )
      {}

      /// Top level call to generate a qtree from a specified region of the input image.
      void generate( BBox2i const& region_bbox, const ProgressCallback &progress_callback ) {
        // Just redirect to the branch function leaving the name blank.
        generate_branch( "", region_bbox, progress_callback );
      }

      /// Generate all images and metadata files (all the way down the tree) for a named region of the input image.
      /// - Note that region_bbox is always in the original source image, not the parent of this particular branch.
      ImageView<PixelT> generate_branch( std::string const& name, BBox2i const& region_bbox, const ProgressCallback &progress_callback ) {
        progress_callback.report_progress(0);
        progress_callback.abort_if_requested();

        ImageView<PixelT> image;
        TileInfo info;
        info.name = name;
        info.region_bbox = region_bbox;

        BBox2i crop_bbox(Vector2i(), qtree->get_dimensions());
        if( ! qtree->get_crop_bbox().empty() ) 
          crop_bbox.crop( qtree->get_crop_bbox() );
        info.image_bbox = info.region_bbox;
        info.image_bbox.crop( crop_bbox );

        if( info.image_bbox.empty() ) {
          if( ! (qtree->get_crop_images() || qtree->get_cull_images()) )
            image.set_size( qtree->get_tile_size(), qtree->get_tile_size() );
          return image;
        }

        if( qtree->m_sparse_image_check && ! qtree->m_sparse_image_check(info.region_bbox) ) 
            return image;

        Vector2i scale = info.region_bbox.size() / qtree->m_tile_size;

        // Call function to compute which children belong to this tile.
        // - Each child contains a name and a bounding box.
        std::vector<std::pair<std::string, BBox2i> > children = qtree->m_branch_func(*qtree, info.name, info.region_bbox);
        
        if( children.empty() ) { // This is the highest resolution level of tiles (bottom of tree)
          image = crop( m_source, info.image_bbox ); // Extract portion of source image
          if( info.image_bbox != info.region_bbox ) { // Pad with zero pixels if needed
            image = edge_extend( image, info.region_bbox - info.image_bbox.min(), ZeroEdgeExtension() );
          }
          if( (info.region_bbox.width() != qtree->m_tile_size) || (info.region_bbox.height() != qtree->m_tile_size) ) {
            image = subsample( image, scale.x(), scale.y() ); // Resample image to the output tile size
          }
        }
        else { // One or more sub-levels below this image, generate and copy from them one at a time
          image.set_size(qtree->m_tile_size,qtree->m_tile_size); // Initialize empty image
          double total_area = (double) info.image_bbox.width() * info.image_bbox.height();
          // For each child
          for( unsigned i=0; i<children.size(); ++i ) { 
            // Double check that the BBox for the child is fully contained in the parent image
            BBox2i image_bbox = children[i].second;
            image_bbox.crop( info.image_bbox ); 
            if( image_bbox.empty() ) 
              continue; // Skip the child if no overlap
            
            // Recursively call this function on this child and get its image
            // - First set up a progress report callback function
            double child_area = (double) image_bbox.width() * image_bbox.height();
            double progress   = progress_callback.progress();
            SubProgressCallback spc( progress_callback, progress, progress + child_area/total_area );
            ImageView<PixelT> child = generate_branch(children[i].first, children[i].second, spc); // (name, BBox, callback)
            if( ! child.is_valid_image() ) 
              continue;
            
            BBox2i dst_bbox = elem_quot( children[i].second - info.region_bbox.min(), scale );            // Compute this child's ROI in the current tile.
            crop(image,dst_bbox) = box_subsample( child, elem_quot(qtree->m_tile_size,dst_bbox.size()) ); // Copy and resample the child image to the destination ROI
          }
        }

        ImageView<PixelT> cropped_image = image;
        if( qtree->m_crop_images || qtree->m_cull_images ) {
        
          BBox2i data_bbox = elem_quot( info.image_bbox-info.region_bbox.min(), scale );
          if( PixelHasAlpha<PixelT>::value )
            data_bbox.crop( nonzero_data_bounding_box( image ) );
            
          if( data_bbox.width() != qtree->m_tile_size || data_bbox.height() != qtree->m_tile_size ) {
            if( data_bbox.empty() ) { 
              cropped_image.reset();
            }
            else { 
              if( qtree->m_crop_images ) {
                cropped_image = crop( image, data_bbox );
              }
            }
            info.image_bbox = elem_prod(data_bbox,scale) + info.region_bbox.min();
          }
        } // End crop or cull images case

        if( qtree->m_file_type == "auto" ) {
          if( is_opaque( cropped_image ) ) 
            info.filetype += ".jpg";  // Use jpg for images with no alpha channel
          else 
            info.filetype += ".png"; // Use png for images with transparency
        }
        else { // User must have submitted the output file type
          info.filetype = "." + qtree->m_file_type;
        }

        // Retrieve the output path for this tile and write it to disk
        info.filepath = qtree->m_image_path_func( *qtree, info.name );
        if( cropped_image.is_valid_image() ) {
          ScopedWatch sw("QuadTreeGenerator::write_tile");
          boost::shared_ptr<DstImageResource> r = qtree->m_tile_resource_func( *qtree, info, cropped_image.format() );
          write_image( *r, cropped_image );
        }
        // Call function to take care of any extra tile metadata tasks
        if( qtree->m_metadata_func ) 
          qtree->m_metadata_func( *qtree, info );

        progress_callback.report_progress(1);
        return image;
      }
    }; // End class Processor

    template<class> friend class Processor;

    // A bunch of public variables
    std::string m_tree_name;
    int32       m_tile_size;
    std::string m_file_type;
    BBox2i      m_crop_bbox;
    bool        m_crop_images;
    bool        m_cull_images;
    Vector2i    m_dimensions;
    boost::shared_ptr<ProcessorBase> m_processor;

    // A bunch of function types
    image_path_func_type    m_image_path_func;
    branch_func_type        m_branch_func;
    tile_resource_func_type m_tile_resource_func;
    metadata_func_type      m_metadata_func;
    sparse_image_check_type m_sparse_image_check;
  };

} // namespace mosaic
} // namespace vw

#endif // __VW_MOSAIC_QUADTREEGENERATOR_H__
