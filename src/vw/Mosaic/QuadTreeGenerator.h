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

#include <boost/function.hpp>
#include <vw/Math/BBox.h>
#include <vw/Image/ImageView.h>
#include <vw/Image/ImageViewRef.h>
#include <vw/Image/EdgeExtension.h>
#include <vw/Image/Manipulation.h>
#include <vw/Image/Algorithms.h>
#include <vw/Image/Filter.h>
#include <vw/Image/ImageIO.h>
#include <vw/Core/ProgressCallback.h>
#include <vw/Core/Stopwatch.h>
#include <vw/Mosaic/SparseTileCheck.h>

namespace vw {
namespace mosaic {

  // This feature, or something like it, should be implemented better and 
  // moved somewhere into the Image module.
  template <class PixelT>
  ImageView<PixelT> box_subsample( ImageView<PixelT> const& image, Vector2i const& scale ) {
    std::vector<double> xkernel(scale.x()), ykernel(scale.y());
    for( int x=0; x<scale.x(); ++x ) xkernel[x] = 1.0 / scale.x();
    for( int y=0; y<scale.y(); ++y ) ykernel[y] = 1.0 / scale.y();
    ImageView<PixelT> result( image.cols()/scale.x(), image.rows()/scale.y() );
    rasterize( subsample( separable_convolution_filter( image, xkernel, ykernel, scale.x()-1, scale.y()-1 ), scale.x(), scale.y() ), result );
    return result;
  }

  class QuadTreeGenerator {
  public:
    struct TileInfo {
      std::string name, filepath, filetype;
      BBox2i image_bbox, region_bbox;
    };

    typedef boost::function<std::string(QuadTreeGenerator const&, std::string const&)> image_path_func_type;
    typedef boost::function<std::vector<std::pair<std::string,BBox2i> >(QuadTreeGenerator const&, std::string const&, BBox2i const&)> branch_func_type;
    typedef boost::function<boost::shared_ptr<ImageResource>(QuadTreeGenerator const&, TileInfo const&, ImageFormat const&)> tile_resource_func_type;
    typedef boost::function<void(QuadTreeGenerator const&, TileInfo const&)> metadata_func_type;
    typedef boost::function<bool(BBox2i const&)> sparse_tile_check_type;

    template <class ImageT>
    QuadTreeGenerator( ImageViewBase<ImageT> const& image, std::string const& tree_name = "output.qtree" )
      : m_tree_name( tree_name ),
        m_tile_size( 256 ),
        m_file_type( "png" ),
        m_crop_bbox(),
        m_crop_images( false ),
        m_dimensions( image.impl().cols(), image.impl().rows() ),
        m_processor( new Processor<typename ImageT::pixel_type>( this, image.impl() ) ),
        m_image_path_func( &simple_image_path ),
        m_branch_func( &default_branch_func ),
        m_tile_resource_func( &default_tile_resource_func ),
        m_metadata_func(),
        m_sparse_tile_check( SparseTileCheck<ImageT>(image.impl()) )
    {}

    virtual ~QuadTreeGenerator() {}

    void generate( const ProgressCallback &progress_callback = ProgressCallback::dummy_instance() );
    
    std::string const& get_name() const {
      return m_tree_name;
    }

    void set_name( std::string const& name ) {
      m_tree_name = name;
    }

    BBox2i const& get_crop_bbox() const {
      return m_crop_bbox;
    }

    void set_crop_bbox( BBox2i const& bbox ) {
      VW_ASSERT( BBox2i(Vector2i(), m_dimensions).contains(bbox),
		 ArgumentErr() << "Requested QuadTree bounding box exceeds source dimensions!" );
      m_crop_bbox = bbox;
    }

    std::string const& get_file_type() const {
      return m_file_type;
    }

    void set_file_type( std::string const& extension ) {
      m_file_type = extension;
    }

    int32 get_tile_size() const {
      return m_tile_size;
    }

    void set_tile_size( int32 size ) {
      m_tile_size = size;
    }

    int32 get_tree_levels() const {
      int32 maxdim = (std::max)( m_dimensions.x(), m_dimensions.y() );
      int32 tree_levels = 1 + int32( ceil( log( maxdim/(double)(m_tile_size) ) / log(2.0) ) );
      if (tree_levels < 1) tree_levels = 1;
      return tree_levels;
    }

    Vector2i const& get_dimensions() const {
      return m_dimensions;
    }

    bool get_crop_images() const {
      return m_crop_images;
    }

    void set_crop_images( bool crop ) {
      m_crop_images = crop;
    }

    void set_image_path_func( image_path_func_type image_path_func ) {
      m_image_path_func = image_path_func;
    }

    std::string image_path( std::string const& name ) const {
      return m_image_path_func( *this, name );
    }

    void set_branch_func( branch_func_type const& branch_func ) {
      m_branch_func = branch_func;
    }

    std::vector<std::pair<std::string,BBox2i> > branches( std::string const& name, BBox2i const& region ) const {
      return m_branch_func( *this, name, region );
    }

    void set_tile_resource_func( tile_resource_func_type const& tile_resource_func ) {
      m_tile_resource_func = tile_resource_func;
    }

    void set_metadata_func( metadata_func_type metadata_func ) {
      m_metadata_func = metadata_func;
    }

    // Makes paths of the form "path/name/r0132.jpg"
    static std::string simple_image_path( QuadTreeGenerator const& qtree, std::string const& name );

    // Makes paths of the form "path/name/r01/r0132.jpg"
    static std::string tiered_image_path( QuadTreeGenerator const& qtree, std::string const& name, int32 levels_per_directory = 3 );

    // Makes paths of the form "path/name/013/0132.jpg" (and "path/name/name.jpg" for the top level)
    static std::string named_tiered_image_path( QuadTreeGenerator const& qtree, std::string const& name, int32 levels_per_directory = 3 );

    // The default quad-tree branching function
    static std::vector<std::pair<std::string,BBox2i> > default_branch_func(QuadTreeGenerator const& qtree, std::string const& name, BBox2i const& region);

    // The default resource function, creates standard disk image resources
    static boost::shared_ptr<ImageResource> default_tile_resource_func( QuadTreeGenerator const& qtree, TileInfo const& info, ImageFormat const& format );

  protected:
    class ProcessorBase {
    protected:
      QuadTreeGenerator *qtree;
    public:
      ProcessorBase( QuadTreeGenerator *qtree ) : qtree(qtree) {}
      virtual ~ProcessorBase() {}
      virtual void generate( BBox2i const& bbox, const ProgressCallback &progress_callback ) = 0;
    };

    template <class PixelT>
    class Processor : public ProcessorBase {
      ImageViewRef<PixelT> m_source;
      
    public:
      template <class ImageT>
      Processor( QuadTreeGenerator *qtree, ImageT const& source )
	: ProcessorBase( qtree ), m_source( source )
      {}
      
      void generate( BBox2i const& region_bbox, const ProgressCallback &progress_callback ) {
	generate_branch( "", region_bbox, progress_callback );
      }

      ImageView<PixelT> generate_branch( std::string const& name, BBox2i const& region_bbox, const ProgressCallback &progress_callback ) {
	progress_callback.report_progress(0);
	progress_callback.abort_if_requested();
  
	ImageView<PixelT> image;
	TileInfo info;
	info.name = name;
	info.region_bbox = region_bbox;

	BBox2i crop_bbox(Vector2i(), qtree->m_dimensions);
	if( ! qtree->get_crop_bbox().empty() ) crop_bbox.crop( qtree->get_crop_bbox() );
	info.image_bbox = info.region_bbox;
	info.image_bbox.crop( crop_bbox );
	if( info.image_bbox.empty() ) return image;
  
	if( qtree->m_sparse_tile_check && ! qtree->m_sparse_tile_check(info.region_bbox) ) return image;

	Vector2i scale = info.region_bbox.size() / qtree->m_tile_size;
 
	std::vector<std::pair<std::string, BBox2i> > children = qtree->m_branch_func(*qtree,info.name,info.region_bbox);
	if( children.empty() ) {
	  image = crop( m_source, info.image_bbox );
	  if( info.image_bbox != info.region_bbox ) {
	    image = edge_extend( image, info.region_bbox - info.image_bbox.min(), ZeroEdgeExtension() );
	  }
	  if( info.region_bbox.width() != qtree->m_tile_size || info.region_bbox.height() != qtree->m_tile_size ) {
	    image = subsample( image, scale.x(), scale.y() );
	  }
	}
	else {
	  image.set_size(qtree->m_tile_size,qtree->m_tile_size);
	  for( unsigned i=0; i<children.size(); ++i ) {
	    BBox2i dst_bbox = elem_quot( children[i].second - info.region_bbox.min(), scale );
	    SubProgressCallback spc(progress_callback, i/(double)children.size(), (i+1)/(double)children.size());
	    ImageView<PixelT> child = generate_branch(children[i].first, children[i].second, spc);
	    if( ! child ) continue;
	    crop(image,dst_bbox) = box_subsample( child, elem_quot(qtree->m_tile_size,dst_bbox.size()) );
	  }
	}
  
	ImageView<PixelT> cropped_image = image;
	if( qtree->m_crop_images ) {
	  BBox2i data_bbox = elem_quot( info.image_bbox-info.region_bbox.min(), scale );
	  if( PixelHasAlpha<PixelT>::value )
	    data_bbox.crop( nonzero_data_bounding_box( image ) );
	  if( data_bbox.width() != qtree->m_tile_size || data_bbox.height() != qtree->m_tile_size ) {
	    if( data_bbox.empty() ) cropped_image.reset();
	    else cropped_image = crop( image, data_bbox );
	    info.image_bbox = elem_prod(data_bbox,scale) + info.region_bbox.min();
	  }
	}
  
	if( qtree->m_file_type == "auto" ) {
	  if( is_opaque( cropped_image ) ) info.filetype += ".jpg";
	  else info.filetype += ".png";
	}
	else {
	  info.filetype = "." + qtree->m_file_type;
	}
  
	info.filepath = qtree->m_image_path_func( *qtree, info.name );
	if( cropped_image ) {
	  ScopedWatch sw("QuadTreeGenerator::write_tile");
	  boost::shared_ptr<ImageResource> r = qtree->m_tile_resource_func( *qtree, info, cropped_image.format() );
	  write_image( *r, cropped_image );
	}
	if( qtree->m_metadata_func ) qtree->m_metadata_func( *qtree, info );
  
	progress_callback.report_progress(1);
	return image;
      }
    };

    template<class> friend class Processor;

    std::string m_tree_name;
    int32 m_tile_size;
    std::string m_file_type;
    BBox2i m_crop_bbox;
    bool m_crop_images;
    Vector2i m_dimensions;
    boost::shared_ptr<ProcessorBase> m_processor;

    image_path_func_type m_image_path_func;
    branch_func_type m_branch_func;
    tile_resource_func_type m_tile_resource_func;
    metadata_func_type m_metadata_func;
    sparse_tile_check_type m_sparse_tile_check;
  };

} // namespace mosaic
} // namespace vw

#endif // __VW_MOSAIC_QUADTREEGENERATOR_H__
