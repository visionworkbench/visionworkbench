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

/// \file DiskImageView.h
///
/// A read-only disk image view.  This is now just a thin
/// wrapper around the more general ImageResourceView.
///
#ifndef __VW_FILEIO_DISKIMAGEVIEW_H__
#define __VW_FILEIO_DISKIMAGEVIEW_H__

#include <vw/FileIO/DiskImageResource.h>
#include <vw/FileIO/TemporaryFile.h>
#include <vw/Image/ImageResourceView.h>
#include <vw/Image/BlockRasterize.h>
#include <vw/Core/Cache.h>

#include <boost/filesystem/operations.hpp>
#include <string>
#include <map>

namespace vw {

  /// A view of an image on disk.
  template <class PixelT>
  class DiskImageView: public ImageViewBase<DiskImageView<PixelT>> {
    typedef BlockRasterizeView<ImageResourceView<PixelT>> impl_type;
    // This is sort of redundant, but holding both the resource and
    // the block rasterize view simplifies construction and access
    // to the underlying resource.
    boost::shared_ptr<DiskImageResource> m_rsrc;
    impl_type m_impl;

  public:
    typedef typename impl_type::pixel_type     pixel_type;
    typedef typename impl_type::result_type    result_type;
    typedef typename impl_type::pixel_accessor pixel_accessor;

    /// Constructs a DiskImageView of the given file on disk
    /// using the specified cache area. NULL cache means skip it.
    DiskImageView( std::string const& filename, Cache* cache = &vw_system_cache() ):
      m_rsrc( DiskImageResource::open( filename ) ),       // Init file interface
      m_impl( boost::shared_ptr<SrcImageResource>(m_rsrc), // Init memory storage
                m_rsrc->block_read_size(), 1, cache ) {
        // Check for type errors now instead of running into them when we access the image
        try {
          check_convertability(m_impl.child().format(), m_rsrc->format());
        }
        catch(vw::Exception& ex) {
          std::string input_error(ex.what());
          vw_throw( NoImplErr() << "DiskImageView constructor: Image file " << filename
          << " does not match storage format in memory!\n"
          << "    The specific error is:\n        " << input_error
          << "\n    The ImageFormat on disk is  : " << m_rsrc->format()
          << "\n    The ImageFormat in memory is: " << m_impl.child().format() << "\n" );
        }
      }

    /// Constructs a DiskImageView of the given resource using the
    /// specified cache area.
    DiskImageView( boost::shared_ptr<DiskImageResource> resource, Cache* cache = &vw_system_cache())
      : m_rsrc( resource ), m_impl( boost::shared_ptr<SrcImageResource>(m_rsrc), m_rsrc->block_read_size(), 1, cache ) {}

    /// Constructs a DiskImageView of the given resource using the
    /// specified cache area.  Takes ownership of the resource object
    /// (i.e. deletes it when it's done using it).
    DiskImageView( DiskImageResource *resource, Cache* cache = &vw_system_cache() )
      : m_rsrc( resource ), 
        m_impl( boost::shared_ptr<SrcImageResource>(m_rsrc), m_rsrc->block_read_size(), 1, cache ) {}

    /// Constructs a DiskImageView of the given resource using the specified
    /// cache area. Does not take ownership, you must ensure resource stays
    /// valid for the lifetime of DiskImageView
    DiskImageView( DiskImageResource &resource, Cache* cache = &vw_system_cache() )
      : m_rsrc( &resource, NOP() ), 
        m_impl( boost::shared_ptr<SrcImageResource>(m_rsrc), m_rsrc->block_read_size(), 1, cache ) {}

    ~DiskImageView() {}

    int32 cols  () const { return m_impl.cols();   }
    int32 rows  () const { return m_impl.rows();   }
    int32 planes() const { return m_impl.planes(); }

    pixel_accessor origin() const { return m_impl.origin(); }
    result_type operator()( int32 x, int32 y, int32 p = 0 ) const { return m_impl(x,y,p); }

    typedef typename impl_type::prerasterize_type prerasterize_type;
    prerasterize_type prerasterize( BBox2i const& bbox ) const { return m_impl.prerasterize( bbox ); }
    template <class DestT> void rasterize( DestT const& dest, BBox2i const& bbox ) const { m_impl.rasterize( dest, bbox ); }

    std::string filename() const { return m_rsrc->filename(); }

  };

} // namespace vw

#endif // __VW_FILEIO_DISKIMAGEVIEW_H__
