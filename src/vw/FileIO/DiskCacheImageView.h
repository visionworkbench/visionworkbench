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

/// \file DiskCacheImageView.h
///
/// An assignable image view that stores the assigned data to a temporary file
/// on disk, and then provides a cached interface (a la DiskImageView) to that
/// data.  The temporary file persists until this object and all copies of this
/// object are destroyed.

#ifndef __VW_FILEIO_DISKCACHEIMAGEVIEW_H__
#define __VW_FILEIO_DISKCACHEIMAGEVIEW_H__

#include <vw/FileIO/DiskImageView.h>
#include <vw/FileIO/DiskImageResource.h>
#include <vw/FileIO/TemporaryFile.h>
#include <vw/Image/ImageResourceView.h>
#include <vw/Image/BlockRasterize.h>
#include <vw/Core/Cache.h>

#include <boost/filesystem/operations.hpp>
#include <string>
#include <map>

namespace vw {

  template <class PixelT>
    class DiskCacheHandle : private boost::noncopyable {
    DiskImageView<PixelT> m_disk_image_view;
    std::string m_filename;

  public:
    template <class ViewT>
    DiskCacheHandle(ImageViewBase<ViewT> const& /*view*/, std::string const& filename) :
      m_disk_image_view(filename), m_filename(filename) {
    }

    ~DiskCacheHandle() {
      VW_OUT(DebugMessage, "fileio") << "DiskCacheImageView: deleting temporary cache file: " << m_filename << "\n";
      boost::filesystem::remove( m_filename );
    }

    inline const DiskImageView<PixelT>& view() const { return m_disk_image_view; }
  };

  template <class PixelT>
  class DiskCacheImageView : public ImageViewBase< DiskCacheImageView<PixelT> > {
  private:
    boost::shared_ptr<DiskCacheHandle<PixelT> > m_handle;
    std::string m_file_type, m_directory;

    template <class ViewT>
    void initialize(ImageViewBase<ViewT> const& view,
                    const ProgressCallback &progress_callback = ProgressCallback::dummy_instance() ) {
      TemporaryFile file(m_directory, false, "vw_cache_", "." + m_file_type);

      VW_OUT(InfoMessage, "fileio") << "Creating disk cache of image in: " << file.filename() << "\n";

      ImageFormat fmt(view.format());
      fmt.pixel_format = PixelFormatID<PixelT>::value;

      boost::scoped_ptr<DiskImageResource> r(DiskImageResource::create( file.filename(), fmt));
      if (r->has_block_write())
        r->set_block_write_size(Vector2i(vw_settings().default_tile_size(), vw_settings().default_tile_size()));
      block_write_image(*r, pixel_cast_rescale<PixelT>(view), progress_callback);
      r->flush();

      m_handle = boost::shared_ptr<DiskCacheHandle<PixelT> >(new DiskCacheHandle<PixelT>(view.impl(), file.filename()));
    }

  public:
    typedef typename DiskImageView<PixelT>::pixel_type     pixel_type;
    typedef typename DiskImageView<PixelT>::result_type    result_type;
    typedef typename DiskImageView<PixelT>::pixel_accessor pixel_accessor;

    /// Create a temporary image view cache file on disk using a
    /// system supplied temporary filename.
    template <class ViewT>
    DiskCacheImageView(ImageViewBase<ViewT> const& view, std::string const& file_type = "tif",
                       const ProgressCallback &progress_callback = ProgressCallback::dummy_instance(),
                       std::string const& directory = "/tmp" ) :
      m_file_type(file_type), m_directory(directory) {
      this->initialize(view.impl(), progress_callback);
    }

    inline int32 cols  () const { return m_handle->view().cols();   }
    inline int32 rows  () const { return m_handle->view().rows();   }
    inline int32 planes() const { return m_handle->view().planes(); }

    inline pixel_accessor origin() const { return m_handle->view().origin(); }

    inline result_type operator()( int32 i, int32 j, int32 p=0 ) const { return m_handle->view()(i, j, p); }

    template <class ViewT>
    DiskCacheImageView& operator=( ImageViewBase<ViewT> const& view ) {
      this->initialize(view.impl());
      return *this;
    }

    /// \cond INTERNAL
    typedef typename DiskImageView<PixelT>::prerasterize_type prerasterize_type;
    inline prerasterize_type prerasterize( BBox2i const& bbox ) const {
      return m_handle->view().prerasterize(bbox);
    }
    template <class DestT> inline void rasterize( DestT const& dest, BBox2i const& bbox ) const {
      vw::rasterize( prerasterize(bbox), dest, bbox );
    }
    /// \endcond
  };
  
} // namespace vw

#endif // __VW_FILEIO_DISKCACHEIMAGEVIEW_H__
