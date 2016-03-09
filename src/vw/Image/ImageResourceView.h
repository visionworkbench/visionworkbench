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


/// \file ImageResourceView.h
///
/// A read-only view of a general image resource.
///
#ifndef __VW_IMAGE_IMAGERESOURCEVIEW_H__
#define __VW_IMAGE_IMAGERESOURCEVIEW_H__

#include <vw/Core/Thread.h>
#include <vw/Image/ImageViewBase.h>
#include <vw/Image/Manipulation.h>
#include <vw/Image/ImageIO.h>

namespace vw {

  /// A view of an image resource.  This class wraps the ImageResource type
  ///  so that image data in a buffer can be treated as a regular ImageView
  ///  style image rather than the buffer style interface provided by ImageResource.
  template <class PixelT>
  class ImageResourceView : public ImageViewBase<ImageResourceView<PixelT> >
  {
  public:
    /// The pixel type of the view.
    typedef PixelT pixel_type;
    typedef PixelT result_type;

    /// The view's pixel accessor type.
    typedef ProceduralPixelAccessor<ImageResourceView<PixelT> > pixel_accessor;

    /// Constructs an ImageResourceView of the given resource.
    ImageResourceView( boost::shared_ptr<SrcImageResource> resource )
      : m_rsrc( resource ), m_planes( m_rsrc->planes() ), m_rsrc_mutex( new Mutex )
    {
      initialize();
    }

    /// Constructs an ImageResourceView of the given resource.  Takes
    /// ownership of the resource object (i.e. deletes it when it's done using it).
    ImageResourceView( SrcImageResource *resource )
      : m_rsrc( resource ), m_planes( m_rsrc->planes() ), m_rsrc_mutex( new Mutex )
    {
      initialize();
    }

    ~ImageResourceView() {}

    // Size queries
    inline int32 cols  () const { return m_rsrc->cols(); }
    inline int32 rows  () const { return m_rsrc->rows(); }
    inline int32 planes() const { return m_planes;       }

    /// Returns the pixel at the given position in the given plane.
    result_type operator()( int32 x, int32 y, int32 plane=0 ) const {
      Mutex::Lock lock(*m_rsrc_mutex);
#if VW_DEBUG_LEVEL > 1
      VW_OUT(VerboseDebugMessage, "image") << "ImageResourceView rasterizing pixel (" << x << "," << y << ")" << std::endl;
#endif
      ImageView<PixelT> buffer(1,1,m_planes);
      read_image( buffer, *m_rsrc, BBox2i(x,y,1,1) );
      return buffer(0,0,plane);
    }

    /// Returns a pixel_accessor pointing to the origin.
    inline pixel_accessor origin() const { return pixel_accessor( *this, 0, 0 ); }

    const SrcImageResource *resource() const { return m_rsrc.get(); }

    typedef CropView<ImageView<pixel_type> > prerasterize_type;
    inline prerasterize_type prerasterize( BBox2i const& bbox ) const {
      ImageView<PixelT> buf( bbox.width(), bbox.height() );
      rasterize( buf, bbox );
      return CropView<ImageView<PixelT> >( buf, BBox2i(-bbox.min().x(),-bbox.min().y(),cols(),rows()) );
    }
    template <class DestT> inline void rasterize( DestT const& dest, BBox2i const& bbox ) const {
      Mutex::Lock lock(*m_rsrc_mutex);
#if VW_DEBUG_LEVEL > 1
      VW_OUT(VerboseDebugMessage, "image") << "ImageResourceView rasterizing bbox " << bbox << std::endl;
#endif
      read_image( dest, *m_rsrc, bbox );
    }

  private:
    void initialize() {
      // If the user has requested a multi-channel pixel type, but the
      // file is a multi-plane, scalar-pixel file, we force a single-plane interpretation.
      if (PixelNumChannels<PixelT>::value > 1 && m_rsrc->pixel_format() == VW_PIXEL_SCALAR) {
        m_planes = 1;
      }

      // On the other hand, the user has requested a scalar pixel type
      // but the file has a multi-channel pixel type, then we force a
      // multi-plane interpretation.
      if (IsScalar<PixelT>::value  && m_rsrc->channels() >= 1 && m_rsrc->planes() == 1) {
        m_planes = m_rsrc->channels();
      }
    }

    boost::shared_ptr<SrcImageResource> m_rsrc;
    int32 m_planes;
    boost::shared_ptr<Mutex> m_rsrc_mutex;
  };

} // namespace vw

#endif // __VW_IMAGE_IMAGERESOURCEVIEW_H__
