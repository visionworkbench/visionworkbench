// __BEGIN_LICENSE__
// Copyright (C) 2006-2009 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifndef __VWV_TEXTURE_CACHE_H__
#define __VWV_TEXTURE_CACHE_H__

#ifdef __linux__
// This is required to get prototypes, according to the opengl linux abi
#define GL_GLEXT_PROTOTYPES 1
#endif

#ifdef __APPLE__
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#else // Linux
#include <GL/gl.h>
#include <GL/glu.h>
#endif 

#include <vw/Core/Cache.h>
#include <vw/Core/Thread.h>
#include <vw/Core/Log.h>
#include <vw/Image/ImageView.h>
#include <vw/Image/ViewImageResource.h>
#include <vw/Image/Algorithms.h>
#include <vw/Image/ImageResourceView.h>
#include <vw/Image/ImageViewRef.h>
#include <vw/Math/BBox.h>
#include <vw/Math/Vector.h>
#include <vw/FileIO.h>

// Forward Declarations
class CachedTextureRenderer;
class TextureFetchTask;
class TextureRequest;

struct TextureRecordBase {
  vw::BBox2i bbox; 
  int lod;
  GLuint texture_id;
  CachedTextureRenderer* parent;
  virtual ~TextureRecordBase() {}
};

// --------------------------------------------------------------
//                     CachedTextureRenderer
// --------------------------------------------------------------

class CachedTextureRenderer {
  std::list<boost::shared_ptr<TextureRequest> > m_incoming_requests;
  vw::Mutex m_incoming_request_mutex;

protected:
  bool m_needs_redraw;

public:

  // These are defined in the subclass: vwv_GlPreviewWidget
  virtual GLuint allocate_texture(vw::ViewImageResource const block) = 0;
  virtual void deallocate_texture(GLuint texture_id) = 0;
  
  virtual ~CachedTextureRenderer() { m_incoming_requests.clear(); }

  virtual void request_allocation(boost::shared_ptr<TextureRecordBase> texture_record, 
                                  vw::ViewImageResource const block);

  virtual void request_deallocation(boost::shared_ptr<TextureRecordBase> texture_record);

  virtual void process_allocation_requests();
};

// --------------------------------------------------------------
//                     GlTextureHandle
// --------------------------------------------------------------

struct GlTextureHandleBase {
  virtual GLuint get_texture_id() const = 0;
  virtual ~GlTextureHandleBase() {}
};

template <class PixelT>
class GlTextureHandle : public GlTextureHandleBase {
  boost::shared_ptr<TextureRecordBase> m_record;
  
public:
  GlTextureHandle(vw::ImageResourceView<PixelT> image, 
                  boost::shared_ptr<TextureRecordBase> record) : m_record(record) {

    // Rasterize the requested block of memory.
    vw::ImageView<PixelT> cropped = crop( image, record->bbox );

    // For debugging:
    //
    // std::ostringstream ostr;
    // ostr << "debug_" << record->lod << "_" << record->bbox.min().x() << "_" 
    //      << record->bbox.min().y() << ".tif";
    // vw::write_image(ostr.str(), cropped);

    vw::ImageView<PixelT> block = subsample( cropped, pow(2,record->lod) );
    m_record->parent->request_allocation( record, vw::ViewImageResource(block) );
    vw::vw_out(vw::VerboseDebugMessage) << "GlTextureHandle requesting allocation (" 
                                        << m_record->texture_id << ") -- " 
                                        << m_record->bbox << " @ " << m_record->lod << "\n";
  }

  virtual GLuint get_texture_id() const { return m_record->texture_id; }

  virtual ~GlTextureHandle() {
    vw::vw_out(vw::VerboseDebugMessage) << "-> GlTextureHandle requesting decallocation (" << m_record->texture_id << ")\n";
    m_record->parent->request_deallocation(m_record);
  }
};

// --------------------------------------------------------------
//                     GlTextureGenerator
// --------------------------------------------------------------

class GlTextureGenerator {
  boost::shared_ptr<vw::ImageResource> m_rsrc;
  boost::shared_ptr<TextureRecordBase> m_record;

public:
  typedef GlTextureHandleBase value_type;
  
  GlTextureGenerator( boost::shared_ptr<vw::ImageResource> const& rsrc, 
                      boost::shared_ptr<TextureRecordBase> record) :
    m_rsrc( rsrc ), m_record(record) {}

  size_t size() const {
    return m_record->bbox.width()/pow(2,m_record->lod) * 
      m_record->bbox.height()/pow(2,m_record->lod) * 
      m_rsrc->planes() * m_rsrc->channels() 
      * 2;  // <-- The textures are stored as 16-bit half's on the
            // GPU, so they are 2 bytes per channel
  }
  
  boost::shared_ptr<GlTextureHandleBase> generate() const {
    // This is yet another block of ugly run-time dispatch code.  Some
    // day we will hopefully unify the whole ImageView/ImageResource
    // divide.  In the meantime, here we go...
    vw::ChannelTypeEnum channel_type = m_rsrc->channel_type();
    vw::PixelFormatEnum pixel_format = m_rsrc->pixel_format();

    switch(pixel_format) {
    case vw::VW_PIXEL_GRAY:
      switch(channel_type) {
      case vw::VW_CHANNEL_UINT8:  
        return boost::shared_ptr<GlTextureHandleBase> ( new GlTextureHandle<vw::PixelGray<vw::uint8> >(m_rsrc, m_record) );
        break;
      case vw::VW_CHANNEL_UINT16:
        return boost::shared_ptr<GlTextureHandleBase> ( new GlTextureHandle<vw::PixelGray<vw::uint16> >(m_rsrc, m_record) );
        break;
      case vw::VW_CHANNEL_INT16: 
        return boost::shared_ptr<GlTextureHandleBase> ( new GlTextureHandle<vw::PixelGray<vw::int16> >(m_rsrc, m_record) );
        break;
      case vw::VW_CHANNEL_FLOAT32: 
        return boost::shared_ptr<GlTextureHandleBase> ( new GlTextureHandle<vw::PixelGray<vw::float32> >(m_rsrc, m_record) );
        break;
      default:
        vw::vw_throw(vw::IOErr() << "GlTextureHandle: generate() failed. Unknown channel type: " << channel_type << ".\n");
      }
    case vw::VW_PIXEL_GRAYA:
      switch(channel_type) {
      case vw::VW_CHANNEL_UINT8:  
        return boost::shared_ptr<GlTextureHandleBase> ( new GlTextureHandle<vw::PixelGrayA<vw::uint8> >(m_rsrc, m_record) );
        break;
      case vw::VW_CHANNEL_UINT16:
        return boost::shared_ptr<GlTextureHandleBase> ( new GlTextureHandle<vw::PixelGrayA<vw::uint16> >(m_rsrc, m_record) );
        break;
      case vw::VW_CHANNEL_INT16: 
        return boost::shared_ptr<GlTextureHandleBase> ( new GlTextureHandle<vw::PixelGrayA<vw::int16> >(m_rsrc, m_record) );
        break;
      case vw::VW_CHANNEL_FLOAT32: 
        return boost::shared_ptr<GlTextureHandleBase> ( new GlTextureHandle<vw::PixelGrayA<vw::float32> >(m_rsrc, m_record) );
        break;
      default:
        vw::vw_throw(vw::IOErr() << "GlTextureHandle: generate() failed. Unknown channel type: " << channel_type << ".\n");
      }
    case vw::VW_PIXEL_RGB:
      switch(channel_type) {
      case vw::VW_CHANNEL_UINT8:  
        return boost::shared_ptr<GlTextureHandleBase> ( new GlTextureHandle<vw::PixelRGB<vw::uint8> >(m_rsrc, m_record) );
        break;
      case vw::VW_CHANNEL_UINT16:
        return boost::shared_ptr<GlTextureHandleBase> ( new GlTextureHandle<vw::PixelRGB<vw::uint16> >(m_rsrc, m_record) );
        break;
      case vw::VW_CHANNEL_INT16: 
        return boost::shared_ptr<GlTextureHandleBase> ( new GlTextureHandle<vw::PixelRGB<vw::int16> >(m_rsrc, m_record) );
        break;
      case vw::VW_CHANNEL_FLOAT32: 
        return boost::shared_ptr<GlTextureHandleBase> ( new GlTextureHandle<vw::PixelRGB<vw::float32> >(m_rsrc, m_record) );
        break;
      default:
        vw::vw_throw(vw::IOErr() << "GlTextureHandle: generate() failed. Unknown channel type: " << channel_type << ".\n");
      }
    case vw::VW_PIXEL_RGBA:
      switch(channel_type) {
      case vw::VW_CHANNEL_UINT8:  
        return boost::shared_ptr<GlTextureHandleBase> ( new GlTextureHandle<vw::PixelRGBA<vw::uint8> >(m_rsrc, m_record) );
        break;
      case vw::VW_CHANNEL_UINT16:
        return boost::shared_ptr<GlTextureHandleBase> ( new GlTextureHandle<vw::PixelRGBA<vw::uint16> >(m_rsrc, m_record) );
        break;
      case vw::VW_CHANNEL_INT16: 
        return boost::shared_ptr<GlTextureHandleBase> ( new GlTextureHandle<vw::PixelRGBA<vw::int16> >(m_rsrc, m_record) );
        break;
      case vw::VW_CHANNEL_FLOAT32: 
        return boost::shared_ptr<GlTextureHandleBase> ( new GlTextureHandle<vw::PixelRGBA<vw::float32> >(m_rsrc, m_record) );
        break;
      default:
        vw::vw_throw(vw::IOErr() << "GlTextureHandle: generate() failed. Unknown channel type: " << channel_type << ".\n");
      }

    case vw::VW_PIXEL_SCALAR:
      switch(channel_type) {
      case vw::VW_CHANNEL_FLOAT32: 
        return boost::shared_ptr<GlTextureHandleBase> ( new GlTextureHandle<vw::float32>(m_rsrc, m_record) );
        break;
      default:
        vw::vw_throw(vw::IOErr() << "GlTextureHandle: generate() failed. Unknown channel type: " << channel_type << ".\n");
      }

    default: 
      vw::vw_throw(vw::IOErr() << "GlTextureHandle: generate() failed.  " 
                   << "Unknown pixel format: " << pixel_format << ".\n");
    }
    // Never reached
    return boost::shared_ptr<GlTextureHandleBase>();
  }
};

// --------------------------------------------------------------
//                     TextureRecord
// --------------------------------------------------------------

struct TextureRecord : public TextureRecordBase {
  vw::Cache::Handle<GlTextureGenerator> handle;
  virtual ~TextureRecord() {}
};

// --------------------------------------------------------------
//                     GlTextureCache
// --------------------------------------------------------------

class GlTextureCache {
  vw::Thread *m_texture_fetch_thread;
  vw::Cache *m_gl_texture_cache_ptr;
  boost::shared_ptr<TextureFetchTask> m_texture_fetch_task;

  // We purge the outgoing request queue whenever there is a change
  // in LOD so that we can immediately begin serving tiles at the
  // new level of detail.  We keep track of the previous LOD here.
  int m_previous_lod;
  
  std::vector<boost::shared_ptr<TextureRecord> > m_texture_records;

  // Communication to/from the texture caching thread is handled using
  // a pair of request queues.  These queues are locked with a mutex.
  std::list<boost::shared_ptr<TextureRecord> > m_outgoing_requests, m_incoming_requests;
  vw::Mutex m_outgoing_requests_mutex, m_incoming_request_mutex;
 
  // Fetch the "best match" texture handle from the list of texture
  // records.  The best match has the bounding box requested and the
  // largest lod possible that does not exceed the requested lod.
  boost::shared_ptr<TextureRecord> get_record(vw::BBox2i bbox, int lod);

public:

  // Constructor/destructor
  GlTextureCache();
  ~GlTextureCache();

  // Register a texture with the cache.    
  void register_texture(boost::shared_ptr<vw::ImageResource> const& rsrc, 
                        vw::BBox2i bbox, int lod, 
                        CachedTextureRenderer* parent) {
    VW_ASSERT(lod >= 0, vw::ArgumentErr() << "GlTextureGenerator : lod must be greater than or equal to 0.");

    vw::vw_out(vw::VerboseDebugMessage) << "GlTextureCache::register_texture() registered bbox " << bbox << " @ lod " << lod << "\n";
    TextureRecord* new_record = new TextureRecord();
    boost::shared_ptr<TextureRecord> new_record_ptr(new_record);
    new_record->bbox = bbox;
    new_record->lod = lod;
    new_record->parent = parent;
    new_record->texture_id = 0;
    new_record->handle = m_gl_texture_cache_ptr->insert( GlTextureGenerator(rsrc, new_record_ptr) );
    m_texture_records.push_back( new_record_ptr );
  }

  // Fetch a texture from the cache.  This is a non-blocking call that
  // will immediately return the GL texture id of the texture *if it
  // is available*.  If the texture is not available, this function
  // will add it to the queue to be rendered by the texture fetch
  // thread and return 0 immediately.
  GLuint get_texture_id(vw::BBox2i bbox, int lod);
};

#endif // __VWV_GL_PREVIEW_WIDGET_TEXTURE_H__



