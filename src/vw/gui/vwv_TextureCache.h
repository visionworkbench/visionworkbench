// __BEGIN_LICENSE__
// 
// Copyright (C) 2008 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration
// (NASA).  All Rights Reserved.
// 
// Copyright 2008 Carnegie Mellon University. All rights reserved.
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
#include <vw/Image/ImageViewRef.h>
#include <vw/Math/BBox.h>

// Forward Declarations
class CachedTextureRenderer;
class TextureFetchTask;
class TextureRequest;
class TextureRecord;

struct TextureRecordBase {
  vw::BBox2i bbox; 
  int lod;
  GLuint texture_id;
  CachedTextureRenderer* parent;
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
  virtual GLuint allocate_texture(vw::ImageView<vw::PixelRGBA<float> > const block) = 0;
  virtual void deallocate_texture(GLuint texture_id) = 0;

  
  virtual ~CachedTextureRenderer() { m_incoming_requests.clear(); }

  virtual GLuint request_allocation(boost::shared_ptr<TextureRecordBase> texture_record, 
                                    vw::ImageView<vw::PixelRGBA<float> > const block);

  virtual void request_deallocation(boost::shared_ptr<TextureRecordBase> texture_record);

  virtual void process_allocation_request();
};

// --------------------------------------------------------------
//                     GlTextureHandle
// --------------------------------------------------------------

class GlTextureHandle {
  boost::shared_ptr<TextureRecordBase> m_record;
  
public:
  template <class ViewT>
  GlTextureHandle(vw::ImageViewBase<ViewT> const& image, boost::shared_ptr<TextureRecordBase> record) : m_record(record) {

    // Rasterize the requested block of memory.
    vw::ImageView<typename ViewT::pixel_type> cropped = crop(image, record->bbox);
    vw::ImageView<typename ViewT::pixel_type> block = subsample(cropped, pow(2,record->lod));
    
    m_record->parent->request_allocation(record, block);
    vw::vw_out(vw::VerboseDebugMessage) << "GlTextureHandle requesting allocation (" 
                                        << m_record->texture_id << ") -- " 
                                        << m_record->bbox << " @ " << m_record->lod << "\n";
  }

  GLuint get_texture_id() const { return m_record->texture_id; }

  ~GlTextureHandle() {
    vw::vw_out(vw::VerboseDebugMessage) << "-> GlTextureHandle requesting decallocation (" << m_record->texture_id << ")\n";
    m_record->parent->request_deallocation(m_record);
  }
};

// --------------------------------------------------------------
//                     GlTextureGenerator
// --------------------------------------------------------------

class GlTextureGenerator {
  vw::ImageViewRef<vw::PixelRGBA<float> > m_image;
  boost::shared_ptr<TextureRecordBase> m_record;

public:
  typedef GlTextureHandle value_type;
  
  template <class ViewT>
  GlTextureGenerator( vw::ImageViewBase<ViewT> const &image, boost::shared_ptr<TextureRecordBase> record) :
    m_image( image.impl() ), m_record(record) {}

  size_t size() const {
    return m_record->bbox.width()/pow(2,m_record->lod) * 
      m_record->bbox.height()/pow(2,m_record->lod) * 
      m_image.planes() * m_image.channels() 
      * 2;  // <-- The textures are stored as 16-bit half's on the GPU, so they are 2 bytes per channel
  }
  
  boost::shared_ptr<GlTextureHandle> generate() const {
    return boost::shared_ptr<GlTextureHandle> ( new GlTextureHandle(m_image, m_record) );
  }
};

// --------------------------------------------------------------
//                     TextureRecord
// --------------------------------------------------------------

struct TextureRecord : public TextureRecordBase {
  vw::Cache::Handle<GlTextureGenerator> handle;
};

// --------------------------------------------------------------
//                     GlTextureCache
// --------------------------------------------------------------

class GlTextureCache {
  vw::Thread *m_texture_fetch_thread;
  vw::Cache *m_gl_texture_cache_ptr;
  boost::shared_ptr<TextureFetchTask> m_texture_fetch_task;
  int m_previous_lod;

  std::vector<boost::shared_ptr<TextureRecord> > m_texture_records;
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
  template<class ViewT>
  void register_texture(vw::ImageViewBase<ViewT> const& image, 
                        vw::BBox2i bbox, int lod, 
                        CachedTextureRenderer* parent) {
    VW_ASSERT(lod >= 0, vw::ArgumentErr() << "GlTextureGenerator : lod must be greater than or equal to 0.");

    vw::vw_out(vw::VerboseDebugMessage) << "GlTextureCache::register_texture() registered bbox " << bbox << " @ lod " << lod << "\n";
    boost::shared_ptr<TextureRecord> new_record( new TextureRecord() );
    new_record->bbox = bbox;
    new_record->lod = lod;
    new_record->parent = parent;
    new_record->texture_id = 0;
    new_record->handle = m_gl_texture_cache_ptr->insert( GlTextureGenerator(image, new_record) );
    m_texture_records.push_back( new_record );
  }

  // Fetch a texture from the cache.  This is a non-blocking call that
  // will immediately return the GL texture id of the texture *if it
  // is available*.  If the texture is not available, this function
  // will add it to the queue to be rendered by the texture fetch
  // thread and return 0 immediately.
  GLuint get_texture_id(vw::BBox2i bbox, int lod);
};

#endif // __VWV_GL_PREVIEW_WIDGET_TEXTURE_H__



