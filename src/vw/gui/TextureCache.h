// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifndef __VW_GUI_TEXTURE_CACHE_H__
#define __VW_GUI_TEXTURE_CACHE_H__

#include <vw/config.h>

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
#include <vw/Image/Algorithms.h>
#include <vw/Image/ImageResourceView.h>
#include <vw/Image/ImageViewRef.h>
#include <vw/Math/BBox.h>
#include <vw/gui/TileGenerator.h>
#include <vw/gui/Tree.h>

namespace vw {
namespace gui {

  // Forward Declarations
  class CachedTextureRenderer;
  class TextureRequest;
  class TextureFetchTask;

  struct TextureRecordBase {
    GLuint texture_id;
    CachedTextureRenderer* requestor;
    virtual ~TextureRecordBase() {}
  };

  // --------------------------------------------------------------
  //                     CachedTextureRenderer
  // --------------------------------------------------------------

  class CachedTextureRenderer {
    std::list<boost::shared_ptr<TextureRequest> > m_requests;
    vw::Mutex m_request_mutex;

  public:

    // These are defined in the subclass: vwv_GlPreviewWidget
    virtual GLuint allocate_texture(boost::shared_ptr<SrcImageResource> tile) = 0;
    virtual void deallocate_texture(GLuint texture_id) = 0;

    virtual ~CachedTextureRenderer() { m_requests.clear(); }

    virtual void request_allocation(boost::shared_ptr<TextureRecordBase> texture_record,
                                    boost::shared_ptr<SrcImageResource> tile);

    virtual void request_deallocation(boost::shared_ptr<TextureRecordBase> texture_record);

    virtual void process_allocation_requests();
  };

  // --------------------------------------------------------------
  //                     GlTextureHandle
  // --------------------------------------------------------------

  struct GlTextureHandleBase {
    virtual GLuint texture_id() const = 0;
    virtual ~GlTextureHandleBase() {}
  };

  class GlTextureHandle : public GlTextureHandleBase {
    boost::shared_ptr<TextureRecordBase> m_record;

  public:
    GlTextureHandle(boost::shared_ptr<TileGenerator> generator,
                    TileLocator tile_info,
                    boost::shared_ptr<TextureRecordBase> record) : m_record(record) {

      // Load the tile into memory by querying the TileGenerator.
      // This can take a little while, which is why we perform this
      // step in a seperate thread.
      boost::shared_ptr<SrcImageResource> tile = generator->generate_tile(tile_info);

      // For debugging:
      //
      // std::ostringstream ostr;
      // ostr << "debug_" << record->lod << "_" << record->bbox.min().x() << "_"
      //      << record->bbox.min().y() << ".tif";
      // vw::write_image(ostr.str(), tile);

      // Send a request to the OpenGL thread to load this tile as a
      // texture on the graphics card.  This will be done on the next
      // rendering pass.
      m_record->requestor->request_allocation( record, tile );
      vw_out(vw::VerboseDebugMessage, "gui") << "GlTextureHandle requesting allocation ("
                                             << m_record->texture_id << ") -- [ "
                                             << tile_info.col << " " << tile_info.row
                                             << " ] @ " << tile_info.level << "\n";
    }

    virtual GLuint texture_id() const { return m_record->texture_id; }

    virtual ~GlTextureHandle() {
      vw_out(vw::VerboseDebugMessage) << "-> GlTextureHandle requesting decallocation ("
                                      << m_record->texture_id << ")\n";

      // Send a request to the OpenGL thread to deallocate this
      // texture on the next rendering pass.
      m_record->requestor->request_deallocation(m_record);
    }
  };

  // --------------------------------------------------------------
  //                     GlTextureGenerator
  // --------------------------------------------------------------

  class GlTextureGenerator {
    boost::shared_ptr<TileGenerator> m_tile_generator;
    TileLocator m_tile_info;
    boost::shared_ptr<TextureRecordBase> m_record;

  public:
    typedef GlTextureHandleBase value_type;

    GlTextureGenerator( boost::shared_ptr<TileGenerator> const& generator,
                        TileLocator const& tile_info,
                        boost::shared_ptr<TextureRecordBase> record) :
      m_tile_generator( generator ), m_tile_info(tile_info), m_record(record) {
    }

    size_t size() const {
      int size =  m_tile_generator->tile_size()[0] * m_tile_generator->tile_size()[1] * num_channels(m_tile_generator->pixel_format()) * channel_size(m_tile_generator->channel_type());
      return size;
    }

    boost::shared_ptr<GlTextureHandleBase> generate() const {
      return boost::shared_ptr<GlTextureHandleBase> (
          new GlTextureHandle(m_tile_generator, m_tile_info, m_record) );
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

    // We purge the outgoing request queue whenever there is a change
    // in LOD so that we can immediately begin serving tiles at the
    // new level of detail.  We keep track of the previous LOD here.
    int m_previous_level;

    // Communication to/from the texture caching thread is handled using
    // a pair of request queues.  These queues are locked with a mutex.
    std::list<boost::shared_ptr<TextureRecord> > m_requests;
    vw::Mutex m_request_mutex;

    // We store texure records in a quad tree structure.  For now we are
    // going to use the tree structure provided by the plate module,
    // although that tree structure should maybe be moved to the core
    // module someday.
    boost::shared_ptr<gui::TreeNode<boost::shared_ptr<TextureRecord> > > m_texture_records;
    vw::Cache* m_gl_texture_cache_ptr;
    boost::shared_ptr<TextureFetchTask> m_texture_fetch_task;
    vw::Thread *m_texture_fetch_thread;

    // Shared ptr to the texture generator
    boost::shared_ptr<TileGenerator> m_tile_generator;

  public:

    // Constructor/destructor
    GlTextureCache(boost::shared_ptr<TileGenerator> tile_generator);
    ~GlTextureCache();

    // Get a handle on the generator being used to produce tiles.
    boost::shared_ptr<TileGenerator> tile_generator() const {
      return m_tile_generator;
    }

    // Clear all entries from the texture cache.
    void clear();

    // Fetch a texture from the cache.  This is a non-blocking call that
    // will immediately return the GL texture id of the texture *if it
    // is available*.  If the texture is not available, this function
    // will add it to the queue to be rendered by the texture fetch
    // thread and return 0 immediately.
    GLuint get_texture_id(vw::gui::TileLocator const& tile_info,
                          CachedTextureRenderer* requestor);
  };

}} // namespace vw::gui

#endif // __VW_GUI_TEXTURE_CACHE_H__
