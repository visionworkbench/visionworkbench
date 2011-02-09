// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <vw/gui/TextureCache.h>
using namespace vw;
using namespace vw::gui;

// --------------------------------------------------------------
//                     TextureRequest
// --------------------------------------------------------------

struct vw::gui::TextureRequest {
  virtual ~TextureRequest() {}
  virtual void process_request() = 0;
};

// Allocate Texture Request
class AllocateTextureRequest : public TextureRequest {
  boost::shared_ptr<TextureRecordBase> m_record;
  boost::shared_ptr<SrcImageResource> m_tile;
  CachedTextureRenderer* m_parent;
public:
  AllocateTextureRequest( boost::shared_ptr<TextureRecordBase> texture_record,
                          boost::shared_ptr<SrcImageResource> tile,
                          CachedTextureRenderer* parent) :
    m_record(texture_record), m_tile(tile), m_parent(parent) {}
  virtual ~AllocateTextureRequest() {}

  virtual void process_request() {
    m_record->texture_id = m_parent->allocate_texture(m_tile);
  }
};

// Deallocate Texture Request
class DeallocateTextureRequest : public TextureRequest {
  boost::shared_ptr<TextureRecordBase> m_record;
  CachedTextureRenderer* m_parent;
public:
  DeallocateTextureRequest( boost::shared_ptr<TextureRecordBase> texture_record,
                            CachedTextureRenderer* parent) :
    m_record(texture_record), m_parent(parent) {}
  virtual ~DeallocateTextureRequest() {}

  virtual void process_request() {
    m_parent->deallocate_texture(m_record->texture_id);
    m_record->texture_id = 0;
  }
};

// --------------------------------------------------------------
//                     CachedTextureRenderer
// --------------------------------------------------------------

void CachedTextureRenderer::request_allocation(boost::shared_ptr<TextureRecordBase> texture_record,
                                               boost::shared_ptr<SrcImageResource> tile) {
  vw::Mutex::Lock lock(m_request_mutex);
  m_requests.push_back( boost::shared_ptr<TextureRequest>(new AllocateTextureRequest(texture_record, tile, this)) );
}

void CachedTextureRenderer::request_deallocation(boost::shared_ptr<TextureRecordBase> texture_record) {
  vw::Mutex::Lock lock(m_request_mutex);
  m_requests.push_back( boost::shared_ptr<TextureRequest>(new DeallocateTextureRequest(texture_record, this)) );
}

void CachedTextureRenderer::process_allocation_requests() {
  vw::Mutex::Lock lock(m_request_mutex);

  boost::shared_ptr<TextureRequest> r;

  while (!m_requests.empty()) {
    r = m_requests.front();
    m_requests.pop_front();
    r->process_request();
  }
}

// --------------------------------------------------------------
//                     TextureFetchTask
// --------------------------------------------------------------

class vw::gui::TextureFetchTask {
  bool terminate;
  vw::Mutex &m_request_mutex;
  std::list<boost::shared_ptr<TextureRecord> > &m_requests;

public:
  TextureFetchTask(vw::Mutex &request_queue_mutex,
                   std::list<boost::shared_ptr<TextureRecord> > &request_queue) :
    terminate(false), m_request_mutex(request_queue_mutex), m_requests(request_queue) {}

  void operator()() {
    while (!terminate) {

      // Drain the request queue
      bool found = false;
      boost::shared_ptr<TextureRecord> r;
      {
        vw::Mutex::Lock lock(m_request_mutex);
        if (!m_requests.empty()) {
          r = m_requests.front();
          m_requests.pop_front();
          found = true;
        }
      }

      if (found) {

        // Force the texture to regenerate.  Doing so will cause the
        // image tile to be loaded into memory, and then a texture
        // allocation request to be generated to be handled later by
        // the OpenGL thread.  This may cause one or more texture
        // deallocation requests to be produced as well if cache tile
        // need to be deallocated to make room for the new tile.
        (*(r->handle)).texture_id();

      } else {

        // If there were no requests, sleep for a short time.
        vw::Thread::sleep_ms(100);

      }
    }
  }

  void kill() { terminate = true; }
};


// --------------------------------------------------------------
//                     GlTextureCache
// --------------------------------------------------------------

vw::gui::GlTextureCache::GlTextureCache(boost::shared_ptr<TileGenerator> tile_generator) :
  m_tile_generator(tile_generator) {

  // Create the texture cache
  int gl_texture_cache_size = 256 * 1024 * 1024; // Use 128-MB of
                                                 // texture cache

  m_gl_texture_cache_ptr = new vw::Cache( gl_texture_cache_size );

  // Start the texture fetch thread
  m_texture_fetch_task.reset(new TextureFetchTask(m_request_mutex, m_requests));
  m_texture_fetch_thread = new vw::Thread( m_texture_fetch_task );

  // Create the texture record tree for storing cache handles and
  // other useful texture-related metadata.
  m_texture_records.reset( new gui::TreeNode<boost::shared_ptr<TextureRecord> >() );
  m_previous_level = 0;
}

vw::gui::GlTextureCache::~GlTextureCache() {
  // Stop the Texture Fetch thread
  m_texture_fetch_task->kill();
  m_texture_fetch_thread->join();
  delete m_texture_fetch_thread;

  // Free up remaining texture handles, and then the cache itself.
  m_requests.clear();
  delete m_gl_texture_cache_ptr;
}

void vw::gui::GlTextureCache::clear() {
  Mutex::Lock lock(m_request_mutex);
  m_requests.clear();

  // Delete all of the existing texture records
  m_texture_records.reset( new gui::TreeNode<boost::shared_ptr<TextureRecord> >() );
  m_previous_level = 0;
}

GLuint vw::gui::GlTextureCache::get_texture_id(vw::gui::TileLocator const& tile_info,
                                               CachedTextureRenderer* requestor) {
  // Bail early if the tile_info request is totally invalid.
  if (!tile_info.is_valid())
    return 0;

  // We purge the outgoing request queue whenever there is a change
  // in LOD so that we can immediately begin serving tiles at the
  // new level of detail.
  if (tile_info.level != m_previous_level) {
    Mutex::Lock lock(m_request_mutex);
    m_requests.clear();
    m_previous_level = tile_info.level;
  }

  try {

    // First, see if the texture record is already in the cache.
    boost::shared_ptr<TextureRecord> rec = m_texture_records->search(tile_info.col,
                                                                     tile_info.row,
                                                                     tile_info.level,
                                                                     tile_info.transaction_id,
                                                                     false);

    // If the shared pointer for this record is empty, then this node
    // was generated as part of a branch that supports a leaf node,
    // but this node does not itself contain any data yet.  To
    // populate it with data, we throw an exception to run the code
    // below.
    if ( !rec )
      vw_throw(gui::TileNotFoundErr() << "invalid record. regenerating...");

    // If the texture_id of this record is 0, then we need to send a
    // request to regenerate the texture.  It will get rendered in the
    // future after it has been loaded.
    if (rec->texture_id == 0) {
      Mutex::Lock lock(m_request_mutex);
      m_requests.push_back( rec );
      return 0;
    }

    // If the texture_id is valid, then the tile must be good.  Return
    // the texture_id to satisfy the request.
    return rec->texture_id;

  } catch (gui::TileNotFoundErr &e) {

    // If the tile isn't found or hasn't been properly initialized
    // yet, we need to add an entry to the cache and then cause the
    // tile to be generated.
    TextureRecord* new_record = new TextureRecord();
    boost::shared_ptr<TextureRecord> new_record_ptr(new_record);
    new_record->texture_id = 0;
    new_record->requestor = requestor;
    new_record->handle = m_gl_texture_cache_ptr->insert( GlTextureGenerator(m_tile_generator,
                                                                            tile_info,
                                                                            new_record_ptr) );

    // Place this cache handle into the tree for later access.
    m_texture_records->insert( new_record_ptr, tile_info.col, tile_info.row,
                               tile_info.level, tile_info.transaction_id );

    Mutex::Lock lock(m_request_mutex);
    m_requests.push_back( new_record_ptr );
    return 0;

  }

  return 0; // never reached
}
