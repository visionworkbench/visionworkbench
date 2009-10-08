// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#ifndef __VW_PLATE_PLATEMANAGER_H__
#define __VW_PLATE_PLATEMANAGER_H__

#include <vw/Image.h>
#include <vw/Math/Vector.h>
#include <vw/Cartography/GeoReference.h>
#include <vw/Cartography/ToastTransform.h>

#include <vw/Plate/Index.h>
#include <vw/Plate/Blob.h>
#include <vw/Plate/PlateFile.h>
#include <vw/Plate/ImageComposite.h>

#include <vector>
#include <fstream>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/convenience.hpp>

// Protocol Buffer
#include <vw/Plate/ProtoBuffers.pb.h>

namespace fs = boost::filesystem;

namespace vw {
namespace platefile {

  // -------------------------------------------------------------------------
  //                            PLATE MANAGER
  // -------------------------------------------------------------------------
  
  class PlateManager {

  protected:
    boost::shared_ptr<PlateFile> m_platefile;
    FifoWorkQueue m_queue;

    // The Tile Entry is used to keep track of the bounding box of
    // tiles and their location in the grid.
    struct TileInfo {
      int i, j;
      BBox2i bbox;
      TileInfo(int i, int j, BBox2i const& bbox) : i(i), j(j), bbox(bbox) {}
    };

    // -------------------------------------------------------------------------
    //                       PLATE FILE TASKS
    // -------------------------------------------------------------------------
    
    template <class ViewT>
    class WritePlateFileTask : public Task {
      boost::shared_ptr<PlateFile> m_platefile;
      TileInfo m_tile_info;
      int m_depth;
      ViewT const& m_view;
      
    public:
      WritePlateFileTask(boost::shared_ptr<PlateFile> platefile, 
                         TileInfo const& tile_info, 
                         int depth, ImageViewBase<ViewT> const& view) : 
        m_platefile(platefile), m_tile_info(tile_info), m_depth(depth), 
        m_view(view.impl()) {}
      
      virtual ~WritePlateFileTask() {}
      virtual void operator() () { 
        std::cout << "\t    Generating tile: [ " << m_tile_info.j << " " << m_tile_info.i 
                  << " @ level " <<  m_depth << "]    BBox: " << m_tile_info.bbox << "\n";

        ImageView<typename ViewT::pixel_type> new_data = crop(m_view, m_tile_info.bbox);
        ImageView<typename ViewT::pixel_type> old_data(new_data.cols(), new_data.rows());
        try {
          m_platefile->read(old_data, m_tile_info.i, m_tile_info.j, m_depth);
        } catch (TileNotFoundErr &e) { 
          // Do nothing... we already have a default constructed empty image above! 
        }
        
        VW_ASSERT(old_data.cols() == new_data.cols() && old_data.rows() == new_data.rows(),
                  LogicErr() << "WritePlateFileTask::operator() -- new tile dimensions do not " 
                  << "match old tile dimensions.");

        vw::platefile::CompositeView<typename ViewT::pixel_type> composite;
        composite.insert(old_data, 0, 0);
        composite.insert(new_data, 0, 0);
        composite.prepare();

        ImageView<typename ViewT::pixel_type> composite_tile = composite;
        if( ! is_transparent(composite_tile) ) 
          m_platefile->write(composite_tile, m_tile_info.i, m_tile_info.j, m_depth);
      }
    };

  public:
  
    PlateManager(boost::shared_ptr<PlateFile> platefile, int num_threads) : 
      m_queue(num_threads), m_platefile(platefile) {}
    
    /// Destructor
    virtual ~PlateManager() {}

    // Read a previously-written tile in from disk.  Cache the most
    // recently accessed tiles, since each will be used roughly four
    // times.
    virtual ImageView<PixelRGBA<uint8> > load_tile( int32 level, int32 x, int32 y ) = 0;

    // Access the tile at the top of the tree.  This will cause a
    // cascade of MipMapping requests that will cause intermediate
    // tiles to be written (via load_tile), and ultimately bring the
    // entire tree up to date.
    virtual void mipmap() { 
      ImageView<PixelRGBA<uint8> > tile = this->load_tile(0,0,0); 
      if (tile && !is_transparent(tile))
        m_platefile->write(tile,0,0,0);
    }
  };


}} // namespace vw::plate

#endif // __VW_PLATE_PLATEMANAGER_H__
