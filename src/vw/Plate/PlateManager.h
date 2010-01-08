// __BEGIN_LICENSE__
// Copyright (C) 2006-2009 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#ifndef __VW_PLATE_PLATEMANAGER_H__
#define __VW_PLATE_PLATEMANAGER_H__

#include <vw/Plate/PlateFile.h>
#include <vw/Plate/ProtoBuffers.pb.h>

#include <list>

namespace vw {
namespace platefile {

  // The Tile Entry is used to keep track of the bounding box of
  // tiles and their location in the grid.
  struct TileInfo {
    int i, j;
    BBox2i bbox;
    TileInfo(int i, int j, BBox2i const& bbox) : i(i), j(j), bbox(bbox) {}
  };

  // Given a bbox, returns a list of smaller bboxes that perfectly
  // tile the space of the larger bbox.
  std::list<vw::BBox2i> bbox_tiles(vw::BBox2i const& bbox, int width, int height);

  // -------------------------------------------------------------------------
  //                              PLATE MANAGER
  // -------------------------------------------------------------------------

  class PlateManager {

  protected:
    boost::shared_ptr<PlateFile> m_platefile;

  public:

    PlateManager(boost::shared_ptr<PlateFile> platefile) : m_platefile(platefile) {}
    
    // ---------------------------- MIPMAPPING --------------------------------

    // mipmap() generates mipmapped (i.e. low resolution) tiles in the mosaic.
    //
    //   starting_level -- select the pyramid level from which to carry out mipmapping
    //   bbox -- bounding box (in terms of tiles) containing the tiles that need 
    //           to be mipmapped at starting_level.  Use to specify affected tiles.
    //   transaction_id -- transaction id to use when reading/writing tiles
    //
    void mipmap(int blob_id, int starting_level, BBox2i const& bbox, int transaction_id,
                const ProgressCallback &progress_callback = ProgressCallback::dummy_instance()) const;

    /// This function generates a specific mipmap tile at the given
    /// col, row, and level, and transaction_id.  It is left to a
    /// subclass of PlateManager to implement.
    virtual void generate_mipmap_tile(int blob_id, int col, int row, int level, int transaction_id) const = 0;

  };    

  // -------------------------------------------------------------------------
  //                            WRITE PLATEFILE TASK
  // -------------------------------------------------------------------------

  template <class ViewT>
  class WritePlateFileTask : public Task {
    int m_blob_id;
    boost::shared_ptr<PlateFile> m_platefile;
    int m_transaction_id;
    TileInfo m_tile_info;
    int m_level;
    ViewT const& m_view;
    bool m_verbose;
    SubProgressCallback m_progress;
      
  public:
    WritePlateFileTask(int blob_id,
                       boost::shared_ptr<PlateFile> platefile, 
                       int transaction_id,
                       TileInfo const& tile_info, 
                       int level, ImageViewBase<ViewT> const& view,
                       bool verbose, int total_num_blocks, 
                       const ProgressCallback &progress_callback = ProgressCallback::dummy_instance()) : m_blob_id(blob_id),
      m_platefile(platefile), m_transaction_id(transaction_id),
      m_tile_info(tile_info), m_level(level), m_view(view.impl()), 
      m_verbose(verbose), m_progress(progress_callback,0.0,1.0/float(total_num_blocks)) {}
      
    virtual ~WritePlateFileTask() {}
    virtual void operator() () {
      if (m_verbose) 
        std::cout << "\t    Generating tile: [ " << m_tile_info.j << " " << m_tile_info.i 
                  << " @ level " <<  m_level << "]    BBox: " << m_tile_info.bbox << "\n";

      // Generate the tile from the image data
      ImageView<typename ViewT::pixel_type> tile = crop(m_view, m_tile_info.bbox);
      
      // If this tile contains no data at all, then we bail early without
      // doing anything.
      if (is_transparent(tile)) {
        m_progress.report_incremental_progress(1.0);
        return;
      }
      
      // If the data is valid, we write it to the platefile.
      //
      // TODO: This is where we could strip the tile of its alpha
      // channel to save space in the placefile.  This will require a
      // view that strips off the alpha channel.
      m_platefile->write_update(m_blob_id, tile, 
                                m_tile_info.i, m_tile_info.j, m_level, m_transaction_id);
      m_progress.report_incremental_progress(1.0);
    }
  };

  

}} // namespace vw::plate

#endif // __VW_PLATE_PLATEMANAGER_H__
