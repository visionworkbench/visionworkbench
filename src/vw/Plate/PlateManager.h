// __BEGIN_LICENSE__
// Copyright (C) 2006-2009 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifndef __VW_PLATE_PLATEMANAGER_H__
#define __VW_PLATE_PLATEMANAGER_H__

#include <vw/Image/ImageView.h>
#include <vw/Image/PixelTypes.h>
#include <vw/Math/Vector.h>
#include <vw/Cartography/GeoReference.h>
#include <vw/Cartography/ToastTransform.h>
#include <vw/Mosaic/ImageComposite.h>

#include <vw/Plate/Index.h>
#include <vw/Plate/Blob.h>
#include <vw/Plate/PlateFile.h>

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

  template <class ImplT>
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
    //                       FIND_NEAREST_TILE
    // -------------------------------------------------------------------------

    // // This function searches up the quad tree, looking for the
    // // highest resolution tiel that is available.  Once it finds this
    // // tile, it crops and supersamples the appropriate region and
    // // returns it.
    // template <class PixelT>
    // ImageView<PixelT> find_nearest_tile()

    // -------------------------------------------------------------------------
    //                       PLATE FILE TASKS
    // -------------------------------------------------------------------------
    
    template <class ViewT>
    class WritePlateFileTask : public Task {
      boost::shared_ptr<PlateFile> m_platefile;
      int m_read_transaction_id, m_write_transaction_id;
      TileInfo m_tile_info;
      int m_depth;
      ViewT const& m_view;
      bool m_verbose;
      SubProgressCallback m_progress;
      
    public:
      WritePlateFileTask(boost::shared_ptr<PlateFile> platefile, 
                         int read_transaction_id,
                         int write_transaction_id,
                         TileInfo const& tile_info, 
                         int depth, ImageViewBase<ViewT> const& view,
                         bool verbose, int total_num_blocks, 
                         const ProgressCallback &progress_callback = ProgressCallback::dummy_instance()) :
        m_platefile(platefile), m_read_transaction_id(read_transaction_id),
        m_write_transaction_id(write_transaction_id),
        m_tile_info(tile_info), m_depth(depth), m_view(view.impl()), 
        m_verbose(verbose), m_progress(progress_callback,0.0,1.0/float(total_num_blocks)) {}
      
      virtual ~WritePlateFileTask() {}
      virtual void operator() () { 
        if (m_verbose) 
          std::cout << "\t    Generating tile: [ " << m_tile_info.j << " " << m_tile_info.i 
                    << " @ level " <<  m_depth << "]    BBox: " << m_tile_info.bbox << "\n";
        ImageView<typename ViewT::pixel_type> new_data = crop(m_view, m_tile_info.bbox);

        // No computation is done here.  This is a shallow copy.
        ImageView<typename ViewT::pixel_type> composite_tile = new_data;
        
        // If the new_data is not transparent, then we need to go
        // fetch the existing data so that we can composite the two
        // tiles.  This search begins at the current level, but if no
        // tile exists at that level, then we search up the tree for a
        // lower resolution tile that we can safely supersample and
        // crop.
        if( ! is_opaque(new_data) ) {
          bool found = false;
          IndexRecord closest_record;
          int search_depth = m_depth;
          int search_col = m_tile_info.i;
          int search_row = m_tile_info.j;

          // Search up the tree, looking for the nearest tile that
          // exists on top of which we can composite the new data.
          while (!found && search_depth >= 0) {
            try {
              closest_record = m_platefile->read_record(search_col, search_row, 
                                                        m_depth, m_read_transaction_id);
              found = true;
            } catch (TileNotFoundErr &e) { 
              // If the tile is not found, the search continues at the next level
              --search_depth;
              search_col /= 2;
              search_row /= 2;
            }
          }
          
          //          if (search_depth >= 0) {
          if (search_depth == m_depth) {
            ImageView<typename ViewT::pixel_type> old_data(new_data.cols(), new_data.rows());
            m_platefile->read(old_data, m_tile_info.i, m_tile_info.j, 
                              search_depth, m_read_transaction_id);
            
            // Supersample and crop

            
            VW_ASSERT(old_data.cols() == new_data.cols() && old_data.rows() == new_data.rows(),
                      LogicErr() << "WritePlateFileTask::operator() -- " 
                      << "new tile dimensions do not match old tile dimensions.");

            vw::mosaic::ImageComposite<typename ViewT::pixel_type> composite;
            composite.insert(old_data, 0, 0);
            composite.insert(new_data, 0, 0);
            composite.set_draft_mode( true );
            composite.prepare();
            composite_tile = composite;
          } else { // if (search_depth >= 0) {
            // std::cout << "Didn't find tile: [ " << m_tile_info.i << " " 
            //           << m_tile_info.j << " @ " << m_depth << "]         but did at [" 
            //           << search_col << " " << search_row << " @ " << search_depth << "]\n";
          }
        }

        if( ! is_transparent(composite_tile) ) 
          m_platefile->write(composite_tile, m_tile_info.i, m_tile_info.j, 
                             m_depth, m_write_transaction_id);

        // Report progress
        m_progress.report_incremental_progress(1.0);
      }
    };

  public:
  
    PlateManager(boost::shared_ptr<PlateFile> platefile, int num_threads) : 
      m_platefile(platefile), m_queue(num_threads)  {}
    
    /// Destructor
    virtual ~PlateManager() {}

    /// \cond INTERNAL
    // Methods to access the derived type
    inline ImplT& impl() { return static_cast<ImplT&>(*this); }
    inline ImplT const& impl() const { return static_cast<ImplT const&>(*this); }
    /// \endcond


    // Access the tile at the top of the tree.  This will cause a
    // cascade of MipMapping requests that will cause intermediate
    // tiles to be written (via load_tile), and ultimately bring the
    // entire tree up to date.
    void mipmap(int read_transaction_id, int write_transaction_id) { 
      PixelFormatEnum pixel_format = m_platefile->pixel_format();
      ChannelTypeEnum channel_type = m_platefile->channel_type();

      // It would seem that these need to be initialized here outside
      // of the case statement, otherwise C++ complains.
      ImageView<PixelGray<uint8> > gray8_tile; 
      ImageView<PixelGray<int16> > gray16_tile;
      ImageView<PixelGray<float> > gray_float32_tile; 
      ImageView<PixelGrayA<uint8> > graya8_tile;
      ImageView<PixelGrayA<float> > graya_float32_tile;
      ImageView<PixelRGBA<uint8> > rgba8_tile; 
      ImageView<PixelRGBA<uint16> > rgba16_tile; 

      switch(pixel_format) {
      case VW_PIXEL_GRAY:
        switch(channel_type) {
        // case VW_CHANNEL_UINT8:  
        //   // impl().load_tile(gray8_tile,0,0,0,read_transaction_id,write_transaction_id);
        //   // if (gray8_tile && !is_transparent(gray8_tile))
        //   //   m_platefile->write(gray8_tile,0,0,0,write_transaction_id);
        //   break;
        // case VW_CHANNEL_INT16:  
        //   // impl().load_tile(gray16_tile,0,0,0,read_transaction_id,write_transaction_id); 
        //   // if (gray16_tile && !is_transparent(gray16_tile))
        //   //   m_platefile->write(gray16_tile,0,0,0,write_transaction_id);
        //   break;
        // case VW_CHANNEL_FLOAT32:  
        //   // impl().load_tile(gray_float32_tile,0,0,0,read_transaction_id,write_transaction_id); 
        //   // if (gray_float32_tile && !is_transparent(gray_float32_tile))
        //   //   m_platefile->write(gray_float32_tile,0,0,0,write_transaction_id);
        //   break;
        default:
          vw_throw(ArgumentErr() << "Platefile contains a channel type not supported by image2plate.\n");
        }
        break;
      case VW_PIXEL_GRAYA:
        switch(channel_type) {
        case VW_CHANNEL_UINT8:  
          impl().load_tile(graya8_tile,0,0,0,read_transaction_id,write_transaction_id); 
          if (graya8_tile && !is_transparent(graya8_tile))
            m_platefile->write(graya8_tile,0,0,0,write_transaction_id);
          break;
        case VW_CHANNEL_FLOAT32:  
          impl().load_tile(graya_float32_tile,0,0,0,read_transaction_id,write_transaction_id); 
          if (graya_float32_tile && !is_transparent(graya_float32_tile))
            m_platefile->write(graya_float32_tile,0,0,0,write_transaction_id);
          break;
        default:
          vw_throw(ArgumentErr() << "Platefile contains a channel type not supported by image2plate.\n");
        }
        break;
      case VW_PIXEL_RGB:
      case VW_PIXEL_RGBA:
      default:
        switch(channel_type) {
        case VW_CHANNEL_UINT8:  
          impl().load_tile(rgba8_tile,0,0,0,read_transaction_id,write_transaction_id); 
          if (rgba8_tile && !is_transparent(rgba8_tile))
            m_platefile->write(rgba8_tile,0,0,0,write_transaction_id);
          break;
        case VW_CHANNEL_UINT16:  
          impl().load_tile(rgba16_tile,0,0,0,read_transaction_id,write_transaction_id); 
          if (rgba16_tile && !is_transparent(rgba16_tile))
            m_platefile->write(rgba16_tile,0,0,0,write_transaction_id);
          break;
        default:
          vw_throw(ArgumentErr() << "Platefile contains a channel type not supported by image2plate.\n");
        }
        break;
      }
      
    }
  };


}} // namespace vw::plate

#endif // __VW_PLATE_PLATEMANAGER_H__
