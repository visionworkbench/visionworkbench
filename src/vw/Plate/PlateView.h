// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifndef __VW_PLATE_PLATEVIEW__
#define __VW_PLATE_PLATEVIEW__

#include <vw/Image/ImageViewBase.h>
#include <vw/Image/ImageViewRef.h>
#include <vw/Plate/PlateFile.h>
#include <vw/Core/Cache.h>

namespace vw {
namespace platefile {

  /// An image view for accessing tiles from a plate file.  Tiles are
  /// cached by this view to increase read speeds.
  template <class PixelT>
  class PlateView : public ImageViewBase<PlateView<PixelT> > {
    boost::shared_ptr<PlateFile> m_platefile;
    int m_current_depth;

  public:
    typedef PixelT pixel_type;
    typedef PixelT result_type;
    typedef ProceduralPixelAccessor<PlateView> pixel_accessor;

    PlateView(std::string plate_filename) {
      
      // Open the platefile
      m_platefile.reset( new PlateFile(plate_filename) );
      vw_out(InfoMessage, "plate") << "PlateView -- opened platefile \"" << plate_filename << "\"\n";
      m_current_depth = m_platefile->depth();
    }

    // Standard ImageView interface methods
    int32 cols() const { 
      return pow(2,m_platefile->depth()) * m_platefile->default_tile_size(); 
    }

    int32 rows() const {
      return pow(2,m_platefile->depth()) * m_platefile->default_tile_size(); 
    }

    int32 planes() const { return 1; }

    void set_depth(int depth) { 
      if (depth < 0 || depth > max_depth()) 
        vw_throw(ArgumentErr() << "PlateView::set_depth() -- invalid depth");
      m_current_depth = depth;
    }

    int max_depth() const {
      return m_platefile->depth();
    }

    inline pixel_accessor origin() const { return pixel_accessor( *this, 0, 0 ); }

    inline pixel_type operator()(float x, float y, int32 p = 0) const {
      vw_throw(NoImplErr() << "PlateView::operator() -- not yet implemented.");
    }

    

    /// \cond INTERNAL
    typedef CropView<TransformView<InterpolationView<EdgeExtensionView<CropView<ImageView<pixel_type> >, ConstantEdgeExtension>, BilinearInterpolation>, ResampleTransform> > prerasterize_type;
    inline prerasterize_type prerasterize(BBox2i bbox) const {

      //std::cout << "\nRasterizing bbox: " << bbox << "\n";
      
      // Compute the bounding box at the current depth.
      int depth_difference = m_platefile->depth() - m_current_depth;
      BBox2i level_bbox = bbox / pow(2,depth_difference);

      //      std::cout << "\tlevel bbox: " << level_bbox << "\n";

      // Grow that bounding box to align with tile boundaries
      BBox2i aligned_level_bbox = level_bbox;
      aligned_level_bbox.min() = ( (level_bbox.min() / m_platefile->default_tile_size()) 
                                  * m_platefile->default_tile_size() );
      aligned_level_bbox.max().x() = ( ceilf( float(level_bbox.max().x()) / 
                                              m_platefile->default_tile_size() ) 
                                       * m_platefile->default_tile_size() );
      aligned_level_bbox.max().y() = ( ceilf( float(level_bbox.max().y()) / 
                                              m_platefile->default_tile_size() )
                                       * m_platefile->default_tile_size() );

      //      std::cout << "\tAligned bbox: " << aligned_level_bbox << "\n";

      // Create an image of the appropriate size to rasterize tiles into.
      ImageView<pixel_type> level_image(aligned_level_bbox.width(), aligned_level_bbox.height());
      
      // Access the tiles needed for this level and mosaic them into the level_image
      int tile_y = aligned_level_bbox.min().y() / m_platefile->default_tile_size();
      int dest_row = 0;
      while ( tile_y < aligned_level_bbox.max().y() / m_platefile->default_tile_size() &&
              tile_y < pow(2, m_current_depth) ) {

        int tile_x = aligned_level_bbox.min().x() / m_platefile->default_tile_size();
        int dest_col = 0;
        while ( tile_x < aligned_level_bbox.max().x() / m_platefile->default_tile_size() &&
                tile_x < pow(2, m_current_depth)) {
          BBox2i tile_bbox(dest_col, dest_row,
                           m_platefile->default_tile_size(), 
                           m_platefile->default_tile_size());

          // std::cout << "\t  Load " << tile_x << " " << tile_y 
          //           << " @ " << m_current_depth << "\n";
          ImageView<PixelT> tile;
          IndexRecord rec;
          try {
            rec = m_platefile->read_record(tile_x, tile_y, m_current_depth);

            m_platefile->read(tile, tile_x, tile_y, m_current_depth, 0);
            crop(level_image, tile_bbox) = tile;

          } catch (TileNotFoundErr &e) {
            // Do nothing... we'll simply skip an empty image below
          } 

          ++tile_x;
          dest_col += m_platefile->default_tile_size();
        }
        ++tile_y;
        dest_row += m_platefile->default_tile_size();
      }



      // Crop the output to the original requested bbox, resample it, and return.
      BBox2i output_bbox( level_bbox.min().x() - aligned_level_bbox.min().x(),
                          level_bbox.min().y() - aligned_level_bbox.min().y(),
                          level_bbox.width(), level_bbox.height() );
      // std::cout << "\tlevel_bbox: " << level_bbox << "\n";
      // std::cout << "\toutput_bbox: " << output_bbox << "\n";
      // std::cout << "\tdepth_dif: " << depth_difference << "\n";;

      // std::ostringstream ostr; 
      // ostr << "level_" << bbox.min().x() << "_" << bbox.min().y() << ".tif";
      // write_image(ostr.str(), level_image);

      // write_image(ostr.str(), transform( crop(level_image, output_bbox), 
      //                                    ResampleTransform( pow(2, depth_difference), pow(2, depth_difference) ),
      //                                    ConstantEdgeExtension(), BilinearInterpolation() ) );

      return crop(transform( crop(level_image, output_bbox), 
                             ResampleTransform( pow(2, depth_difference), 
                                                pow(2, depth_difference) ),
                             ConstantEdgeExtension(), BilinearInterpolation() ),
                  BBox2i(-bbox.min().x(), -bbox.min().y(), bbox.width(), bbox.height()));
    }

    template <class DestT> inline void rasterize(DestT const& dest, BBox2i bbox) const {
      vw::rasterize(prerasterize(bbox), dest, bbox);
    }
    /// \endcond

  };

}} // namespace vw::platefile

#endif // __VW_PLATE_PLATEVIEW__
