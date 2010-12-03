// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifndef __VW_PLATE_PLATE_CARREE_PLATEMANAGER_H__
#define __VW_PLATE_PLATE_CARREE_PLATEMANAGER_H__

#include <vw/Plate/PlateManager.h>

namespace vw {
namespace platefile {

  template <class PixelT>
  class PlateCarreePlateManager : public PlateManager<PixelT> {
  protected:
    // Transforms an input image with input georef to plate carree at
    // most ideal matching pyramid level.
    void transform_image( cartography::GeoReference const& georef,
                          ImageViewRef<PixelT>& image,
                          TransformRef& txref, int& level ) const;

    // This handles the border crossing of +-180 degrees.
    virtual void affected_tiles( BBox2i const& image_size,
                                 TransformRef const& tx, int tile_size,
                                 int level, std::list<TileInfo>& tiles ) const;

  public:
    PlateCarreePlateManager(boost::shared_ptr<PlateFile> platefile) :
      PlateManager<PixelT>(platefile) {}

    // Provide a georeference that represents a pyramid level
    cartography::GeoReference georeference( int level ) const;
  };

}} // namespace vw::plate

#endif // __VW_PLATE_PLATE_CARREE_PLATEMANAGER_H__
