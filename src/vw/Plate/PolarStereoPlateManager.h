// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#ifndef __VW_PLATE_POLAR_STEREO_PLATEMANAGER_H__
#define __VW_PLATE_POLAR_STEREO_PLATEMANAGER_H__

#include <vw/Plate/PlateManager.h>

namespace vw {
namespace platefile {

  template <class PixelT>
  class PolarStereoPlateManager : public PlateManager<PixelT> {
  protected:
    // Transforms an input image with input georef to polar stereo
    // graphic at most ideal matching pyramid level.
    void transform_image( cartography::GeoReference const& georef,
                          ImageViewRef<PixelT>& image,
                          TransformRef& txref, int& level ) const;

  public:
    PolarStereoPlateManager(boost::shared_ptr<PlateFile> platefile):
      PlateManager<PixelT>(platefile) {}

    // Provide a georeference that represents a pyramid level
    cartography::GeoReference georeference( int level, bool north_pole,
                                            cartography::Datum const& datum ) const;
    cartography::GeoReference georeference( int level ) const;

    /// This function generates a specific mipmap tile at the given
    /// col, row, and level, and transaction_id.
    void generate_mipmap_tile(int col, int row, int level,
                              TransactionOrNeg transaction_id, bool preblur) const;
  };

}} // end vw::platefile

#endif//__VW_PLATE_POLAR_STEREO_PLATEMANAGER_H__
