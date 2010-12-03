// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifndef __VW_PLATE_TOAST_PLATEMANAGER_H__
#define __VW_PLATE_TOAST_PLATEMANAGER_H__

#include <vw/Plate/PlateManager.h>

namespace vw {
namespace platefile {

  template <class PixelT>
  class ToastPlateManager : public PlateManager<PixelT> {
  protected:

    // Transforms an input image with input georef to toast projection
    // at the most ideal matching pyramid level.
    virtual void transform_image( cartography::GeoReference const& georef,
                                  ImageViewRef<PixelT>& image,
                                  TransformRef& txref, int& level ) const;

 public:
    ToastPlateManager(boost::shared_ptr<PlateFile> platefile) :
      PlateManager<PixelT>(platefile) {}

    // Non-working .. just matching base type
    cartography::GeoReference georeference( int level ) const;

    void generate_mipmap_tile(ImageView<PixelT>& dest,
        int col, int row, int level, TransactionOrNeg transaction_id, bool preblur) const;

    ImageView<PixelT> fetch_child_tile(
        int x, int y, int level, TransactionOrNeg transaction_id) const;
  };

}} // namespace vw::plate

#endif // __VW_PLATE_TOAST_PLATEMANAGER_H__
