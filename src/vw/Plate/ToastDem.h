// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifndef __VW_PLATE_TOAST_DEM_H__
#define __VW_PLATE_TOAST_DEM_H__

#include <vw/Plate/PlateFile.h>

namespace vw {
namespace platefile {

  // -------------------------------------------------------------------------
  //                              TOAST DEM UTILITIES
  // -------------------------------------------------------------------------

  struct ToastDemWriter {
    virtual ~ToastDemWriter() {}
    virtual void operator()(const boost::shared_array<uint8> data, uint64 data_size, int32 dem_col, int32 dem_row, int32 dem_level, int32 transaction_id) const = 0;
  };

  // Writer is called for every DEM tile (probably will be called more
  // than once for a given image tile)
  bool make_toast_dem_tile(const ToastDemWriter& writer,
      const PlateFile& platefile, int32 col, int32 row, int32 level, int32 transaction_id);

  void save_toast_dem_tile(std::string base_output_name, boost::shared_ptr<PlateFile> platefile,
                           int32 col, int32 row, int32 level, int32 transaction_id);


}} // namespace vw::platefile

#endif // __VW_PLATE_TOAST_DEM_H__
