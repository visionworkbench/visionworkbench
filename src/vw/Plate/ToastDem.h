// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifndef __VW_PLATE_TOAST_DEM_H__
#define __VW_PLATE_TOAST_DEM_H__

#include <vw/Plate/FundamentalTypes.h>
#include <boost/shared_array.hpp>

namespace vw {
namespace platefile {

  class PlateFile;

  // -------------------------------------------------------------------------
  //                              TOAST DEM UTILITIES
  // -------------------------------------------------------------------------

  struct ToastDemWriter {
    virtual ~ToastDemWriter() {}
    virtual void operator()(const boost::shared_array<uint8> data, size_t data_size,
                            int32 dem_col, int32 dem_row, int32 dem_level,
                            Transaction output_transaction_id) const = 0;
  };

  // Writer is called for every DEM tile (probably will be called more
  // than once for a given image tile)
  bool make_toast_dem_tile(const ToastDemWriter& writer,
                           const PlateFile& platefile, int32 col, int32 row, int32 level,
                           int32 level_difference,
                           TransactionOrNeg input_transaction_id, Transaction output_transaction_id);

  void save_toast_dem_tile(std::string base_output_name, boost::shared_ptr<PlateFile> platefile,
                           int32 col, int32 row, int32 level, int32 level_difference, TransactionOrNeg transaction_id);

  boost::shared_array<uint8> toast_dem_null_tile(uint64& output_tile_size);


}} // namespace vw::platefile

#endif // __VW_PLATE_TOAST_DEM_H__
