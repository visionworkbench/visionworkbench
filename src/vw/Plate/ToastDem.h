// __BEGIN_LICENSE__
//  Copyright (c) 2006-2013, United States Government as represented by the
//  Administrator of the National Aeronautics and Space Administration. All
//  rights reserved.
//
//  The NASA Vision Workbench is licensed under the Apache License,
//  Version 2.0 (the "License"); you may not use this file except in
//  compliance with the License. You may obtain a copy of the License at
//  http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
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
    virtual void operator()(const boost::shared_array<uint8> data, uint64 data_size,
                            int32 dem_col, int32 dem_row, int32 dem_level) const = 0;
  };

  // Writer is called for every DEM tile (probably will be called more
  // than once for a given image tile)
  bool make_toast_dem_tile(const ToastDemWriter& writer,
                           const PlateFile& platefile, int32 col, int32 row, int32 level,
                           int32 level_difference, TransactionOrNeg input_transaction_id);

  void save_toast_dem_tile(std::string base_output_name, boost::shared_ptr<PlateFile> platefile,
                           int32 col, int32 row, int32 level, int32 level_difference, TransactionOrNeg transaction_id);

  boost::shared_array<uint8> toast_dem_null_tile(uint64& output_tile_size);


}} // namespace vw::platefile

#endif // __VW_PLATE_TOAST_DEM_H__
