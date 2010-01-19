// __BEGIN_LICENSE__
// Copyright (C) 2006-2009 United States Government as represented by
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

  void save_toast_dem_tile(std::string base_output_name, boost::shared_ptr<PlateFile> platefile, 
                           int32 col, int32 row, int32 level, int32 transaction_id);


}} // namespace vw::platefile

#endif // __VW_PLATE_TOAST_DEM_H__
