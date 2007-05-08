// __BEGIN_LICENSE__
// 
// Copyright (C) 2006 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration
// (NASA).  All Rights Reserved.
// 
// Copyright 2006 Carnegie Mellon University. All rights reserved.
// 
// This software is distributed under the NASA Open Source Agreement
// (NOSA), version 1.3.  The NOSA has been approved by the Open Source
// Initiative.  See the file COPYING at the top of the distribution
// directory tree for the complete NOSA document.
// 
// THE SUBJECT SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY WARRANTY OF ANY
// KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT
// LIMITED TO, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL CONFORM TO
// SPECIFICATIONS, ANY IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR
// A PARTICULAR PURPOSE, OR FREEDOM FROM INFRINGEMENT, ANY WARRANTY THAT
// THE SUBJECT SOFTWARE WILL BE ERROR FREE, OR ANY WARRANTY THAT
// DOCUMENTATION, IF PROVIDED, WILL CONFORM TO THE SUBJECT SOFTWARE.
// 
// __END_LICENSE__
#include <vw/Math/Matrix.h>
#include <vw/Core/Exception.h>
#include <vw/Cartography/DiskImageResourceGeoReferenceHelperGDAL.h>
#include <vw/FileIO/DiskImageResourceGDAL.h>
#include <vw/Cartography/GeoReference.h>

//FIXME: can probably trim down these headers
// GDAL Headers
#include "gdal.h"
#include "gdal_priv.h"
#include "cpl_string.h"
#include "ogr_spatialref.h"
#include "ogr_api.h"

// For boost::split
#include <boost/algorithm/string.hpp>

namespace vw {
namespace cartography {

  /*static*/ void DiskImageResourceGeoReferenceHelperGDAL::read_georeference( FileMetadata* georef_, DiskImageResource* r_ ) {
    DiskImageResourceGDAL* r = (DiskImageResourceGDAL*)r_;
    GeoReference* georef = (GeoReference*)georef_;
    GDALDataset* dataset = (GDALDataset*)(r->dataset());
    
    if (!dataset) 
      vw_throw( LogicErr() << "DiskImageResourceGeoReferenceHelperGDAL: Could not read georeference. No file has been opened." );
    
    if( dataset->GetProjectionRef() != NULL ) {
      georef->set_wkt_str(dataset->GetProjectionRef());
    }
    
    double geo_transform[6];
    Matrix<double,3,3> transform;
    if( dataset->GetGeoTransform( geo_transform ) == CE_None ) {
      transform.set_identity();
      transform(0,0) = geo_transform[1];
      transform(0,1) = geo_transform[2];
      transform(0,2) = geo_transform[0];
      transform(1,0) = geo_transform[4];
      transform(1,1) = geo_transform[5];
      transform(1,2) = geo_transform[3];
      georef->set_transform(transform);

      // Determine the pixel interpretation for the image.  See the
      // comments in GeoReference for more information.  If nothing is
      // encoded in the file, the default is to assume PixelAsArea.
      georef->set_pixel_interpretation(GeoReference::PixelAsArea);
      char **metadata = dataset->GetMetadata();
      if( CSLCount(metadata) > 0 ) {
        for( int i = 0; metadata[i] != NULL; i++ ) {
          std::vector<std::string> split_vec;
          boost::split(split_vec, metadata[i], boost::is_any_of("=") );
          if (split_vec[0] == GDALMD_AREA_OR_POINT && split_vec.size() >= 2)  
            if (boost::trim_copy(split_vec[1]) == GDALMD_AOP_POINT)
              georef->set_pixel_interpretation(GeoReference::PixelAsPoint);
        }
      }


    }
  }
  
  /*static*/ void DiskImageResourceGeoReferenceHelperGDAL::write_georeference( DiskImageResource* r_, FileMetadata const* georef_ ) {
    DiskImageResourceGDAL* r = (DiskImageResourceGDAL*)r_;
    GeoReference const* georef = (GeoReference const*)georef_;
    GDALDataset* dataset = (GDALDataset*)(r->dataset());
    
    if (!dataset) 
      vw_throw( LogicErr() << "DiskImageResourceGeoReferenceHelperGDAL: Could not write georeference. No file has been opened." );

    // Store the transform matrix
    double geo_transform[6] = { georef->transform()(0,2), georef->transform()(0,0), georef->transform()(0,1), 
                                georef->transform()(1,2), georef->transform()(1,0), georef->transform()(1,1) };
    dataset->SetGeoTransform( geo_transform );
    dataset->SetProjection( georef->wkt_str().c_str() );

    // Set the pixel interpretation for the image.  See the comments
    // in GeoReference for more information.
    if (georef->pixel_interpretation() == GeoReference::PixelAsArea) 
      dataset->SetMetadataItem(GDALMD_AREA_OR_POINT, GDALMD_AOP_AREA);
    else // PixelAsPoint
      dataset->SetMetadataItem(GDALMD_AREA_OR_POINT, GDALMD_AOP_POINT);
  }

}} // namespace vw::cartography
