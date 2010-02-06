// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <vw/Math/Matrix.h>
#include <vw/Core/Exception.h>
#include <vw/Cartography/GeoReferenceResourceGDAL.h>
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

  void read_gdal_georeference( GeoReference& georef, DiskImageResourceGDAL const& resource ) {
    boost::shared_ptr<GDALDataset>dataset = resource.get_dataset_ptr();
    if (!dataset) 
      vw_throw( LogicErr() << "read_gdal_georeference: Could not read georeference. No file has been opened." );
    
    // Pull the projection and datum information out of the file if available
    if( dataset->GetProjectionRef() != NULL ) {
      georef.set_wkt(dataset->GetProjectionRef());
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
      georef.set_transform(transform);

      // Determine the pixel interpretation for the image.  See the
      // comments in GeoReference for more information.  If nothing is
      // encoded in the file, the default is to assume PixelAsArea.
      georef.set_pixel_interpretation(GeoReference::PixelAsArea);
      char **metadata = dataset->GetMetadata();
      if( CSLCount(metadata) > 0 ) {
        for( int i = 0; metadata[i] != NULL; i++ ) {
          std::vector<std::string> split_vec;
          boost::split(split_vec, metadata[i], boost::is_any_of("=") );
          if (split_vec[0] == GDALMD_AREA_OR_POINT && split_vec.size() >= 2)  
            if (boost::trim_copy(split_vec[1]) == GDALMD_AOP_POINT)
              georef.set_pixel_interpretation(GeoReference::PixelAsPoint);
        }
      }
    }

    // Georeference functions need not be invertible.  When we perform a reverse 
    // lookup (e.g. during a geotransformation) we rely on PROJ.4 to pick one 
    // possible value.  However, the georeference might actually place the image 
    // at another (equivalent) location in the projected space.  This is hard to 
    // recover from at run-time, and we currently have no generic solution to this 
    // problem.  In the mean time, we at least test whether the georefernce is 
    // likely to run into this problem (and warn the user) by checking whether 
    // forward- and reverse-projecting the origin pixel lands us back at the origin.
    Vector2 origin = georef.lonlat_to_pixel( georef.pixel_to_lonlat( Vector2() ) );
    if( origin.x()*origin.x() + origin.y()*origin.y() > 0.1 ) {
      vw_out(WarningMessage) << "read_gdal_georeference(): WARNING! Resource file " <<
	resource.filename() << " contains a non-normal georeference.  Bad things "
	"may happen: puppies dying, rainbows fading, mortage foreclosures, etc." << std::endl;
    }
  }

  
  void write_gdal_georeference( DiskImageResourceGDAL& resource, GeoReference const& georef ) {

    boost::shared_ptr<GDALDataset> dataset = resource.get_dataset_ptr();
    if (!dataset) 
      vw_throw( LogicErr() << "GeoReferenceHelperGDAL: Could not write georeference. No file has been opened." );

    // Store the transform matrix
    double geo_transform[6] = { georef.transform()(0,2), georef.transform()(0,0), georef.transform()(0,1), 
                                georef.transform()(1,2), georef.transform()(1,0), georef.transform()(1,1) };
    dataset->SetGeoTransform( geo_transform );

    // This is a little ridiculous, but GDAL can't write geotiffs
    // without a string of Well Known Text (WKT), so we must conjure
    // up an OGRSpatialReference here and use it to convert from
    // proj.4 to WKT.
    // 
    // However, we first set the datum parameters in the
    // ORGSpatialReference object directly.
    OGRSpatialReference gdal_spatial_ref;
    gdal_spatial_ref.importFromProj4(georef.proj4_str().c_str());
    gdal_spatial_ref.SetGeogCS( "Geographic Coordinate System",
                                georef.datum().name().c_str(),
                                georef.datum().spheroid_name().c_str(),
                                georef.datum().semi_major_axis(),
                                georef.datum().inverse_flattening(),
                                georef.datum().meridian_name().c_str(),
                                georef.datum().meridian_offset() );

    char* wkt_str_tmp;
    gdal_spatial_ref.exportToWkt(&wkt_str_tmp);
    std::string wkt_str = wkt_str_tmp;
    OGRFree(wkt_str_tmp);

    dataset->SetProjection( wkt_str.c_str() );

    // Set the pixel interpretation for the image.  See the comments
    // in GeoReference for more information.
    if (georef.pixel_interpretation() == GeoReference::PixelAsArea) 
      dataset->SetMetadataItem(GDALMD_AREA_OR_POINT, GDALMD_AOP_AREA);
    else // PixelAsPoint
      dataset->SetMetadataItem(GDALMD_AREA_OR_POINT, GDALMD_AOP_POINT);
  }

}} // namespace vw::cartography
