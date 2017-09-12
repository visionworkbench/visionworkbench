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

  bool read_gdal_georeference( GeoReference& georef,
                               DiskImageResourceGDAL const& resource ) {
    boost::shared_ptr<GDALDataset>dataset = resource.get_dataset_ptr();
    if (!dataset)
      vw_throw( LogicErr() << "read_gdal_georeference: Could not read georeference. No file has been opened." );

    // Pull the projection and datum information out of the file if available
    if( dataset->GetProjectionRef() != NULL )
      georef.set_wkt(dataset->GetProjectionRef());

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
    } else {
      return false;
    }

    // Update the longitude center of the georef using the image size.
    int cols = resource.format().cols;
    int rows = resource.format().rows; 
    BBox2 image_bbox(0,0,cols, rows);
    georef.update_lon_center(image_bbox);
    

    // Georeference functions need not be invertible.  When we perform
    // a reverse lookup (e.g. during a geotransformation) we rely on
    // PROJ.4 to pick one possible value.  However, the georeference
    // might actually place the image at another (equivalent) location
    // in the projected space.  This is hard to recover from at
    // run-time, and we currently have no generic solution to this
    // problem.  In the mean time, we at least test whether the
    // georeference is likely to run into this problem (and warn the
    // user) by checking whether forward- and reverse-projecting the
    // four corner pixels lands us back at the same pixel.
    
    std::vector<Vector2> test_pixels(4);
    test_pixels[0] = Vector2(0,     0);
    test_pixels[1] = Vector2(0,     rows-1);
    test_pixels[2] = Vector2(cols-1,0);
    test_pixels[3] = Vector2(cols-1,rows-1);
    
    for (int i=0; i<4; ++i) {
      double error = 9999;
      bool have_error = false;
      try {
        error = georef.test_pixel_reprojection_error(test_pixels[i]);
      } catch (std::exception &e) {
        vw_out(WarningMessage) << e.what() << std::endl;
        have_error = true;
      }
      if ( error > 0.1 || have_error ) {
        vw_out(WarningMessage) << "read_gdal_georeference(): WARNING! Resource file " <<
          resource.filename() << " contains a non-normal georeference." << std::endl;
      }
    }
    
    return true;
  }

  void write_gdal_georeference( DiskImageResourceGDAL& resource,
                                GeoReference const& georef ) {

    boost::shared_ptr<GDALDataset> dataset = resource.get_dataset_ptr();
    if (!dataset)
      vw_throw( LogicErr() << "GeoReferenceHelperGDAL: Could not write georeference. No file has been opened." );

    // Store the transform matrix
    double geo_transform[6] = { georef.transform()(0,2), georef.transform()(0,0),
                                georef.transform()(0,1), georef.transform()(1,2),
                                georef.transform()(1,0), georef.transform()(1,1) };
    dataset->SetGeoTransform( geo_transform );

    // This is a little ridiculous, but GDAL can't write geotiffs
    // without a string of Well Known Text (WKT), so we must conjure
    // up an OGRSpatialReference here and use it to convert from proj.4 to WKT.
    std::string wkt_str = georef.get_wkt();
    dataset->SetProjection( wkt_str.c_str() );

    // Set the pixel interpretation for the image.  See the comments
    // in GeoReference for more information.
    if (georef.pixel_interpretation() == GeoReference::PixelAsArea)
      dataset->SetMetadataItem(GDALMD_AREA_OR_POINT, GDALMD_AOP_AREA);
    else // PixelAsPoint
      dataset->SetMetadataItem(GDALMD_AREA_OR_POINT, GDALMD_AOP_POINT);
  }

  // Read an arbitrary name = value pair from the geoheader.
  bool read_gdal_string( DiskImageResourceGDAL const& resource,
                         std::string const& str_name,
                         std::string & str_val ) {

    boost::shared_ptr<GDALDataset>dataset = resource.get_dataset_ptr();
    if (!dataset)
      vw_throw( LogicErr() << "read_gdal_string: Could not read string. "
                << "No file has been opened." );

    char **metadata = dataset->GetMetadata();
    if( CSLCount(metadata) > 0 ) {
      for( int i = 0; metadata[i] != NULL; i++ ) {
        std::vector<std::string> split_vec;
        boost::split(split_vec, metadata[i], boost::is_any_of("=") );
        if (split_vec[0] == str_name && split_vec.size() >= 2){
          str_val = split_vec[1];
          return true;
        }
      }
    }
    return false;
  }

  bool read_gdal_strings( DiskImageResourceGDAL const& resource, 
                          std::map<std::string, std::string>& value_pairs) {

    boost::shared_ptr<GDALDataset>dataset = resource.get_dataset_ptr();
    if (!dataset)
      vw_throw( LogicErr() << "read_gdal_string: Could not read string. "
                << "No file has been opened." );

    char **metadata = dataset->GetMetadata();
    if( CSLCount(metadata) > 0 ) {
      for( int i = 0; metadata[i] != NULL; i++ ) {
        std::vector<std::string> split_vec;
        boost::split(split_vec, metadata[i], boost::is_any_of("=") );
        if (split_vec.size() >= 2){
          value_pairs[split_vec[0]] = split_vec[1];
        }
      }
    }
    return false;
  }

  // Write an arbitrary name = value pair in the geoheader.
  void write_gdal_string( DiskImageResourceGDAL& resource,
                          std::string const& str_name,
                          std::string const& str_val ) {

    boost::shared_ptr<GDALDataset> dataset = resource.get_dataset_ptr();
    if (!dataset)
      vw_throw( LogicErr() << "write_gdal_string: Could not write string. "
                           << "No file has been opened." );

    dataset->SetMetadataItem(str_name.c_str(), str_val.c_str());
  }

}} // namespace vw::cartography
