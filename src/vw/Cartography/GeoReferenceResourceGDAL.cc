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
    GDALDataset *dataset = static_cast<GDALDataset*>(resource.get_read_dataset_ptr());    
    if (!dataset) 
      vw_throw( LogicErr() << "read_gdal_georeference: Could not read georeference. No file has been opened." );
    
    // Pull the projection and datum information out of the file if available
    if( dataset->GetProjectionRef() != NULL ) {

      // Create a GDAL spatial ref object
      OGRSpatialReference gdal_spatial_ref;
      const char* wkt = dataset->GetProjectionRef();
      char** wkt_ptr = (char**)(&wkt);
      gdal_spatial_ref.importFromWkt(wkt_ptr);

      
      // Read projection information out of the file
      char* proj_str_tmp;
      gdal_spatial_ref.exportToProj4(&proj_str_tmp);
      std::string proj4_str = proj_str_tmp;
      CPLFree( proj_str_tmp );
      // For debugging:
      //      vw_out(0) << "PROJ in --> " << proj4_str << "\n";

      std::vector<std::string> input_strings;
      std::vector<std::string> output_strings;
      std::vector<std::string> datum_strings;
      std::string trimmed_proj4_str = boost::trim_copy(proj4_str);
      boost::split( input_strings, trimmed_proj4_str, boost::is_any_of(" ") ); 
      for (unsigned int i = 0; i < input_strings.size(); ++i) {

        // Pick out the parts of the projection string that pertain to
        // map projections.  We essentially want to eliminate all of
        // the strings that have to do with the datum, since those are
        // handled by interacting directly with the
        // OGRSpatialReference below. This is sort of messy, but it's
        // the easiest way to do this, as far as I can tell.
        if ((input_strings[i].find("+proj=") == 0) || 
            (input_strings[i].find("+x_0=") == 0) || 
            (input_strings[i].find("+y_0=") == 0) ||
            (input_strings[i].find("+lon") == 0) || 
            (input_strings[i].find("+lat") == 0) || 
            (input_strings[i].find("+k=") == 0) || 
            (input_strings[i].find("+lat_ts=") == 0) || 
            (input_strings[i].find("+ns") == 0) || 
            (input_strings[i].find("+no_cut") == 0) || 
            (input_strings[i].find("+h=") == 0) || 
            (input_strings[i].find("+W=") == 0) || 
            (input_strings[i].find("+units=") == 0) ||
            (input_strings[i].find("+zone=") == 0)) {
          output_strings.push_back(input_strings[i]);
        } else if ((input_strings[i].find("+ellps=") == 0) ||
                   (input_strings[i].find("+datum=") == 0)) {
          // We put these in the proj4_str for the Datum class.
          datum_strings.push_back(input_strings[i]);
        }
      }
      std::ostringstream strm;
      for (unsigned int i = 0; i < output_strings.size(); ++i) {
        strm << output_strings[i] << " ";
      }
      // For debugging:
      //      vw_out(0) << "     out --> " << strm.str() << "\n";

      // If the file contains no projection related information, we
      // supply proj.4 with a "default" interpretation that the file
      // is in geographic (unprojected) coordinates.
      if (output_strings.size() == 0) 
        georef.set_proj4_projection_str("+proj=longlat");
      else
        georef.set_proj4_projection_str(strm.str());
      
      // Read in the datum information
      Datum datum;
      const char* datum_name = gdal_spatial_ref.GetAttrValue("DATUM");
      if (datum_name) { datum.name() = datum_name; }
      const char* spheroid_name = gdal_spatial_ref.GetAttrValue("SPHEROID");
      if (spheroid_name) { datum.spheroid_name() = spheroid_name; }
      const char* meridian_name = gdal_spatial_ref.GetAttrValue("PRIMEM");
      if (meridian_name) { datum.meridian_name() = meridian_name; }
      OGRErr e1, e2;
      double semi_major = gdal_spatial_ref.GetSemiMajor(&e1);
      double semi_minor = gdal_spatial_ref.GetSemiMinor(&e2);
      if (e1 != OGRERR_FAILURE && e2 != OGRERR_FAILURE) { 
        datum.set_semi_major_axis(semi_major);
        datum.set_semi_minor_axis(semi_minor);
      }
      datum.meridian_offset() = gdal_spatial_ref.GetPrimeMeridian();
      // Set the proj4 string for datum.
      std::stringstream datum_proj4_ss;
      for(unsigned i=0; i < datum_strings.size(); i++)
          datum_proj4_ss << datum_strings[i] << ' ';
      // Add the current proj4 string in the case that our ellipse/datum 
      // values are empty.
      if(boost::trim_copy(datum_proj4_ss.str()) == "")
        datum_proj4_ss << datum.proj4_str();
      datum.proj4_str() = boost::trim_copy(datum_proj4_ss.str());
      georef.set_datum(datum);
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
  }

  
  void write_gdal_georeference( DiskImageResourceGDAL& resource, GeoReference const& georef ) {

    GDALDataset *dataset = static_cast<GDALDataset*>(resource.get_write_dataset_ptr());
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
