// __BEGIN_LICENSE__
//
// Copyright (C) 2006 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration
// (NASA).  All Rights Reserved.
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

/// \file DiskImageResourceGDAL.cc
/// 
/// Provides support for several geospatially referenced file formats
/// via GDAL.
///

#ifdef _MSC_VER
#pragma warning(disable:4244)
#pragma warning(disable:4267)
#pragma warning(disable:4996)
#endif

#include <vw/config.h>

#ifdef VW_HAVE_PKG_GDAL

#include <vw/Cartography/DiskImageResourceGDAL.h>
#include <vw/Cartography/GeoReference.h>

// GDAL Headers
#include "gdal.h"
#include "gdal_priv.h"
#include "cpl_string.h"
#include "ogr_spatialref.h"
#include "ogr_api.h"

#include <vector>
#include <vw/Core/Exception.h>
#include <vw/Image/PixelTypes.h>
#include <boost/algorithm/string.hpp>

namespace vw {
namespace cartography {

  /// \cond INTERNAL
  // Type conversion class for selecting the GDAL data type for
  // various VW channel formats.
  struct vw_channel_id_to_gdal_pix_fmt { 
    static GDALDataType value(ChannelTypeEnum vw_type) {
      switch (vw_type) {
      case VW_CHANNEL_UINT8:   return GDT_Byte;
      case VW_CHANNEL_INT16:   return GDT_Int16;
      case VW_CHANNEL_UINT16:  return GDT_UInt16;
      case VW_CHANNEL_INT32:   return GDT_Int32;
      case VW_CHANNEL_UINT32:  return GDT_UInt32;
      case VW_CHANNEL_FLOAT32: return GDT_Float32;
      case VW_CHANNEL_FLOAT64: return GDT_Float64;
      default: 
        vw_throw( IOErr() << "DiskImageResourceGDAL: Unsupported channel type (" << vw_type << ")." );
        return (GDALDataType)0; // never reached
      }
    }
  };


  // Type conversion class for selecting the GDAL data type for
  // various VW channel formats.
  struct gdal_pix_fmt_to_vw_channel_id {
    static ChannelTypeEnum value(GDALDataType gdal_type) {
      switch (gdal_type) {
      case GDT_Byte:    return VW_CHANNEL_UINT8;
      case GDT_Int16:   return VW_CHANNEL_INT16; 
      case GDT_UInt16:  return VW_CHANNEL_UINT16;
      case GDT_Int32:   return VW_CHANNEL_INT32; 
      case GDT_UInt32:  return VW_CHANNEL_UINT32;
      case GDT_Float32: return VW_CHANNEL_FLOAT32; 
      case GDT_Float64: return VW_CHANNEL_FLOAT64; 
      default:
        vw_throw( IOErr() << "DiskImageResourceGDAL: Unsupported channel type (" << gdal_type << ")." );
        return (ChannelTypeEnum)0; // never reached
      }
    }
  };

  struct gdal_color_interp_for_vw_pixel_format {
    static GDALColorInterp value(int band, int vw_pixel_type) {
      if (vw_pixel_type == VW_PIXEL_GRAY) {
        switch (band) {
        case 0: return GCI_GrayIndex;
        default: vw_throw( LogicErr() << "DiskImageResourceGDAL: requested band does not exist for selected color table." );
        }
      } else if (vw_pixel_type == VW_PIXEL_GRAYA) {
        switch (band) {
        case 0: return GCI_GrayIndex;
        case 1: return GCI_AlphaBand;
        default: vw_throw( LogicErr() << "DiskImageResourceGDAL: requested band does not exist for selected color table." );
        }
      } else if (vw_pixel_type == VW_PIXEL_RGB) {
        switch (band) {
        case 0: return GCI_RedBand;
        case 1: return GCI_GreenBand;
        case 2: return GCI_BlueBand;
        default: vw_throw( LogicErr() << "DiskImageResourceGDAL: requested band does not exist for selected color table." );
        }
      } else if (vw_pixel_type == VW_PIXEL_RGBA) {
        switch (band) {
        case 0: return GCI_RedBand;
        case 1: return GCI_GreenBand;
        case 2: return GCI_BlueBand;
        case 3: return GCI_AlphaBand;
        default: vw_throw( LogicErr() << "DiskImageResourceGDAL: requested band does not exist for selected color table." );
        }
      } else if (vw_pixel_type == VW_PIXEL_HSV) {
        switch (band) {
        case 0: return GCI_HueBand;
        case 1: return GCI_SaturationBand;
        case 2: return GCI_LightnessBand;
        default: vw_throw( LogicErr() << "DiskImageResourceGDAL: requested band does not exist for selected color table." );
        }
      } else {
        return GCI_GrayIndex;
      }
    }
  };

  static std::string file_extension( std::string const& filename ) {
    std::string::size_type dot = filename.find_last_of('.');
    std::string extension = filename.substr( dot );
    boost::to_lower( extension );
    return extension;
  }
  
  // GDAL Support many file formats not specifically enabled here.
  // Consult http://www.remotesensing.org/gdal/formats_list.html for a
  // list of available formats.
  struct gdal_file_format_from_filename {
    static std::string format(std::string const& filename) {

      if (file_extension(filename) == ".tif")      // GeoTiff
        return "GTiff";  

      if (file_extension(filename) == ".grd")      // GMT compatible netcdf
        return "GMT";  

      else if (file_extension(filename) == ".dem") // ENVI Labelled raster
        return  "ENVI" ; 

      else if (file_extension(filename) == ".bil") // ESRI .hdr labelled
        return "EHdr" ;

      else if (file_extension(filename) == ".jpg") // JPEG JFIF
        return "JPEG" ;  

      else if (file_extension(filename) == ".jp2" || 
               file_extension(filename) == ".j2k") // JPEG 2000
        return "JPEG2000"; 

      else if (file_extension(filename) == ".png") // PNG
        return "PNG";     

      else {
        vw_throw( IOErr() << "DiskImageResourceGDAL: \"" << file_extension(filename) << "\" is an unsupported file extension." );
        return std::string(); // never reached
      }
    }
  };
  /// \endcond

  void DiskImageResourceGDAL::write_georeference( GeoReference const& georef ) {

    if (!m_dataset) 
      vw_throw( LogicErr() << "DiskImageResourceGDAL: Could not write georeference. No file has been opened." );
    GDALDataset* dataset = (GDALDataset*)m_dataset;

    // Store the transform matrix
    double geo_transform[6] = { georef.transform()(0,2), georef.transform()(0,0), georef.transform()(0,1), 
                                georef.transform()(1,2), georef.transform()(1,0), georef.transform()(1,1) };
    dataset->SetGeoTransform( geo_transform );
    dataset->SetProjection( georef.wkt_str().c_str() );
    
  }

  void DiskImageResourceGDAL::read_georeference( GeoReference& georef ) {
    if (!m_dataset) 
      vw_throw( LogicErr() << "DiskImageResourceGDAL: Could not read georeference. No file has been opened." );
    GDALDataset* dataset = (GDALDataset*)m_dataset;
    
    if( dataset->GetProjectionRef()  != NULL ) {
      georef.set_wkt_str(dataset->GetProjectionRef());
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
    }
  }  
  
  /// Bind the resource to a file for reading.  Confirm that we can
  /// open the file and that it has a sane pixel format.
  void DiskImageResourceGDAL::open( std::string const& filename )
  {
    // Register all GDAL file readers and writers and open the data set.
    GDALAllRegister();
    GDALDataset *dataset = (GDALDataset *) GDALOpen( filename.c_str(), GA_ReadOnly );
    if( dataset == NULL ) {
      vw_throw( IOErr() << "DiskImageResourceGDAL: Failed to read " << filename << "." );
    }      

    m_filename = filename;
    m_dataset = (void*) dataset;
    m_format.cols = dataset->GetRasterXSize();
    m_format.rows = dataset->GetRasterYSize();
    double geo_transform[6];
   
   // <test code> 
    vw_out(DebugMessage) << "\n\tMetadata description: " << dataset->GetDescription() << std::endl;
    char** metadata = dataset->GetMetadata();
    vw_out(DebugMessage) << "\tCount: " << CSLCount(metadata) << std::endl;
    for (int i = 0; i < CSLCount(metadata); i++) {
      vw_out(DebugMessage) << "\t\t" << CSLGetField(metadata,i) << std::endl;
    }

    vw_out(DebugMessage) << "\tDriver: " << 
      dataset->GetDriver()->GetDescription() <<
      dataset->GetDriver()->GetMetadataItem( GDAL_DMD_LONGNAME ) << std::endl;
   
    vw_out(DebugMessage) << "\tSize is " <<
      dataset->GetRasterXSize() << "x" <<
      dataset->GetRasterYSize() << "x" <<
      dataset->GetRasterCount() << std::endl;
      
   if( dataset->GetGeoTransform( geo_transform ) == CE_None ) {
     m_geo_transform.set_identity();
     m_geo_transform(0,0) = geo_transform[1];
     m_geo_transform(0,1) = geo_transform[2];
     m_geo_transform(0,2) = geo_transform[0];
     m_geo_transform(1,0) = geo_transform[4];
     m_geo_transform(1,1) = geo_transform[5];
     m_geo_transform(1,2) = geo_transform[3];
     vw_out(DebugMessage) << "\tAffine transform: " << m_geo_transform << std::endl;
   }
   // <test code> 
   
   // We do our best here to determine what pixel format the GDAL image is in.  
   if ( dataset->GetRasterCount() == 1 && 
        std::string(GDALGetColorInterpretationName(dataset->GetRasterBand(1)->GetColorInterpretation())) == std::string("Gray")) {
     m_format.pixel_format = VW_PIXEL_GRAY;     
     m_format.planes = 1;
   } else if ( dataset->GetRasterCount() == 2 && 
               std::string(GDALGetColorInterpretationName(dataset->GetRasterBand(1)->GetColorInterpretation())) == std::string("Gray") &&
               std::string(GDALGetColorInterpretationName(dataset->GetRasterBand(2)->GetColorInterpretation())) == std::string("Alpha")) {
     m_format.pixel_format = VW_PIXEL_GRAYA;
     m_format.planes = 1;
   } else if ( dataset->GetRasterCount() == 3 && 
               std::string(GDALGetColorInterpretationName(dataset->GetRasterBand(1)->GetColorInterpretation())) == std::string("Red") &&
               std::string(GDALGetColorInterpretationName(dataset->GetRasterBand(2)->GetColorInterpretation())) == std::string("Green") &&
               std::string(GDALGetColorInterpretationName(dataset->GetRasterBand(3)->GetColorInterpretation())) == std::string("Blue")) {
     m_format.pixel_format = VW_PIXEL_RGB;
     m_format.planes = 1;
   } else if ( dataset->GetRasterCount() == 4 && 
               std::string(GDALGetColorInterpretationName(dataset->GetRasterBand(1)->GetColorInterpretation())) == std::string("Red") &&
               std::string(GDALGetColorInterpretationName(dataset->GetRasterBand(2)->GetColorInterpretation())) == std::string("Green") &&
               std::string(GDALGetColorInterpretationName(dataset->GetRasterBand(3)->GetColorInterpretation())) == std::string("Blue") &&
               std::string(GDALGetColorInterpretationName(dataset->GetRasterBand(4)->GetColorInterpretation())) == std::string("Alpha")) {
     m_format.pixel_format = VW_PIXEL_RGBA;
     m_format.planes = 1;
   } else {
     m_format.pixel_format = VW_PIXEL_SCALAR;    
     m_format.planes = dataset->GetRasterCount();
   }
   m_format.channel_type = gdal_pix_fmt_to_vw_channel_id::value(dataset->GetRasterBand(1)->GetRasterDataType());
  }

  /// Bind the resource to a file for writing.  
  void DiskImageResourceGDAL::create( std::string const& filename, 
                                      GenericImageFormat const& format )
  {
    VW_ASSERT(format.planes == 1 || format.pixel_format==VW_PIXEL_SCALAR,
              NoImplErr() << "DiskImageResourceGDAL: Cannot create " << filename << "\n\t"
              << "The image cannot have both multiple channels and multiple planes.\n");

    // Store away relevent information into the internal data
    // structure for this DiskImageResource
    m_filename = filename;
    m_format = format;
    int num_bands = std::max( format.planes, num_channels( format.pixel_format ) );

    // Register the various file types with GDAL 
    GDALAllRegister();  
    
    // Open the appropriate GDAL I/O driver, depending on the fileFormat
    // argument specified by the user.
    std::string gdal_format_string = gdal_file_format_from_filename::format(filename);
    vw_out(DebugMessage) << "Creating a new file with the following type: " << gdal_format_string << "   ";
    GDALDriver *driver = GetGDALDriverManager()->GetDriverByName(gdal_format_string.c_str());  
    if( driver == NULL )
      vw_throw( vw::IOErr() << "Error opening selected GDAL file I/O driver." );
    
    char** metadata = driver->GetMetadata();
    if( !CSLFetchBoolean( metadata, GDAL_DCAP_CREATE, FALSE ) )
      vw_throw( vw::IOErr() << "Selected GDAL driver not supported." );
    
    char **options = NULL;
    GDALDataType gdal_pix_fmt = vw_channel_id_to_gdal_pix_fmt::value(format.channel_type);
    GDALDataset *dataset = driver->Create( filename.c_str(), cols(), rows(), num_bands, gdal_pix_fmt, options );

    m_dataset = (void*) dataset;
  }


  /// Read the disk image into the given buffer.
  void DiskImageResourceGDAL::read_generic( GenericImageBuffer const& dest ) const
  {
    VW_ASSERT( dest.format.cols==cols() && dest.format.rows==rows(),
               IOErr() << "Buffer has wrong dimensions in GDAL read." );

    VW_ASSERT( channels() == 1 || planes()==1,
               LogicErr() << "DiskImageResourceGDAL: cannot read an image that has both multiple channels and multiple planes." );
 
    if (!m_dataset) 
      vw_throw( LogicErr() << "DiskImageResourceGDAL: Could not read file. No file has been opened." );

    GDALRasterBand  *band;
    int             blocksize_x, blocksize_y;
    int             bGotMin, bGotMax;
    double          adfMinMax[2];
    
    uint8 *data = new uint8[channel_size(channel_type()) * rows() * cols() * planes() * channels()];
    GenericImageBuffer src;
    src.data = data;
    src.format = m_format;
    src.cstride = channels() * channel_size(src.format.channel_type);
    src.rstride = cols() * channels() * channel_size(src.format.channel_type);
    src.pstride = rows() * cols() * channels() * channel_size(src.format.channel_type);

    int b = 1;
    for ( unsigned p = 0; p < planes(); ++p ) {
      for ( unsigned c = 0; c < channels(); ++c ) {
        band = ((GDALDataset *)m_dataset)->GetRasterBand(b++);

        // <Test code>
        vw_out(DebugMessage) << "\n\tMetadata description: " << band->GetDescription() << std::endl;
        char** metadata = band->GetMetadata();
        vw_out(DebugMessage) << "\tCount: " << CSLCount(metadata) << std::endl;
        for (int i = 0; i < CSLCount(metadata); i++) {
          vw_out(DebugMessage) << "\t\t" << CSLGetField(metadata,i) << std::endl;
        }
        
        band->GetBlockSize( &blocksize_x, &blocksize_y );
        vw_out(DebugMessage)
          << "\tBlock=" << blocksize_x << "x" << blocksize_y
          << " Type=" << GDALGetDataTypeName(band->GetRasterDataType())
          << " ColorInterp=" << GDALGetColorInterpretationName(band->GetColorInterpretation())
          << std::endl;
        
        adfMinMax[0] = band->GetMinimum( &bGotMin );
        adfMinMax[1] = band->GetMaximum( &bGotMax );
        if( ! (bGotMin && bGotMax) )
          GDALComputeRasterMinMax((GDALRasterBandH)band, TRUE, adfMinMax);   
        vw_out(DebugMessage)
          << "\tMin=" << adfMinMax[0] << ", Max=" << adfMinMax[1] << std::endl;
        
        //       if( band->GetOverviewCount() > 0 )
        //         printf( "\tBand has %d overviews.\n", band->GetOverviewCount() );
        
        //       if( band->GetColorTable() != NULL )
        //         printf( "\tBand has a color table with %d entries.\n", 
        //                 band->GetColorTable()->GetColorEntryCount() );
        // <Test code>
        
        // Read in the data one scanline at a time and copy that data into an ImageView
        GDALDataType gdal_pix_fmt = vw_channel_id_to_gdal_pix_fmt::value(channel_type());
        band->RasterIO( GF_Read, 0, 0, m_format.cols, m_format.rows, 
                        (uint8*)src(0,0,p) + channel_size(channel_type())*c, m_format.cols, m_format.rows, gdal_pix_fmt, 
                        channel_size(channel_type())*channels(), channel_size(channel_type())*channels()*cols() );
      }
    }
    convert( dest, src );
    delete [] data;
  }


  // Write the given buffer into the disk image.
  void DiskImageResourceGDAL::write_generic( GenericImageBuffer const& src )
  {
    VW_ASSERT( src.format.cols==cols() && src.format.rows==rows(),
               IOErr() << "Buffer has wrong dimensions in GDAL write." );

    // This is pretty simple since we always write 8-bit integer files.
    // Note that we handle multi-channel images with interleaved planes. 
    // We've already ensured that either planes==1 or channels==1.
    uint8 *data = new uint8[channel_size(channel_type()) * rows() * cols() * planes() * channels()];
    GenericImageBuffer dst;
    dst.data = data;
    dst.format = m_format;
    dst.cstride = channels() * channel_size(dst.format.channel_type);
    dst.rstride = cols() * channels() * channel_size(dst.format.channel_type);
    dst.pstride = rows() * cols() * channels() * channel_size(dst.format.channel_type);
    convert( dst, src );

    // Write the data to the selected raster band. 
    int num_bands = std::max(planes(), num_channels(pixel_format()));
    vw_out(DebugMessage)
      << "\tWriting geo-referenced file " << m_filename.c_str()
      << " (" << cols() << " x " << rows() << ") with " << num_bands << " band(s)." << std::endl;
    
    int b = 1;
    for (unsigned int p = 0; p < dst.format.planes; p++) {
      for (unsigned int c = 0; c < num_channels(dst.format.pixel_format); c++) {
        GDALRasterBand *band = ((GDALDataset*)m_dataset)->GetRasterBand(b++);
	//        band->SetColorInterpretation(gdal_color_interp_for_vw_pixel_format::value(p, dst.format.pixel_format));
        GDALDataType gdal_pix_fmt = vw_channel_id_to_gdal_pix_fmt::value(channel_type());
        band->RasterIO( GF_Write, 0, 0, dst.format.cols, dst.format.rows, 
                        (uint8*)dst(0,0,p) + channel_size(dst.format.channel_type)*c, 
                        dst.format.cols, dst.format.rows, gdal_pix_fmt, 
                        channel_size(dst.format.channel_type)*channels(), channel_size(dst.format.channel_type)*channels()*cols() );
      }  
    }
    
    delete [] data;
  }

  void DiskImageResourceGDAL::flush() { 
    if (m_dataset) 
      delete (GDALDataset*)m_dataset;
  }

  // A FileIO hook to open a file for reading
  vw::DiskImageResource* DiskImageResourceGDAL::construct_open( std::string const& filename ) {
    return new DiskImageResourceGDAL( filename );
  }

  // A FileIO hook to open a file for writing
  vw::DiskImageResource* DiskImageResourceGDAL::construct_create( std::string const& filename,
                                                                  GenericImageFormat const& format ) {
    return new DiskImageResourceGDAL( filename, format );
  }
  

}} // namespace vw::cartography

#endif // HAVE_PKG_GDAL
