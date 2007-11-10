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
#include <iostream>
#ifdef VW_HAVE_PKG_GDAL

#include <vw/FileIO/DiskImageResourceGDAL.h>
#include <vw/FileIO/JP2.h>

// GDAL Headers
#include "gdal.h"
#include "gdal_priv.h"
#include "cpl_string.h"
#include "ogr_spatialref.h"
#include "ogr_api.h"

#include <vector>
#include <list>
#include <vw/Core/Exception.h>
#include <vw/Image/PixelTypes.h>
#include <boost/algorithm/string.hpp>

namespace {

  bool has_gmljp2(vw::JP2File* f)
  {
    vw::JP2Box* b;
    vw::JP2Box* b2;
    vw::JP2SuperBox::JP2BoxIterator pos;
    bool retval = false;

    while((b = f->find_box(0x61736F63 /*"asoc"*/, &pos)))
    {
      // GMLJP2 outer Association box must contain a Label box (with label "gml.data")
      // as its first sub-box
      b2 = ((vw::JP2SuperBox*)b)->find_box(0);
      if(b2 && b2->box_type() == 0x6C626C20 /*"lbl\040"*/ && strncmp((char*)(((vw::JP2DataBox*)b2)->data()), "gml.data", ((vw::JP2DataBox*)b2)->bytes_dbox()) == 0)
      {
        retval = true;
        break;
      }
    }

    return retval;
  }

}

namespace vw {

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

  static std::string file_extension( std::string const& filename ) {
    std::string::size_type dot = filename.find_last_of('.');
    std::string extension = filename.substr( dot );
    boost::to_lower( extension );
    return extension;
  }

  static bool is_jp2( std::string const& filename ) {
    std::string extension = file_extension( filename );
    return (extension == ".jp2" || extension == ".j2k");
  }
  
  // GDAL Support many file formats not specifically enabled here.
  // Consult http://www.remotesensing.org/gdal/formats_list.html for a
  // list of available formats.
  struct gdal_file_format_from_filename {
    static std::list<std::string> format(std::string const& filename) {
      std::list<std::string> retval;

      if (file_extension(filename) == ".tif" ||
          file_extension(filename) == ".tiff") {     // GeoTiff
        retval.push_back("GTiff");
        return retval;
      }  

      if (file_extension(filename) == ".grd") {      // GMT compatible netcdf
        retval.push_back("GMT");
        return retval;  
      }

      else if (file_extension(filename) == ".dem") { // ENVI Labelled raster
        retval.push_back("ENVI");
        return  retval ; 
      }

      else if (file_extension(filename) == ".bil") { // ESRI .hdr labelled
        retval.push_back("EHdr");
        return retval ;
      }

      else if (file_extension(filename) == ".jpg" ||
               file_extension(filename) == ".jpeg") { // JPEG JFIF
        retval.push_back("JPEG");
        return retval ;  
      }

      else if (is_jp2(filename)) {                    // JPEG 2000
        retval.push_back("JP2ECW");
        retval.push_back("JPEG2000");
        return retval; 
      }

      else if (file_extension(filename) == ".png") { // PNG
        retval.push_back("PNG");
        return retval;     
      }

      else if (file_extension(filename) == ".img") { // Imagine
        retval.push_back("HFA");
        return retval;     
      }

      else {
        vw_throw( IOErr() << "DiskImageResourceGDAL: \"" << file_extension(filename) << "\" is an unsupported file extension." );
        return retval; // never reached
      }
    }
  };

  static void convert_jp2(std::string const& filename) {
    uint8* d_original = 0;
    uint8* d_converted = 0;
    FILE* fp;
    uint64 nbytes_read, nbytes_converted, nbytes_written;
    uint64 i;
    int c;
    bool has_gml;
    int retval;
  
    if (!(fp = fopen(filename.c_str(), "r")))
      vw_throw( IOErr() << "convert_jp2: Failed to read " << filename << "." );
  
    for (nbytes_read = 0; fgetc(fp) != EOF; nbytes_read++);
    rewind(fp);
  
    d_original = new uint8[nbytes_read];
  
    for (i = 0; i < nbytes_read && (c = fgetc(fp)) != EOF; d_original[i] = (uint8)c, i++);

    if (!(i == nbytes_read  && fgetc(fp) == EOF))
      vw_throw( IOErr() << "convert_jp2: Size of " << filename << " has changed." );

    fclose(fp);
  
    JP2File f(d_original, nbytes_read);
    //f.print();
    has_gml = has_gmljp2(&f);
    retval = f.convert_to_jpx();
    if (retval != 0)
      vw_throw( IOErr() << "convert_jp2: Failed to convert " << filename << " to jp2-compatible jpx." );
    if (has_gml) {
      JP2ReaderRequirementsList req;
      // 67 is the (standard) requirement number for GMLJP2
      req.push_back(std::make_pair(67, false));
      retval = f.add_requirements(req);
      if (retval != 0)
        vw_throw( IOErr() << "convert_jp2: Failed to add GMLJP2 requirement to jp2-compatible jpx " << filename << "." );
    }
    //f.print();

    nbytes_converted = f.bytes();
    d_converted = new uint8[nbytes_converted];
    f.serialize(d_converted);

    if (!(fp = fopen(filename.c_str(), "w")))
      vw_throw( IOErr() << "convert_jp2: Failed to open " << filename << " for writing." );

    nbytes_written = fwrite(d_converted, 1, nbytes_converted, fp);

    fclose(fp);

    if (nbytes_written != nbytes_converted)
      vw_throw( IOErr() << "convert_jp2: Failed to write " << filename << "." );
  
    if (d_original)
      delete[] d_original;
    if (d_converted)
      delete[] d_converted;
  }
  /// \endcond
  
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
    for( int i=1; i<=dataset->GetRasterCount(); ++i )
    if ( dataset->GetRasterCount() == 1 && 
         dataset->GetRasterBand(1)->GetColorInterpretation() == GCI_GrayIndex ) {
      m_format.pixel_format = VW_PIXEL_GRAY;     
      m_format.planes = 1;
    } else if ( dataset->GetRasterCount() == 2 && 
                dataset->GetRasterBand(1)->GetColorInterpretation() == GCI_GrayIndex /* &&
                dataset->GetRasterBand(2)->GetColorInterpretation() == GCI_AlphaBand */ ) {
      m_format.pixel_format = VW_PIXEL_GRAYA;
      m_format.planes = 1;
    } else if ( dataset->GetRasterCount() == 3 && 
                dataset->GetRasterBand(1)->GetColorInterpretation() == GCI_RedBand &&
                dataset->GetRasterBand(2)->GetColorInterpretation() == GCI_GreenBand &&
                dataset->GetRasterBand(3)->GetColorInterpretation() == GCI_BlueBand) {
      m_format.pixel_format = VW_PIXEL_RGB;
      m_format.planes = 1;
    } else if ( dataset->GetRasterCount() == 4 && 
                dataset->GetRasterBand(1)->GetColorInterpretation() == GCI_RedBand &&
                dataset->GetRasterBand(2)->GetColorInterpretation() == GCI_GreenBand &&
                dataset->GetRasterBand(3)->GetColorInterpretation() == GCI_BlueBand &&
                dataset->GetRasterBand(4)->GetColorInterpretation() == GCI_AlphaBand) {
      m_format.pixel_format = VW_PIXEL_RGBA;
     m_format.planes = 1;
    } else {
      m_format.pixel_format = VW_PIXEL_SCALAR;    
      m_format.planes = dataset->GetRasterCount();
    }
    m_format.channel_type = gdal_pix_fmt_to_vw_channel_id::value(dataset->GetRasterBand(1)->GetRasterDataType());
    // Special limited read-only support for palette-based images
    if( dataset->GetRasterCount() == 1 && 
        dataset->GetRasterBand(1)->GetColorInterpretation() == GCI_PaletteIndex &&
        m_format.channel_type == VW_CHANNEL_UINT8 ) {
      m_format.pixel_format = VW_PIXEL_RGBA;
      m_format.planes = 1;
      GDALColorTable *color_table = dataset->GetRasterBand(1)->GetColorTable();
      int num_entries = color_table->GetColorEntryCount();
      m_palette.resize( num_entries );
      GDALColorEntry color;
      for( int i=0; i<num_entries; ++i ) {
        color_table->GetColorEntryAsRGB( i, &color );
        m_palette[i] = PixelRGBA<uint8>( color.c1, color.c2, color.c3, color.c4 );
      }
    }
  }

  /// Bind the resource to a file for writing.  
  void DiskImageResourceGDAL::create( std::string const& filename, 
                                      ImageFormat const& format )
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
    std::list<std::string> gdal_format_string = gdal_file_format_from_filename::format(filename);
    GDALDriver *driver = NULL;
    std::list<std::string>::iterator i;
    bool unsupported_driver = false;
    for( i = gdal_format_string.begin(); i != gdal_format_string.end() && driver == NULL; i++ ) {
      vw_out(DebugMessage) << "Attempting to creating a new file with the following type: " << (*i) << std::endl;
      driver = GetGDALDriverManager()->GetDriverByName((*i).c_str());
      if( driver == NULL )
        continue;
        
      char** metadata = driver->GetMetadata();
      if( !CSLFetchBoolean( metadata, GDAL_DCAP_CREATE, FALSE ) ) {
        vw_out(DebugMessage) << "GDAL driver " << (*i) << " does not support create." << std::endl;
        driver = NULL;
        unsupported_driver = true;
      }
    }
    if( driver == NULL ) {
      if( unsupported_driver )
        vw_throw( vw::IOErr() << "Selected GDAL driver not supported." );
      else
        vw_throw( vw::IOErr() << "Error opening selected GDAL file I/O driver." );
    }
    
    char **options = NULL;
    options = CSLSetNameValue( options, "ALPHA", "YES" );
    options = CSLSetNameValue( options, "INTERLEAVE", "PIXEL" );
    if( m_format.pixel_format == VW_PIXEL_RGB || m_format.pixel_format == VW_PIXEL_RGBA ) {
      options = CSLSetNameValue( options, "PHOTOMETRIC", "RGB" );
    }
    GDALDataType gdal_pix_fmt = vw_channel_id_to_gdal_pix_fmt::value(format.channel_type);
    GDALDataset *dataset = driver->Create( filename.c_str(), cols(), rows(), num_bands, gdal_pix_fmt, options );
    CSLDestroy( options );

    m_dataset = (void*) dataset;
  }


  /// Read the disk image into the given buffer.
  void DiskImageResourceGDAL::read( ImageBuffer const& dest, BBox2i const& bbox ) const
  {
    VW_ASSERT( channels() == 1 || planes()==1,
               LogicErr() << "DiskImageResourceGDAL: cannot read an image that has both multiple channels and multiple planes." );
 
    if (!m_dataset) 
      vw_throw( LogicErr() << "DiskImageResourceGDAL: Could not read file. No file has been opened." );

    int             blocksize_x, blocksize_y;
    int             bGotMin, bGotMax;
    double          adfMinMax[2];
    
    uint8 *data = new uint8[channel_size(channel_type()) * bbox.width() * bbox.height() * planes() * channels()];
    ImageBuffer src;
    src.data = data;
    src.format = m_format;
    src.format.cols = bbox.width();
    src.format.rows = bbox.height();
    src.cstride = channels() * channel_size(src.format.channel_type);
    src.rstride = src.cstride * src.format.cols;
    src.pstride = src.rstride * src.format.rows;

    if( m_palette.size() == 0 ) {
      for ( int32 p = 0; p < planes(); ++p ) {
        for ( int32 c = 0; c < channels(); ++c ) {
          GDALRasterBand  *band = ((GDALDataset *)m_dataset)->GetRasterBand(c+1);
          GDALDataType gdal_pix_fmt = vw_channel_id_to_gdal_pix_fmt::value(channel_type());
          band->RasterIO( GF_Read, bbox.min().x(), bbox.min().y(), bbox.width(), bbox.height(),
                          (uint8*)src(0,0,p) + channel_size(src.format.channel_type)*c, 
                          src.format.cols, src.format.rows, gdal_pix_fmt, src.cstride, src.rstride );
        }
      }
    }
    else { // palette conversion
      GDALRasterBand  *band = ((GDALDataset *)m_dataset)->GetRasterBand(1);
      uint8 *index_data = new uint8[bbox.width() * bbox.height()];
      band->RasterIO( GF_Read, bbox.min().x(), bbox.min().y(), bbox.width(), bbox.height(),
                      index_data, bbox.width(), bbox.height(), GDT_Byte, 1, bbox.width() );
      PixelRGBA<uint8> *rgba_data = (PixelRGBA<uint8>*) data;
      for( int i=0; i<bbox.width()*bbox.height(); ++i )
        rgba_data[i] = m_palette[index_data[i]];
      delete [] index_data;
    }
    convert( dest, src );
    delete [] data;
  }


  // Write the given buffer into the disk image.
  void DiskImageResourceGDAL::write( ImageBuffer const& src, BBox2i const& bbox )
  {
    // This is pretty simple since we always write 8-bit integer files.
    // Note that we handle multi-channel images with interleaved planes. 
    // We've already ensured that either planes==1 or channels==1.
    uint8 *data = new uint8[channel_size(channel_type()) * bbox.width() * bbox.height() * planes() * channels()];
    ImageBuffer dst;
    dst.data = data;
    dst.format = m_format;
    dst.format.cols = bbox.width();
    dst.format.rows = bbox.height();
    dst.cstride = channels() * channel_size(dst.format.channel_type);
    dst.rstride = dst.format.cols * dst.cstride;
    dst.pstride = dst.format.rows * dst.rstride;
    convert( dst, src );

    GDALDataType gdal_pix_fmt = vw_channel_id_to_gdal_pix_fmt::value(channel_type());
    for (int32 p = 0; p < dst.format.planes; p++) {
      for (int32 c = 0; c < num_channels(dst.format.pixel_format); c++) {
        GDALRasterBand *band = ((GDALDataset*)m_dataset)->GetRasterBand(c+1);
        band->RasterIO( GF_Write, bbox.min().x(), bbox.min().y(), bbox.width(), bbox.height(),
                        (uint8*)dst(0,0,p) + channel_size(dst.format.channel_type)*c, 
                        dst.format.cols, dst.format.rows, gdal_pix_fmt, dst.cstride, dst.rstride );
      }  
    }
    
    //FIXME: if we allow partial writes, m_convert_jp2 should probably only be set on complete writes
    if (is_jp2(m_filename))
      m_convert_jp2 = true;
    
    delete [] data;
  }

  Vector2i DiskImageResourceGDAL::native_block_size() const {
    // GDAL assumes a single-row stripsize even for file formats 
    // like PNG for which it does not support true strip access.
    // Thus, we check the file driver type before accepting GDAL's 
    // block size as our own.
    GDALDataset *dataset = (GDALDataset*)m_dataset;
    if ( dataset->GetDriver() != GetGDALDriverManager()->GetDriverByName("GTiff") ) {
      return Vector2i(cols(),rows());
    }
    else {
      GDALRasterBand *band = dataset->GetRasterBand(1);
      int xsize, ysize;
      band->GetBlockSize(&xsize,&ysize);
      return Vector2i(xsize,ysize);
    }
  }

  void DiskImageResourceGDAL::flush() { 
    if (m_dataset) {
      delete (GDALDataset*)m_dataset;
      if (m_convert_jp2)
        convert_jp2(m_filename);
    }
  }

  // A FileIO hook to open a file for reading
  vw::DiskImageResource* DiskImageResourceGDAL::construct_open( std::string const& filename ) {
    return new DiskImageResourceGDAL( filename );
  }

  // A FileIO hook to open a file for writing
  vw::DiskImageResource* DiskImageResourceGDAL::construct_create( std::string const& filename,
                                                                  ImageFormat const& format ) {
    return new DiskImageResourceGDAL( filename, format );
  }
  

} // namespace vw

#endif // HAVE_PKG_GDAL
