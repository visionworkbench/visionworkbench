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

/// \file DiskImageResourceHDF.cc
/// 
/// Provides support for the HDF file format.
///
/// Currently this support is very limited.  We only support reading 
/// image data stored in SDS swaths, a la HDF-EOS.  Even that we only 
/// support when nothing too fancy is going on.  Pretty much this has 
/// only been tested for reading, color-correcting, and reprojecting 
/// true-color MODIS imagery.
///

#ifdef _MSC_VER
#pragma warning(disable:4244)
#pragma warning(disable:4267)
#pragma warning(disable:4996)
#endif

#include <mfhdf.h>

#include <vw/Core/Exception.h>
#include <vw/Core/Debugging.h>

#include <vw/FileIO/DiskImageResourceHDF.h>

// *******************************************************************
// The HDF library interface class
// *******************************************************************

namespace vw {
  static ChannelTypeEnum hdf_to_vw_type( ::int32 data_type ) {
    switch( data_type ) {
    case DFNT_CHAR8:
    case DFNT_UCHAR8:
      return VW_CHANNEL_CHAR;
    case DFNT_INT8:
      return VW_CHANNEL_INT8;
    case DFNT_UINT8:
      return VW_CHANNEL_UINT8;
    case DFNT_INT16:
      return VW_CHANNEL_INT16;
    case DFNT_UINT16:
      return VW_CHANNEL_UINT16;
    case DFNT_INT32:
      return VW_CHANNEL_INT32;
    case DFNT_UINT32:
      return VW_CHANNEL_UINT32;
    case DFNT_FLOAT32:
      return VW_CHANNEL_FLOAT32;
    case DFNT_FLOAT64:
      return VW_CHANNEL_FLOAT64;
    default:
      vw_throw( IOErr() << "Unsupported SDS channel type in HDF file (" << data_type << ")!" );
    }
    return VW_CHANNEL_UNKNOWN;
  }
}

class vw::DiskImageResourceInfoHDF {
public:
  struct SDSInfo {
    ::int32 id;
    char name[65];
    ::int32 rank;
    ::int32 dim_sizes[MAX_VAR_DIMS];
    ::int32 data_type;
    ::int32 n_attrs;
  };

  struct PlaneInfo {
    ::int32 sds;
    ::int32 band;
  };

  DiskImageResourceHDF &resource;
  ::int32 sd_id;
  std::vector<SDSInfo> sds_info;
  std::vector<PlaneInfo> plane_info;

  DiskImageResourceInfoHDF( std::string const& filename, DiskImageResourceHDF& resource  ) : resource(resource), sd_id(FAIL) {
    
    if( (sd_id = SDstart( filename.c_str(), DFACC_READ )) == FAIL )
      vw_throw( IOErr() << "Unable to open HDF file \"" << filename << "\"!" );

    ::int32 n_datasets, n_fileattrs;
    if( SDfileinfo( sd_id, &n_datasets, &n_fileattrs ) == FAIL )
      vw_throw( IOErr() << "Unable to read info from HDF file \"" << filename << "\"!" );

    vw_out(VerboseDebugMessage) << "HDF SD file \"" << filename << "\": " << n_datasets << " datasets, " << n_fileattrs << " fileattrs" << std::endl;

    sds_info.resize( n_datasets );
    for( ::int32 i=0; i<n_datasets; ++i ) {

      if( (sds_info[i].id = SDselect( sd_id, i )) == FAIL )
        vw_throw( IOErr() << "Unable to select SDS in HDF file \"" << filename << "\"!" );

      if( SDgetinfo( sds_info[i].id, sds_info[i].name, &sds_info[i].rank, sds_info[i].dim_sizes, &sds_info[i].data_type, &sds_info[i].n_attrs ) == FAIL )
        vw_throw( IOErr() << "Unable to read SDS info from HDF file \"" << filename << "\"!" );

      if( SDiscoordvar( sds_info[i].id ) )
        vw_out(VerboseDebugMessage) << "  DIM ";
      else
        vw_out(VerboseDebugMessage) << "  SDS ";

      vw_out(VerboseDebugMessage) << i << ": \"" << sds_info[i].name << "\" ";
      for( int j=sds_info[i].rank-1; ; ) {
        vw_out(VerboseDebugMessage) << sds_info[i].dim_sizes[j];
        if( j-- == 0 ) break;
        vw_out(VerboseDebugMessage) << "x";
      }
      vw_out(VerboseDebugMessage) << " (data type " << sds_info[i].data_type << ", " << sds_info[i].n_attrs << " attrs)" << std::endl;

      if( ! SDiscoordvar( sds_info[i].id ) ) {
        for( int j=sds_info[i].rank-1; j>=0; --j ) {
          ::int32 dim_id = SDgetdimid( sds_info[i].id, j );
          if( dim_id == FAIL )
            vw_throw( IOErr() << "Unable to get dimension ID in HDF file \"" << filename << "\"!" );
          char dim_name[65];
          ::int32 dim_size;
          ::int32 dim_data_type;
          ::int32 dim_n_attrs;
          if( SDdiminfo( dim_id, dim_name, &dim_size, &dim_data_type, &dim_n_attrs ) == FAIL )
            vw_throw( IOErr() << "Unable to get dimension info in HDF file \"" << filename << "\"!" );
          // FIXME: We really don't support this yet.
          // if( dim_data_type != 0 )
          //   vw_throw( NoImplErr() << "SDS dimension scales not supported (HDF file \"" << filename << "\")!" );
          vw_out(VerboseDebugMessage) << "    Dim " << j << ": \"" << dim_name << "\" (size " << dim_size << ", " << dim_n_attrs << " attrs)" << std::endl;
        }
      }

      SDendaccess( sds_info[i].id );
    }
  }
  
  ~DiskImageResourceInfoHDF() {
    if( sd_id != FAIL ) SDend( sd_id );
  }

  ImageFormat select_sds_planes( std::vector<DiskImageResourceHDF::SDSBand> const& sds_planes ) {
    ::int32 cols=0, rows=0, data_type=0;
    std::vector<PlaneInfo> new_plane_info(sds_planes.size());
    // For each requested plane
    for( unsigned plane=0; plane<sds_planes.size(); ++plane ) {
      // Look for an SDS with the given name
      unsigned sds = 0;
      for( ; sds < sds_info.size(); ++sds ) {
        if( sds_info[sds].name == sds_planes[plane].name ) {
          // Make sure the requested band is valid
          ::int32 rank = sds_info[sds].rank;
          ::int32 band = sds_planes[plane].band;
          if( ! (rank==2 && band==0) && ! (rank==3 && band<sds_info[sds].dim_sizes[0]) )
            vw::vw_throw( vw::IOErr() << "Invalid SDS band requested from HDF file \"" << resource.filename() << "\"!" );
          // Get the plane's properties
          ::int32 plane_cols=0, plane_rows=0;
          ::int32 plane_type = sds_info[sds].data_type;
          if( rank == 2 ) {
            plane_cols = sds_info[sds].dim_sizes[1];
            plane_rows = sds_info[sds].dim_sizes[0];
          }
          else if( rank==3 ) {
            plane_cols = sds_info[sds].dim_sizes[2];
            plane_rows = sds_info[sds].dim_sizes[1];
          }
          else vw_throw( NoImplErr() << "Unsupported SDS rank in HDF file \"" << resource.filename() << "\"!" );
          // Make sure all planes have the same properties
          if( plane == 0 ) {
            cols = plane_cols;
            rows = plane_rows;
            data_type = plane_type;
          }
          else {
            if( plane_cols != cols || plane_rows != rows || plane_type != data_type )
              vw_throw( IOErr() << "Requested planes have mismatched dimensions in HDF file \"" << resource.filename() << "\"!" );
          }
          new_plane_info[plane].sds = sds;
          new_plane_info[plane].band = sds_planes[plane].band;
          break;
        }
      }
      if( sds == sds_info.size() )
        vw_throw( IOErr() << "Requested SDS not found in HDF file \"" << resource.filename() << "\"!" );
    }
    ImageFormat new_format;
    new_format.channel_type = hdf_to_vw_type( data_type );
    new_format.cols = cols;
    new_format.rows = rows;
    new_format.planes = sds_planes.size();
    new_format.pixel_format = VW_PIXEL_SCALAR;
    plane_info = new_plane_info;
    vw_out(VerboseDebugMessage) << "Configured resource: " << new_format.cols << "x" << new_format.rows << "x" << new_format.planes << std::endl;
    return new_format;
  }

  ImageFormat select_sds( std::string const& name ) {
    std::vector<DiskImageResourceHDF::SDSBand> sds_planes;
    // Look for an SDS with the given name
    for( unsigned sds = 0; sds < sds_info.size(); ++sds ) {
      if( sds_info[sds].name == name ) {
        if( sds_info[sds].rank == 2 ) {
          sds_planes.push_back( DiskImageResourceHDF::SDSBand( name ) );
        }
        else if( sds_info[sds].rank == 3 ) {
          for( int i=0; i<sds_info[sds].dim_sizes[0]; ++i )
            sds_planes.push_back( DiskImageResourceHDF::SDSBand( name, i ) );
        }
        else vw::vw_throw( vw::IOErr() << "Invalid SDS rank in HDF file \"" << resource.filename() << "\"!" );
        return select_sds_planes( sds_planes );
      }
    }
    vw::vw_throw( vw::IOErr() << "SDS not found in HDF file \"" << resource.filename() << "\"!" );
    return ImageFormat();
  }

  void read( ImageBuffer const& dstbuf, BBox2i const& bbox ) const {
    ImageBuffer srcbuf;
    boost::scoped_array<uint8> buffer( new uint8[ bbox.width() * bbox.height() * resource.planes() * channel_size( resource.channel_type() ) ] );
    srcbuf.data = buffer.get();
    srcbuf.format.cols = bbox.width();
    srcbuf.format.rows = bbox.height();
    srcbuf.format.planes = resource.planes();
    srcbuf.format.channel_type = resource.channel_type();
    srcbuf.format.pixel_format = VW_PIXEL_SCALAR;
    srcbuf.cstride = channel_size( srcbuf.format.channel_type );
    srcbuf.rstride = srcbuf.cstride * bbox.width();
    srcbuf.pstride = srcbuf.rstride * bbox.height();
    // For each requested plane...
    for( uint32 p=0; p<srcbuf.format.planes; ++p ) {
      // Select the SDS
      ::int32 sds_id = SDselect( sd_id, plane_info[p].sds );
      if( sds_id == FAIL ) vw_throw( IOErr() << "Unable to select SDS in HDF file \"" << resource.filename() << "\"!" );

      if( sds_info[plane_info[p].sds].rank == 2 ) {
        ::int32 start[2] = { bbox.min().y(), bbox.min().x() };
        ::int32 edges[2] = { bbox.height(), bbox.width() };
        if ( SDreaddata( sds_id, start, NULL, edges, (uint8*)srcbuf.data + p*srcbuf.pstride ) == FAIL )
          vw_throw( IOErr() << "Unable to read data from HDF file \"" << resource.filename() << "\"!" );
      }
      else if( sds_info[plane_info[p].sds].rank == 3 ) {
        ::int32 start[3] = { plane_info[p].band,   bbox.min().y(), bbox.min().x() };
        ::int32 edges[3] = { 1, bbox.height(), bbox.width() };
        if ( SDreaddata( sds_id, start, NULL, edges, (uint8*)srcbuf.data + p*srcbuf.pstride ) == FAIL )
          vw_throw( IOErr() << "Unable to read data from HDF file \"" << resource.filename() << "\"!" );
      }
      else vw_throw( IOErr() << "Invalid SDS rank in HDF file \"" << resource.filename() << "\"!" );

      SDendaccess( sds_id );
    }
    convert( dstbuf, srcbuf );
  }

  void get_sds_fillvalue( std::string const& sds_name, float32& result ) const {
    ::int32 sds_index = SDnametoindex( sd_id, sds_name.c_str() );
    if( sds_index == FAIL ) vw_throw( NotFoundErr() << "SDS not found!" );
    if( sds_info[sds_index].data_type != DFNT_FLOAT32 )
      vw_throw( IOErr() << "SDS has incorrect data type (expected " << DFNT_FLOAT32 << ", found " << sds_info[sds_index].data_type << ")!" );
    ::int32 sds_id = SDselect( sd_id, sds_index );
    if( sds_id == FAIL ) vw_throw( IOErr() << "Unable to select SDS in HDF file \"" << resource.filename() << "\"!" );
    if( SDgetfillvalue( sds_id, &result ) == FAIL )
      vw_throw( IOErr() << "Unable to read SDS fill value!" );
    SDendaccess( sds_id );
  }

  std::vector<DiskImageResourceHDF::AttrInfo> get_sds_attrs( std::string const& sds_name ) const {
    ::int32 sds_index = SDnametoindex( sd_id, sds_name.c_str() );
    if( sds_index == FAIL ) vw_throw( NotFoundErr() << "SDS not found!" );
    ::int32 sds_id = SDselect( sd_id, sds_index );
    if( sds_id == FAIL ) vw_throw( IOErr() << "Unable to select SDS in HDF file \"" << resource.filename() << "\"!" );
    std::vector<DiskImageResourceHDF::AttrInfo> attrs;
    for( int i=0; i<sds_info[sds_index].n_attrs; ++i ) {
      char attr_name[MAX_NC_NAME];
      ::int32 data_type, n_values;
      ::int32 status = SDattrinfo( sds_id, i, attr_name, &data_type, &n_values);
      if( status == FAIL ) {
        SDendaccess( sds_id );
        vw_throw( IOErr() << "Unable to get SDS attribute info!" );
      }
      DiskImageResourceHDF::AttrInfo info;
      info.name = attr_name;
      info.type = hdf_to_vw_type( data_type );
      info.size = n_values;
      attrs.push_back( info );
    }
    SDendaccess( sds_id );
    return attrs;
  }

  void get_sds_attr( std::string const& sds_name, std::string const& attr_name, std::vector<uint16>& result ) const {
    ::int32 sds_index = SDnametoindex( sd_id, sds_name.c_str() );
    if( sds_index == FAIL ) vw_throw( NotFoundErr() << "SDS not found!" );
    ::int32 sds_id = SDselect( sd_id, sds_index );
    if( sds_id == FAIL ) vw_throw( IOErr() << "Unable to select SDS in HDF file \"" << resource.filename() << "\"!" );
    ::int32 attr_index = SDfindattr( sds_id, attr_name.c_str() );
    if( attr_index == FAIL ) vw_throw( NotFoundErr() << "SDS attribute not found!" );
    char name[65];
    ::int32 data_type, count;
    ::int32 status = SDattrinfo( sds_id, attr_index, name, &data_type, &count );
    if( status == FAIL ) vw_throw( IOErr() << "Unable to get SDS attribute info in HDF file \"" << resource.filename() << "\"!" );
    if( data_type != DFNT_UINT16 ) vw_throw( NotFoundErr() << "SDS attribute with requested type not found! (expected " << DFNT_UINT16 << ", found " << data_type << ")" );
    result.resize( count );
    status = SDreadattr( sds_id, attr_index, &result[0] );
    if( status == FAIL ) vw_throw( IOErr() << "Unable to read SDS attribute in HDF file \"" << resource.filename() << "\"!" );
    vw_out(VerboseDebugMessage) << "Attribute \"" << name << "\": data type " << data_type << ", length " << count << std::endl;
    SDendaccess( sds_id );
  }

  void get_sds_attr( std::string const& sds_name, std::string const& attr_name, uint16& result ) const {
    std::vector<uint16> result_vec;
    get_sds_attr( sds_name, attr_name, result_vec );
    if( result_vec.size() != 1 ) vw_throw( NotFoundErr() << "SDS attribute with requested size not found!" );
    result = result_vec[0];
  }

  void get_sds_attr( std::string const& sds_name, std::string const& attr_name, std::vector<float32>& result ) const {
    ::int32 sds_index = SDnametoindex( sd_id, sds_name.c_str() );
    if( sds_index == FAIL ) vw_throw( NotFoundErr() << "SDS not found!" );
    ::int32 sds_id = SDselect( sd_id, sds_index );
    if( sds_id == FAIL ) vw_throw( IOErr() << "Unable to select SDS in HDF file \"" << resource.filename() << "\"!" );
    ::int32 attr_index = SDfindattr( sds_id, attr_name.c_str() );
    if( attr_index == FAIL ) vw_throw( NotFoundErr() << "SDS attribute not found!" );
    char name[65];
    ::int32 data_type, count;
    ::int32 status = SDattrinfo( sds_id, attr_index, name, &data_type, &count );
    if( status == FAIL ) vw_throw( IOErr() << "Unable to get SDS attribute info in HDF file \"" << resource.filename() << "\"!" );
    if( data_type != DFNT_FLOAT32 ) vw_throw( NotFoundErr() << "SDS attribute with requested type not found! (expected " << DFNT_FLOAT32 << ", found " << data_type << ")" );
    result.resize( count );
    status = SDreadattr( sds_id, attr_index, &result[0] );
    if( status == FAIL ) vw_throw( IOErr() << "Unable to read SDS attribute in HDF file \"" << resource.filename() << "\"!" );
    vw_out(VerboseDebugMessage) << "Attribute \"" << name << "\": data type " << data_type << ", length " << count << std::endl;
    SDendaccess( sds_id );
  }

  void get_sds_attr( std::string const& sds_name, std::string const& attr_name, float32& result ) const {
    std::vector<float32> result_vec;
    get_sds_attr( sds_name, attr_name, result_vec );
    if( result_vec.size() != 1 ) vw_throw( NotFoundErr() << "SDS attribute with requested size not found!" );
    result = result_vec[0];
  }

  void get_sds_attr( std::string const& sds_name, std::string const& attr_name, std::vector<float64>& result ) const {
    ::int32 sds_index = SDnametoindex( sd_id, sds_name.c_str() );
    if( sds_index == FAIL ) vw_throw( NotFoundErr() << "SDS not found!" );
    ::int32 sds_id = SDselect( sd_id, sds_index );
    if( sds_id == FAIL ) vw_throw( IOErr() << "Unable to select SDS in HDF file \"" << resource.filename() << "\"!" );
    ::int32 attr_index = SDfindattr( sds_id, attr_name.c_str() );
    if( attr_index == FAIL ) vw_throw( NotFoundErr() << "SDS attribute not found!" );
    char name[65];
    ::int32 data_type, count;
    ::int32 status = SDattrinfo( sds_id, attr_index, name, &data_type, &count );
    if( status == FAIL ) vw_throw( IOErr() << "Unable to get SDS attribute info in HDF file \"" << resource.filename() << "\"!" );
    if( data_type != DFNT_FLOAT64 ) vw_throw( NotFoundErr() << "SDS attribute with requested type not found! (expected " << DFNT_FLOAT64 << ", found " << data_type << ")" );
    result.resize( count );
    status = SDreadattr( sds_id, attr_index, &result[0] );
    if( status == FAIL ) vw_throw( IOErr() << "Unable to read SDS attribute in HDF file \"" << resource.filename() << "\"!" );
    vw_out(VerboseDebugMessage) << "Attribute \"" << name << "\": data type " << data_type << ", length " << count << std::endl;
    SDendaccess( sds_id );
  }

  void get_sds_attr( std::string const& sds_name, std::string const& attr_name, float64& result ) const {
    std::vector<float64> result_vec;
    get_sds_attr( sds_name, attr_name, result_vec );
    if( result_vec.size() != 1 ) vw_throw( NotFoundErr() << "SDS attribute with requested size not found!" );
    result = result_vec[0];
  }

  void get_sds_attr( std::string const& sds_name, std::string const& attr_name, std::string& result ) const {
    ::int32 sds_index = SDnametoindex( sd_id, sds_name.c_str() );
    if( sds_index == FAIL ) vw_throw( NotFoundErr() << "SDS not found!" );
    ::int32 sds_id = SDselect( sd_id, sds_index );
    if( sds_id == FAIL ) vw_throw( IOErr() << "Unable to select SDS in HDF file \"" << resource.filename() << "\"!" );
    ::int32 attr_index = SDfindattr( sds_id, attr_name.c_str() );
    if( attr_index == FAIL ) vw_throw( NotFoundErr() << "SDS attribute not found!" );
    char name[65];
    ::int32 data_type, count;
    ::int32 status = SDattrinfo( sds_id, attr_index, name, &data_type, &count );
    if( status == FAIL ) vw_throw( IOErr() << "Unable to get SDS attribute info in HDF file \"" << resource.filename() << "\"!" );
    if( data_type != DFNT_CHAR8 && data_type != DFNT_CHAR8 )
      vw_throw( NotFoundErr() << "SDS attribute with requested type not found! (expected 8-bit type, found type " << data_type << ")" );
    result.resize( count );
    status = SDreadattr( sds_id, attr_index, &result[0] );
    if( count>0 && result[count-1]==0 ) result.resize(count-1); // Remove excess trailing zero
    if( status == FAIL ) vw_throw( IOErr() << "Unable to read SDS attribute in HDF file \"" << resource.filename() << "\"!" );
    vw_out(VerboseDebugMessage) << "Attribute \"" << name << "\": data type " << data_type << ", length " << count << std::endl;
    SDendaccess( sds_id );
  }
};

// *********************************************************************
// Public interface function definitions
// *********************************************************************

vw::DiskImageResourceHDF::DiskImageResourceHDF( std::string const& filename )
  : DiskImageResource( filename ),
    m_info( new DiskImageResourceInfoHDF( filename, *this ) )
{}

vw::DiskImageResourceHDF::~DiskImageResourceHDF() {}

void vw::DiskImageResourceHDF::open( std::string const& filename ) {
  boost::shared_ptr<DiskImageResourceInfoHDF> new_info( new DiskImageResourceInfoHDF( filename, *this ) );
  m_info = new_info;
}

void vw::DiskImageResourceHDF::read( ImageBuffer const& dstbuf, BBox2i const& bbox ) const {
  m_info->read( dstbuf, bbox );
}

void vw::DiskImageResourceHDF::select_sds_planes( std::vector<vw::DiskImageResourceHDF::SDSBand> const& sds_planes ) {
  m_format = m_info->select_sds_planes( sds_planes );
}

vw::DiskImageResourceHDF& vw::DiskImageResourceHDF::select_sds( std::string const& name ) {
  m_format = m_info->select_sds( name );
  return *this;
}

void vw::DiskImageResourceHDF::get_sds_fillvalue( std::string const& sds_name, float32& result ) const {
  m_info->get_sds_fillvalue( sds_name, result );
}

std::vector<vw::DiskImageResourceHDF::AttrInfo> vw::DiskImageResourceHDF::get_sds_attrs( std::string const& sds_name ) const {
  return m_info->get_sds_attrs( sds_name );
}

void vw::DiskImageResourceHDF::get_sds_attr( std::string const& sds_name, std::string const& attr_name, std::vector<uint16>& result ) const {
  m_info->get_sds_attr( sds_name, attr_name, result );
}

void vw::DiskImageResourceHDF::get_sds_attr( std::string const& sds_name, std::string const& attr_name, uint16& result ) const {
  m_info->get_sds_attr( sds_name, attr_name, result );
}

void vw::DiskImageResourceHDF::get_sds_attr( std::string const& sds_name, std::string const& attr_name, std::vector<float32>& result ) const {
  m_info->get_sds_attr( sds_name, attr_name, result );
}

void vw::DiskImageResourceHDF::get_sds_attr( std::string const& sds_name, std::string const& attr_name, float32& result ) const {
  m_info->get_sds_attr( sds_name, attr_name, result );
}

void vw::DiskImageResourceHDF::get_sds_attr( std::string const& sds_name, std::string const& attr_name, std::vector<float64>& result ) const {
  m_info->get_sds_attr( sds_name, attr_name, result );
}

void vw::DiskImageResourceHDF::get_sds_attr( std::string const& sds_name, std::string const& attr_name, float64& result ) const {
  m_info->get_sds_attr( sds_name, attr_name, result );
}

void vw::DiskImageResourceHDF::get_sds_attr( std::string const& sds_name, std::string const& attr_name, std::string& result ) const {
  m_info->get_sds_attr( sds_name, attr_name, result );
}
