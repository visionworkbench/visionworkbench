// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// \file DiskImageResourceHDF.h
///
/// Provides support for the HDF file format.
///
#ifndef __VW_FILEIO_DISKIMAGERESOUCEHDF_H__
#define __VW_FILEIO_DISKIMAGERESOUCEHDF_H__

#include <string>
#include <boost/shared_ptr.hpp>

#include <vw/FileIO/DiskImageResource.h>

namespace vw {

  class DiskImageResourceInfoHDF;

  class DiskImageResourceHDF : public DiskImageResource {
  public:

    // The standard DiskImageResource interface:

    DiskImageResourceHDF( std::string const& filename );

    DiskImageResourceHDF( std::string const& filename, ImageFormat const& /*format*/ ) : DiskImageResource( filename ) {
      vw_throw( NoImplErr() << "Creating HDF files is not yet supported!" );
    }

    virtual ~DiskImageResourceHDF();

    /// Returns the type of disk image resource.
    static std::string type_static() { return "HDF"; }

    /// Returns the type of disk image resource.
    virtual std::string type() { return type_static(); }

    virtual void read( ImageBuffer const& buf, BBox2i const& bbox ) const;

    virtual void write( ImageBuffer const& /*dest*/, BBox2i const& /*bbox*/ ) {
      vw_throw( NoImplErr() << "Writing to HDF files is not yet supported!" );
    }

    void open( std::string const& filename );

    void create( std::string const& /*filename*/, ImageFormat const& /*format*/ ) {
      vw_throw( NoImplErr() << "Creating HDF files is not yet supported!" );
    }

    static DiskImageResource* construct_open( std::string const& filename );

    virtual bool has_block_write()  const {return false;}
    virtual bool has_nodata_write() const {return false;}
    virtual bool has_block_read()   const {return false;}
    virtual bool has_nodata_read()  const {return false;}

    // The HDF-specific interface:

    struct SDSInfo {
      std::string name;
      ChannelTypeEnum type;
      int32 rank;
      std::vector<int32> dim_sizes;
      bool coord;
      int32 n_attrs;
    };

    typedef std::vector<SDSInfo>::const_iterator sds_iterator;
    sds_iterator sds_begin() const;
    sds_iterator sds_end() const;

    struct AttrInfo {
      std::string name;
      ChannelTypeEnum type;
      int size;
    };

    struct SDSBand {
      std::string name;
      int band;
      SDSBand( std::string const& name, int band=0 ) : name(name), band(band) {}
    };

    void select_sds_planes( std::vector<SDSBand> const& sds_planes );

    DiskImageResourceHDF& select_sds( std::string const& name );

    void get_sds_fillvalue( std::string const& sds_name, float32& result ) const;

    std::vector<AttrInfo> get_sds_attrs( std::string const& sds_name ) const;

    void get_sds_attr( std::string const& sds_name, std::string const& attr_name, std::vector<int8>& result ) const;
    void get_sds_attr( std::string const& sds_name, std::string const& attr_name, int8& result ) const;
    void get_sds_attr( std::string const& sds_name, std::string const& attr_name, std::vector<uint8>& result ) const;
    void get_sds_attr( std::string const& sds_name, std::string const& attr_name, uint8& result ) const;
    void get_sds_attr( std::string const& sds_name, std::string const& attr_name, std::vector<int16>& result ) const;
    void get_sds_attr( std::string const& sds_name, std::string const& attr_name, int16& result ) const;
    void get_sds_attr( std::string const& sds_name, std::string const& attr_name, std::vector<uint16>& result ) const;
    void get_sds_attr( std::string const& sds_name, std::string const& attr_name, uint16& result ) const;
    void get_sds_attr( std::string const& sds_name, std::string const& attr_name, std::vector<int32>& result ) const;
    void get_sds_attr( std::string const& sds_name, std::string const& attr_name, int32& result ) const;
    void get_sds_attr( std::string const& sds_name, std::string const& attr_name, std::vector<uint32>& result ) const;
    void get_sds_attr( std::string const& sds_name, std::string const& attr_name, uint32& result ) const;
    void get_sds_attr( std::string const& sds_name, std::string const& attr_name, std::vector<float32>& result ) const;
    void get_sds_attr( std::string const& sds_name, std::string const& attr_name, float32& result ) const;
    void get_sds_attr( std::string const& sds_name, std::string const& attr_name, std::vector<float64>& result ) const;
    void get_sds_attr( std::string const& sds_name, std::string const& attr_name, float64& result ) const;
    void get_sds_attr( std::string const& sds_name, std::string const& attr_name, std::string& result ) const;

  private:
    boost::shared_ptr<DiskImageResourceInfoHDF> m_info;
  };

} // namespace vw

#endif // __VW_FILEIO_DISKIMAGERESOUCEHDF_H__
