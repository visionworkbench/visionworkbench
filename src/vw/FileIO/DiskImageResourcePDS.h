// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// \file DiskImageResourcePDS.h
///
/// Provides support for some NASA mission data from the Planetary
/// Data System (PDS).
///
#ifndef __VW_FILEIO_DISKIMAGERESOUCEPDS_H__
#define __VW_FILEIO_DISKIMAGERESOUCEPDS_H__

#include <map>
#include <string>
#include <fstream>

#include <vw/FileIO/DiskImageResource.h>

namespace vw {

  class DiskImageResourcePDS : public DiskImageResource {
  public:

    DiskImageResourcePDS( std::string const& filename )
      : DiskImageResource( filename ),
        m_image_data_offset( 0 ),
        m_invalid_as_alpha( false )
    {
      m_pds_data_filename = ""; // For PDS images that have a seperate file with image data.
      open( filename );
    }

    DiskImageResourcePDS( std::string const& filename,
                          ImageFormat const& /*format*/ )
      : DiskImageResource( filename )
    {
      vw_throw( NoImplErr() << "The PDS driver does not yet support creation of PDS files." );
    }

    virtual ~DiskImageResourcePDS() {}

    /// Returns the type of disk image resource.
    static std::string type_static() { return "PDS"; }

    /// Returns the type of disk image resource.
    virtual std::string type() { return type_static(); }

    virtual void read( ImageBuffer const& dest, BBox2i const& bbox ) const;
    virtual void write( ImageBuffer const& dest, BBox2i const& bbox );
    virtual void flush() {}

    /// Query for a string value in the PDS header.  Places the value
    /// in the result field and returns true if the value is found,
    /// otherwise returns false.
    bool query( std::string const& key, std::string& result ) const {
      std::map<std::string, std::string>::const_iterator query_object = m_header_entries.find(key);
      if (query_object != m_header_entries.end()) {
        result = (*query_object).second;
        return true;
      }
      return false;
    }

    /// Query for a floating-point value in the PDS header.  Places
    /// the value in the result field and returns true if the value is
    /// found, otherwise returns false.
    bool query( std::string const& key, float& result ) const {
      std::string result_str;
      if( query( key, result_str ) ) {
        result = (float)atof( result_str.c_str() );
        return true;
      }
      return false;
    }

    /// Search for a string value that might match several keys.  The
    /// first key with a value in the table will be used to extract
    /// that value, so order does matter.  Places the value in the
    /// result field and returns true if a matching field is found,
    /// otherwise returns false.
    bool query( std::vector<std::string> const& keys, std::string& result ) const {
      for( unsigned i = 0; i < keys.size(); ++i ) {
        std::map<std::string, std::string>::const_iterator query_object = m_header_entries.find(keys[i]);
        if( query_object != m_header_entries.end() ) {
          result = (*query_object).second;
          return true;
        }
      }
      return false;
    }

    /// Search for a floating-point value that might match several
    /// keys.  The first key with a value in the table will be used to
    /// extract that value, so order does matter.  Places the value in
    /// the result field and returns true if a matching field is
    /// found, otherwise returns false.
    bool query( std::vector<std::string> const& keys, float& result ) const {
      std::string result_str;
      if( query( keys, result_str ) ) {
        result = (float)atof( result_str.c_str() );
        return true;
      }
      return false;
    }

    /// Configure the resource to convert data validity information
    /// to an alpha channel (i.e. make missing data transparent).
    void treat_invalid_data_as_alpha();


    void open( std::string const& filename );

    void create( std::string const& filename,
                 ImageFormat const& format );

    static DiskImageResource* construct_open( std::string const& filename );

    static DiskImageResource* construct_create( std::string const& filename,
                                                ImageFormat const& format );

    virtual bool has_block_write()  const {return false;}
    virtual bool has_nodata_write() const {return false;}
    virtual bool has_block_read()   const {return false;}
    virtual bool has_nodata_read()  const {return false;}

  private:
    void parse_pds_header(std::vector<std::string> const& header);
    PixelFormatEnum planes_to_pixel_format(int32 planes) const;
    std::map<std::string, std::string> m_header_entries;
    int m_image_data_offset;
    int m_native_num_planes;
    bool m_invalid_as_alpha;
    bool m_file_is_msb_first;
    std::string m_pds_data_filename;
    enum { BAND_SEQUENTIAL, SAMPLE_INTERLEAVED, LINE_INTERLEAVED } m_band_storage;
  };

} // namespace vw

#endif // __VW_FILEIO_DISKIMAGERESOUCEPDS_H__
