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

#include <vw/Core/Exception.h>
#include <vw/Image/PixelTypes.h>
#include <vw/FileIO/DiskImageResource.h>

namespace vw {

  class DiskImageResourcePDS : public DiskImageResource {
  public:

    DiskImageResourcePDS( std::string const& filename )
      : DiskImageResource( filename ),
        m_image_data_offset( 0 ),
        m_invalid_as_alpha( false )
    {
      open( filename );
    }
    
    DiskImageResourcePDS( std::string const& filename, 
                          ImageFormat const& /*format*/ )
      : DiskImageResource( filename )
    {
      vw_throw( NoImplErr() << "The PDS driver does not yet support creation of PDS files." );
    }
    
    virtual ~DiskImageResourcePDS() {}
    
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

  private:
    void parse_pds_header(std::vector<std::string> const& header);
    std::map<std::string, std::string> m_header_entries;
    int m_image_data_offset;
    bool m_invalid_as_alpha;
  };

} // namespace vw

#endif // __VW_FILEIO_DISKIMAGERESOUCEPDS_H__
