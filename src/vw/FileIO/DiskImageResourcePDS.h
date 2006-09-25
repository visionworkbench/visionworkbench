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

/// \file DiskImageResourcePDS.h
/// 
/// Provides support for some NASA mission data from the Planetary
/// Data System (PDS).
///
#ifndef __VW_FILEIO_DISK_IMAGE_RESOUCE_PDS_H__
#define __VW_FILEIO_DISK_IMAGE_RESOUCE_PDS_H__

#include <map>
#include <string>
#include <fstream>

#include <vw/Core/Exception.h>
#include <vw/Image/PixelTypes.h>
#include <vw/Image/GenericImageBuffer.h>
#include <vw/FileIO/DiskImageResource.h>

namespace vw {

  class DiskImageResourcePDS : public DiskImageResource {
  public:

    DiskImageResourcePDS( std::string const& filename )
    {
      open( filename );
    }
    
    DiskImageResourcePDS( std::string const& filename, 
                          GenericImageFormat const& format )
    {
      create( filename, format );
    }
    
    virtual ~DiskImageResourcePDS();
    
    virtual void read( GenericImageBuffer const& dest ) const;
    virtual void write( GenericImageBuffer const& dest );
    virtual void flush();

    /// Query for a value in the PDS header.  The returned value will
    /// be a string regardless of whether or not the value is a
    /// numerical type (i.e. it's up to you to convert to numerical
    /// types from a string where appropriate).  If no value is found
    /// for a matching key, a vw::NotFoundErr will be thrown.
    std::string query(std::string const& key) {
      std::map<std::string, std::string>::iterator query_object = m_header_entries.find(key);
      if (query_object != m_header_entries.end()) 
        return (*query_object).second;
      else 
        throw NotFoundErr() << "DiskImageResourcePDS: no matching value found for \"" << key << "\""; 
    }

    /// Search for a value that might match several keys.  The first
    /// key with a value in the table will be used to extract that
    /// value, so order does matter.  If no value is found for a
    /// matching key, a vw::NotFoundErr will be thrown.
    std::string query(std::vector<std::string> const& keys) {
      for (int i = 0; i < keys.size(); i++) {
        std::map<std::string, std::string>::iterator query_object = m_header_entries.find(keys[i]);
        if (query_object != m_header_entries.end()) 
          return (*query_object).second;
      }
      throw NotFoundErr() << "DiskImageResourcePDS: no matching value found for the keys provided."; 
    }


    void open( std::string const& filename );

    void create( std::string const& filename,
                 GenericImageFormat const& format );

    static DiskImageResource* construct_open( std::string const& filename );

    static DiskImageResource* construct_create( std::string const& filename,
                                                GenericImageFormat const& format );

  private:
    void parse_pds_header(std::vector<std::string> const& header);
    std::string m_filename;
    std::map<std::string, std::string> m_header_entries;
    int m_image_data_offset;
  };

} // namespace vw

#endif // __VW_FILEIO_DISK_IMAGE_RESOUCE_PDS_H__
