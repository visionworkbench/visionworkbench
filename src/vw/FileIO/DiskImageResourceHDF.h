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

/// \file DiskImageResourceHDF.h
/// 
/// Provides support for the HDF file format.
///
#ifndef __VW_FILEIO_DISKIMAGERESOUCEHDF_H__
#define __VW_FILEIO_DISKIMAGERESOUCEHDF_H__

#include <vw/config.h>

#include <string>
#include <boost/shared_ptr.hpp>

#include <vw/Image/PixelTypes.h>
#include <vw/FileIO/DiskImageResource.h>

namespace vw {

  class DiskImageResourceInfoHDF;

  class DiskImageResourceHDF : public DiskImageResource {
  public:

    // The standard DiskImageResource interface:

    DiskImageResourceHDF( std::string const& filename );

    DiskImageResourceHDF( std::string const& filename, ImageFormat const& format ) : DiskImageResource( filename ) {
      vw_throw( NoImplErr() << "Creating HDF files is not yet supported!" );
    }

    virtual ~DiskImageResourceHDF();
    
    /// Returns the type of disk image resource.
    static std::string type_static() { return "HDF"; }

    /// Returns the type of disk image resource.
    virtual std::string type() { return type_static(); }
    
    virtual void read( ImageBuffer const& buf, BBox2i const& bbox ) const;

    virtual void write( ImageBuffer const& dest, BBox2i const& bbox ) {
      vw_throw( NoImplErr() << "Writing to HDF files is not yet supported!" );
    }

    void open( std::string const& filename );

    void create( std::string const& filename, ImageFormat const& format ) {
      vw_throw( NoImplErr() << "Creating HDF files is not yet supported!" );
    }

    // The HDF-specific interface:

    struct SDSBand {
      std::string name;
      int band;
      SDSBand( std::string const& name, int band=0 ) : name(name), band(band) {}
    };

    struct AttrInfo {
      std::string name;
      ChannelTypeEnum type;
      int size;
    };

    void select_sds_planes( std::vector<SDSBand> const& sds_planes );

    DiskImageResourceHDF& select_sds( std::string const& name );

    void get_sds_fillvalue( std::string const& sds_name, float32& result ) const;

    std::vector<AttrInfo> get_sds_attrs( std::string const& sds_name ) const;

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
