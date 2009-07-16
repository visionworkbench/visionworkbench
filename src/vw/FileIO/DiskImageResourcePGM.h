// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// \file DiskImageResourcePGM.h
/// 
/// Provides support for the PGM file format.
///
#ifndef __VW_FILEIO_DISKIMAGERESOURCEPGM_H__
#define __VW_FILEIO_DISKIMAGERESOURCEPGM_H__

#include <string>
#include <fstream>
#include <boost/shared_ptr.hpp>

#include <vw/FileIO/DiskImageResource.h>

namespace vw {

  class DiskImageResourcePGM : public DiskImageResource { 
  public:
    
    // Standard DiskImageResource interface:
    
    DiskImageResourcePGM( std::string const& filename ); // Reading
    
    DiskImageResourcePGM( std::string const& filename,
			  ImageFormat const& format ); // Writing

    virtual ~DiskImageResourcePGM() {}

    // Returns the type of disk image resource.
    static std::string type_static(){ return "PGM"; }

    // Returns the type of disk image resource.
    virtual std::string type() { return type_static(); }
    
    virtual void read(ImageBuffer const& buf, BBox2i const& bbox ) const;
    virtual void write( ImageBuffer const& dest, BBox2i const& bbox );
    virtual void flush() {}

    void open( std::string const& filename );

    void create( std::string const& filename,
                 ImageFormat const& format );

    static DiskImageResource* construct_open( std::string const& filename );

    static DiskImageResource* construct_create( std::string const& filename,
                                                ImageFormat const& format );

  private:
    DiskImageResourcePGM(DiskImageResourcePGM&);
    DiskImageResourcePGM& operator=(DiskImageResourcePGM&);
    fpos_t m_image_data_position;
    
    void read_reset() const;
  };

} // namespace VW

#endif//__VW_FILEIO_DISKIMAGERESOURCEPGM_H__
