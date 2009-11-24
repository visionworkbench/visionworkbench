// __BEGIN_LICENSE__
// Copyright (C) 2006-2009 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#include <vw/Plate/PlateFile.h>
using namespace vw::platefile;
using namespace vw;

PlateFile::PlateFile(std::string url) {
    m_index = boost::shared_ptr<Index>( new LocalIndex(url) );
    vw_out(DebugMessage, "platefile") << "Re-opened plate file: \"" << url << "\"\n";
}

PlateFile::PlateFile(std::string url, std::string type, std::string description,
                     int tile_size, std::string tile_filetype, 
                     PixelFormatEnum pixel_format, ChannelTypeEnum channel_type) {
  
  // Plate files are stored as an index file and one or more data
  // blob files in a directory.  We create that directory here if
  // it doesn't already exist.
  if( !exists( fs::path( url, fs::native ) ) ) {
    fs::create_directory(url);

    IndexHeader hdr;
    hdr.set_type(type);
    hdr.set_description(description);
    hdr.set_tile_size(tile_size);
    hdr.set_tile_filetype(tile_filetype);
    hdr.set_pixel_format(pixel_format);
    hdr.set_channel_type(channel_type);

    m_index = boost::shared_ptr<Index>( new LocalIndex(url, hdr) );
    vw_out(DebugMessage, "platefile") << "Creating new plate file: \"" << url << "\"\n";

    // However, if it does exist, then we attempt to open the
    // platefile that is stored there.
  } else {
    m_index = boost::shared_ptr<Index>( new LocalIndex(url) );
    vw_out(DebugMessage, "platefile") << "Re-opened plate file: \"" << url << "\"\n";
  }
      
}
