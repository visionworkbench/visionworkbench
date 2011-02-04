// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// An abstract base class referring to an encoded image
#ifndef __VW_FILEIO_MEMORYIMAGERESOURCE_H__
#define __VW_FILEIO_MEMORYIMAGERESOURCE_H__

#include <vw/Image/ImageResource.h>

namespace vw {

  class SrcMemoryImageResource : public SrcImageResource {
    public:
      // constructs the appropriate subclass for the type
      // Does not take ownership of data
      static SrcMemoryImageResource* open( const std::string& type, const uint8* data, size_t len);
      // Takes ownership of data
      static SrcMemoryImageResource* open( const std::string& type, boost::shared_array<const uint8> data, size_t len);
  };

  class DstMemoryImageResource : public DstImageResource {
    public:
      // constructs the appropriate subclass for the type
      static DstMemoryImageResource* create( const std::string& type, const ImageFormat& format );
      virtual const uint8* data() const = 0;
      virtual size_t size() const = 0;
  };


} // namespace vw

#endif
