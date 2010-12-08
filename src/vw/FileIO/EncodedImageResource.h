// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// An abstract base class referring to an encoded image
#ifndef __VW_FILEIO_ENCODEDIMAGERESOURCE_H__
#define __VW_FILEIO_ENCODEDIMAGERESOURCE_H__

#include <vw/Image/ImageResource.h>

namespace vw {

  template <class T>
  class VarArray;

  class SrcEncodedImageResource : public SrcImageResource {
    public:
      // Write encoded data to be decoded by a read() (zerocopy)
      virtual void set_encoded_data( const uint8* data, size_t len) = 0;
      // constructs the appropriate subclass for the type, then calls write_encoded
      static SrcEncodedImageResource* open( const uint8* data, size_t len, const std::string& type);
  };

  class DstEncodedImageResource : public DstImageResource {
    public:
      // Copy encoded data into data, resizing as necessary
      virtual void set_decoded_data( VarArray<uint8>& data) = 0;
      // constructs the appropriate subclass for the type, then calls write() using the data and the format
      static DstEncodedImageResource* create( const uint8* data, size_t len, const ImageFormat& format, const std::string& type);
  };


} // namespace vw

#endif
