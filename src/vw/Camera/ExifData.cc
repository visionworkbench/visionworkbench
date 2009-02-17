// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <cstdio>
#include <vw/Core/FundamentalTypes.h>
#include <vw/Core/Exception.h>
#include <vw/Camera/ExifData.h>
#include <boost/algorithm/string/predicate.hpp>

namespace vw {
namespace camera {

  typedef enum {
    UBYTE = 1,
    ASCII,
    USHORT,
    ULONG,
    URATIONAL,
    SBYTE,
    UNDEFINED,
    SSHORT,
    SLONG,
    SRATIONAL,
    SINGLE,
    DOUBLE
  } ExifDataFormat;
  
  const int M_SOI = 0xD8; //Start Of Image
  const int M_EOI = 0xD9; //End Of Image
  const int M_SOS = 0xDA; //Start Of Scan (begins compressed data)
  const int M_EXIF = 0xE1; //Exif marker
 
const int BytesPerFormat[] = {0,1,1,2,4,8,1,1,2,4,8,4,8};


}} // namespace vw::camera

// --------------------------------------------------------------
//                   ExifData
// --------------------------------------------------------------

vw::camera::ExifData::~ExifData() {
  typedef std::map<unsigned int, ExifTagData>::iterator iterator;
  for (iterator tag = tags.begin(); tag != tags.end(); tag++) {
    if ((*tag).second.type == StringType) {
      free((*tag).second.value.s);
    }
  }
}

// Convert a 16 bit unsigned value to file's native byte order
void vw::camera::ExifData::Put16u(void * Short, unsigned short PutValue) {
  if (MotorolaOrder){
    ((uint8 *)Short)[0] = (uint8)(PutValue>>8);
    ((uint8 *)Short)[1] = (uint8)PutValue;
  } else {
    ((uint8 *)Short)[0] = (uint8)PutValue;
    ((uint8 *)Short)[1] = (uint8)(PutValue>>8);
  }
}

// Convert a 16 bit unsigned value from file's native byte order
int vw::camera::ExifData::Get16u(void * Short) {
  if (MotorolaOrder){
    return (((uint8 *)Short)[0] << 8) | ((uint8 *)Short)[1];
  } else {
    return (((uint8 *)Short)[1] << 8) | ((uint8 *)Short)[0];
  }
}

// Convert a 32 bit signed value from file's native byte order
int vw::camera::ExifData::Get32s(void * Long) {
  if (MotorolaOrder){
    return  ((( char *)Long)[0] << 24) | (((uint8 *)Long)[1] << 16)
            | (((uint8 *)Long)[2] << 8 ) | (((uint8 *)Long)[3] << 0 );
  } else {
    return  ((( char *)Long)[3] << 24) | (((uint8 *)Long)[2] << 16)
            | (((uint8 *)Long)[1] << 8 ) | (((uint8 *)Long)[0] << 0 );
  }
}

// Convert a 32 bit unsigned value to file's native byte order
void vw::camera::ExifData::Put32u(void * Value, unsigned PutValue) {
  if (MotorolaOrder){
    ((uint8 *)Value)[0] = (uint8)(PutValue>>24);
    ((uint8 *)Value)[1] = (uint8)(PutValue>>16);
    ((uint8 *)Value)[2] = (uint8)(PutValue>>8);
    ((uint8 *)Value)[3] = (uint8)PutValue;
  } else {
    ((uint8 *)Value)[0] = (uint8)PutValue;
    ((uint8 *)Value)[1] = (uint8)(PutValue>>8);
    ((uint8 *)Value)[2] = (uint8)(PutValue>>16);
    ((uint8 *)Value)[3] = (uint8)(PutValue>>24);
  }
}

// Convert a 32 bit unsigned value from file's native byte order
unsigned vw::camera::ExifData::Get32u(void * Long) {
  return (unsigned)Get32s(Long) & 0xffffffff;
}

double vw::camera::ExifData::convert_any_format(void * ValuePtr, int Format) {
  double Value = 0;
  int Num, Den;

  switch(Format) {
    case SBYTE:
      Value = *(signed char *)ValuePtr;
      break;

    case UBYTE:
      Value = *(uint8 *)ValuePtr;
      break;

    case USHORT:
      Value = Get16u(ValuePtr);
      break;

    case ULONG:
      Value = Get32u(ValuePtr);
      break;

    case URATIONAL:
    case SRATIONAL: 
      Num = Get32s(ValuePtr);
      Den = Get32s(4+(char *)ValuePtr);
      if (Den == 0) {
	Value = 0;
      } else {
	Value = (double)Num/Den;
      }
      break;

    case SSHORT:
      Value = (signed short)Get16u(ValuePtr);
      break;

    case SLONG:
      Value = Get32s(ValuePtr);
      break;

    // Not sure if this is correct (never seen float used in Exif format)
    case SINGLE:
      Value = (double)*(float *)ValuePtr;
      break;

    case DOUBLE:
      Value = *(double *)ValuePtr;
      break;

    default:
      printf("Warning: illegal format code %d", Format);
  }
  return Value;
}

unsigned char * vw::camera::ExifData::dir_entry_addr(unsigned char * start, int entry) {
  return (start + 2 + 12 * entry);
}

void vw::camera::ExifData::process_exif_dir(unsigned char * DirStart, unsigned char * OffsetBase, 
				unsigned ExifLength, int NestingLevel) {
  VW_ASSERT( NestingLevel <= 4, IOErr() << "Maximum directory nesting exceeded (corrupt Exif header)." );

  int NumDirEntries = Get16u(DirStart);
  //printf("Number of directory entries: %i\n", NumDirEntries);

  uint8* DirEnd = dir_entry_addr(DirStart, NumDirEntries);
  if (DirEnd + 4 > (OffsetBase + ExifLength)) {
    VW_ASSERT( DirEnd+2 == OffsetBase+ExifLength || DirEnd == OffsetBase+ExifLength,
	       IOErr() << "Illegally sized directory." );
  }

  for (int de = 0; de < NumDirEntries; de++) {
    unsigned char * ValuePtr;
    uint8 * DirEntry = dir_entry_addr(DirStart, de);
    
    int Tag = Get16u(DirEntry);
    int Format = Get16u(DirEntry+2);
    int Components = Get32u(DirEntry+4);

    if ( (Format > ExifData::NUM_FORMATS) || (Format <= 0) ) {
      printf("Warning: illegal number format %d for tag %04x\n", Format, Tag);
      continue;
    }
    
    if ((unsigned)Components > 0x10000) {
      printf("Warning: illegal number of components %d for tag %04x\n", Components, Tag);
      continue;
    }

    int ByteCount = Components * BytesPerFormat[Format];

    if (ByteCount > 4) {
      unsigned OffsetVal = Get32u(DirEntry + 8);
      // If its bigger than 4 bytes, the dir entry contains an offset.
      if (OffsetVal + ByteCount > ExifLength){
	// Bogus pointer offset and / or bytecount value
	printf("Warning: illegal value pointer for tag %04x\n", Tag);
	continue;
      }
      ValuePtr = OffsetBase + OffsetVal;
    } else {
      // 4 bytes or less and value is in the dir entry itself
      ValuePtr = DirEntry + 8;
    }

    // Store tag
    switch (Format) {
      case ASCII:
      case UNDEFINED:
	// Store as string data
	tags[Tag].type = StringType;
	// Next line might not work if a tag like MakerNote includes '\0' characters
	//tags[Tag].value.s = strndup((const char *)ValuePtr, ByteCount);
	tags[Tag].value.s = (char *)malloc(ByteCount + 1);
	memcpy(tags[Tag].value.s, ValuePtr, ByteCount);
	tags[Tag].value.s[ByteCount] = '\0';
	break;

      case UBYTE:
      case USHORT:
      case ULONG:
      case SBYTE:
      case SSHORT:
      case SLONG:
	// Store as integer data
	tags[Tag].type = IntType;
	tags[Tag].value.i = (int)convert_any_format(ValuePtr, Format);
	break;

    default:
      // Store as floating point data
      tags[Tag].type = DoubleType;
      tags[Tag].value.d = convert_any_format(ValuePtr, Format);
    }
    
    // Do any special processing for specific tags
    switch (Tag) {
    case 0x8769:  // TAG_ExifOffset
    case 0xA005:  // TAG_InteroperabilityOffset:
      unsigned char * SubdirStart = OffsetBase + Get32u(ValuePtr);
      if (SubdirStart < OffsetBase || SubdirStart > OffsetBase+ExifLength){
        printf("Warning: illegal exif or interop offset directory link");
      } else {
        process_exif_dir(SubdirStart, OffsetBase, ExifLength, NestingLevel+1);
      }
      continue;
      
      // Process MakerNote
      /*
        case TAG_MakerNote:
        continue;
      */
      
      // Process GPS info
      // If you need GPS info, probably the Cartography module is a surer bet.
      /*
        case TAG_GPSInfo:
        unsigned char * SubdirStart = OffsetBase + Get32u(ValuePtr);
        if (SubdirStart < OffsetBase || SubdirStart > OffsetBase+ExifLength){
        printf("Warning: illegal GPS directory link");
        } else {
        process_gps_info(SubdirStart, ByteCount, OffsetBase, ExifLength);
        }
        continue;
      */
    }
  }

  // In addition to linking to subdirectories via exif tags, 
  // there's also a potential link to another directory at the end of each
  // directory.
  if (dir_entry_addr(DirStart, NumDirEntries) + 4 <= OffsetBase + ExifLength){
    unsigned Offset = Get32u(DirStart+2+12*NumDirEntries);
    if (Offset){
      uint8* SubdirStart = OffsetBase + Offset;
      if ((SubdirStart > OffsetBase + ExifLength) || (SubdirStart < OffsetBase)) {
	if ((SubdirStart > OffsetBase) && (SubdirStart < OffsetBase + ExifLength + 20)) {
	  // let this pass silently
	} else {
	  printf("Warning: illegal subdirectory link\n");
	}
      } else {
	if (SubdirStart <= OffsetBase + ExifLength){
	  process_exif_dir(SubdirStart, OffsetBase, ExifLength, NestingLevel+1);
	}
      }
    }
  }
}

int vw::camera::ExifData::process_tiff_header(unsigned char * buffer) {
  // Bytes 0-1 of TIFF header indicate byte order
  if (memcmp(buffer, "II", 2) == 0) {
    MotorolaOrder = 0; //Intel order
  } else {
    VW_ASSERT( memcmp(buffer, "MM", 2) == 0, IOErr() << "Invalid Exif alignment marker." );
    MotorolaOrder = 1; //Motorola order
  }

  // Sanity check of bytes 2-3 (arbitrarily always 42 according to standard)
  VW_ASSERT(Get16u(buffer+2) == 0x2a, IOErr() << "Invalid Exif start." );

  // Bytes 4-7 contain offset of first IFD
  int first_offset = Get32u(buffer+4);
  if (first_offset < 8 || first_offset > 16){
    printf("Warning: suspicious offset of first IFD value.\n");
  }
  return first_offset;
}

void vw::camera::ExifData::process_exif(unsigned char * ExifSection, unsigned int length) {

  // Check the EXIF header component
  static uint8 ExifHeader[] = "Exif\0\0";
  VW_ASSERT( memcmp(ExifSection+2, ExifHeader, 6) == 0, IOErr() << "Incorrect Exif header." );

  int first_offset = process_tiff_header(ExifSection+8);

  // First directory starts 16 bytes in.  All offset are relative to 8 bytes in.
  process_exif_dir(ExifSection + 8 + first_offset, ExifSection + 8, length - 8, 0);
}

bool vw::camera::ExifData::read_tiff_ifd(FILE* infile) {
  // Obtain file size
  fseek(infile, 0, SEEK_END);
  size_t lSize = ftell(infile);
  rewind(infile);

  // Read complete file into buffer (inefficient, but allows us to use same process_exif_dir function
  // unchanged for both jpg and tiff).
  unsigned char * buffer = (unsigned char *) malloc(lSize);
  VW_ASSERT( buffer != NULL, NullPtrErr() << "Could not allocate memory.");

  if (fread(buffer, 1, lSize, infile) != lSize)
    vw_throw(IOErr() << "File changed size!");

  int first_offset = process_tiff_header(buffer);

  process_exif_dir(buffer + first_offset, buffer, lSize, 0);

  free(buffer);

  return true;
}

bool vw::camera::ExifData::read_jpeg_sections(FILE* infile) {
  unsigned int pos = 0;
  int a = fgetc(infile);

  if (a != 0xff || fgetc(infile) != M_SOI){
    return false;
  }
  pos+=2;

  while (true) {
    int marker = 0;
    
    for (int i = 0; i < 7; i++){
      marker = fgetc(infile);
      pos++;
      if (marker != 0xff) break;
      
      VW_ASSERT( i < 6, IOErr() << "Too many padding bytes." );
    }
    
    // Read the length of the section.
    int lh = fgetc(infile);
    int ll = fgetc(infile);
    int itemlen = (lh << 8) | ll;
    
    VW_ASSERT( itemlen >= 2, IOErr() << "Invalid JPEG marker." );
    
    uint8* data = (uint8 *)malloc(itemlen);
    VW_ASSERT( data != NULL, NullPtrErr() << "Could not allocate memory." );
    
    // Store first two pre-read bytes.
    data[0] = (uint8)lh;
    data[1] = (uint8)ll;
    
    int got = fread(data+2, 1, itemlen-2, infile); // Read the whole section.
    pos += itemlen;
    VW_ASSERT( got == itemlen - 2, IOErr() << "Premature end of file." );
    
    switch(marker){
      case M_SOS:   // stop before hitting compressed data 
	free(data);
	return false;
      
      case M_EOI:   // end of input
	free(data);
	return false;
      
      case M_EXIF:
	// Make sure section is marked "Exif", as some software may use
	// marker 31 for other purposes.
	if (memcmp(data+2, "Exif", 4) == 0) {
          ExifLocation = pos - itemlen + 8;
	  process_exif(data, itemlen);
	  free(data);
	  return true;
	} else {
	  free(data);
	  return false;
	}
      
      default:
	// Skip any other sections.
	free(data);
	break;
    }
  }
  return false;
}

bool vw::camera::ExifData::import_data(std::string const &filename) {
  tags.clear();
  FILE * infile = fopen(filename.c_str(), "rb"); // Unix ignores 'b', windows needs it.
  
  VW_ASSERT( infile != NULL, IOErr() << "Cannot open file.");
  bool ret = false;
  
  // Identify file type (using suffixes)
  
  if (boost::algorithm::iends_with(filename, ".jpg") ||
      boost::algorithm::iends_with(filename, ".jpeg")) {
    // Scan the JPEG headers
    ret = read_jpeg_sections(infile);
  } else if (boost::algorithm::iends_with(filename, ".tif") ||
             boost::algorithm::iends_with(filename, ".tiff")) {
    // Process TIFF IFD structure
    ret = read_tiff_ifd(infile);
  } else {
    VW_ASSERT( 0, IOErr() << "Cannot determine file type.");
  }
  fclose(infile);
  
  return ret;
}

bool vw::camera::ExifData::get_tag_value(const uint16 tag, int &value) const {
  std::map<unsigned int, ExifTagData>::const_iterator tag_iter = tags.find(tag);
  if (tag_iter == tags.end()) return false;
  switch ((*tag_iter).second.type) {
    case IntType:
      value = (*tag_iter).second.value.i;
      break;
    case DoubleType:
      value = (int)(*tag_iter).second.value.d;
      break;
    default:
      return false;
  }
  return true;
}

bool vw::camera::ExifData::get_tag_value(const uint16 tag, double &value) const {
  std::map<unsigned int, ExifTagData>::const_iterator tag_iter = tags.find(tag);
  if (tag_iter == tags.end()) return false;
  switch ((*tag_iter).second.type) {
    case IntType:
      value = (double)(*tag_iter).second.value.i;
      break;
    case DoubleType:
      value = (*tag_iter).second.value.d;
      break;
    default:
      return false;
  }
  return true;
}

bool vw::camera::ExifData::get_tag_value(const uint16 tag, std::string &value) const {
  std::map<unsigned int, ExifTagData>::const_iterator tag_iter = tags.find(tag);
  if (tag_iter == tags.end()) return false;
  if ((*tag_iter).second.type != StringType) return false;
  value = (*tag_iter).second.value.s;
  return true;
}

unsigned int vw::camera::ExifData::get_exif_location() const {
  return ExifLocation;
}

void vw::camera::ExifData::print_debug() {
  typedef std::map<unsigned int, ExifTagData>::const_iterator iterator;
  for (iterator tag = tags.begin(); tag != tags.end(); tag++) {
    printf("Tag %04x: ", (*tag).first);
    switch ((*tag).second.type) {
    case IntType:
      printf("%i\n", (*tag).second.value.i);
      break;
    case DoubleType:
      printf("%f\n", (*tag).second.value.d);
      break;
    case StringType:
      printf("%s\n", (*tag).second.value.s);
      break;
    }
  }
}
