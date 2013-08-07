// __BEGIN_LICENSE__
//  Copyright (c) 2006-2013, United States Government as represented by the
//  Administrator of the National Aeronautics and Space Administration. All
//  rights reserved.
//
//  The NASA Vision Workbench is licensed under the Apache License,
//  Version 2.0 (the "License"); you may not use this file except in
//  compliance with the License. You may obtain a copy of the License at
//  http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
// __END_LICENSE__

/// \file ExifData.h Contains supporting code for Exif.h

// Note: The following code for parsing EXIF information in the Vision
// Workbench camera module was adapted from JHead : the EXIF Jpeg
// header manipulation tool written by Matthias Wandel
// (http://www.sentex.net/~mwandel/jhead/).  Here is the JHead
// copyright notice:
//
//    Jhead is public domain software - that is, you can do whatever
//    you want with it, and include it in software that is licensed under
//    the GNU or the BSD license, or whatever other licence you chose,
//    including proprietary closed source licenses.  Although not part
//    of the license, I do expect common courtesy, please.
//
//    -Matthias Wandel
//

#ifndef __VW_CAMERA_EXIF_DATA_H__
#define __VW_CAMERA_EXIF_DATA_H__

#include <vw/Core/FundamentalTypes.h>
#include <map>
#include <cstdio>

namespace vw {
namespace camera {


  typedef enum { IntType, DoubleType, StringType } ExifDataType;

  typedef struct {
    ExifDataType type;
    union {
      int i;
      double d;
      char* s;
    } value;
  } ExifTagData;

  class ExifData {
  private:
    static const int NUM_FORMATS = 12;

    std::map<unsigned int, ExifTagData> tags;
    bool MotorolaOrder;
    unsigned int ExifLocation;

    int process_tiff_header(const unsigned char * buffer);
    bool read_jpeg_sections(FILE* infile);
    bool read_tiff_ifd(const std::string &filename);
    void process_exif(unsigned char * ExifSection, unsigned int length);
    void process_exif_dir(const unsigned char * DirStart, const unsigned char * OffsetBase,
                          unsigned ExifLength, int NestingLevel);
    const unsigned char * dir_entry_addr(const unsigned char * start, int entry);
    void Put16u(void * Short, unsigned short PutValue);
    int Get16u(const void * Short);
    int Get32s(const void * Long);
    void Put32u(void * Value, unsigned PutValue);
    unsigned Get32u(const void * Long);
    double convert_any_format(const void * ValuePtr, int Format);

  public:
    ExifData() {}
    ~ExifData();

    bool import_data(std::string const &filename);

    bool get_tag_value(const uint16 tag, int &value) const;
    bool get_tag_value(const uint16 tag, double &value) const;
    bool get_tag_value(const uint16 tag, std::string &value) const;
    unsigned int get_exif_location() const;

    void print_debug();
  };

}} // namespace vw::camera

#endif // __VW_CAMERA_EXIF_DATA_H__
