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
#ifndef __VW_CAMERA_EXIF_DATA_H__
#define __VW_CAMERA_EXIF_DATA_H__

#include <vw/Core/FundamentalTypes.h>
#include <map> 

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
    
    int process_tiff_header(unsigned char * buffer);
    bool read_jpeg_sections(FILE* infile);
    bool read_tiff_ifd(FILE* infile);
    void process_exif(unsigned char * ExifSection, unsigned int length);
    void process_exif_dir(unsigned char * DirStart, unsigned char * OffsetBase, 
                          unsigned ExifLength, int NestingLevel);
    unsigned char * dir_entry_addr(unsigned char * start, int entry);
    void Put16u(void * Short, unsigned short PutValue);
    int Get16u(void * Short);
    int Get32s(void * Long);
    void Put32u(void * Value, unsigned PutValue);
    unsigned Get32u(void * Long);
    double convert_any_format(void * ValuePtr, int Format);
    
  public:
    ExifData() {}
    ~ExifData();
    
    bool import_data(std::string const &filename);
    
    bool get_tag_value(const uint16 tag, int &value);
    bool get_tag_value(const uint16 tag, double &value);
    bool get_tag_value(const uint16 tag, std::string &value);
    unsigned int get_exif_location();
    
    void print_debug();
  };

}} // namespace vw::camera

#endif // __VW_CAMERA_EXIF_DATA_H__
