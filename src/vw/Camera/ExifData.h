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

  /*
   * Convenience constants for use with get_tag_value. Names are of the form
   * TAG_<Field name>, with the original capitalization of the field name
   * retained. Some obscure tags may not be listed, but can still be
   * accessed (if stored by the camera) by their tag number with get_tag_value.
   */
  static const uint16 TAG_InteropIndex                = 0x0001;
  static const uint16 TAG_InteropVersion              = 0x0002;
  static const uint16 TAG_ImageWidth                  = 0x0100;
  static const uint16 TAG_ImageLength                 = 0x0101;
  static const uint16 TAG_BitsPerSample               = 0x0102;
  static const uint16 TAG_Compression                 = 0x0103;
  static const uint16 TAG_PhotometricInterpretation   = 0x0106;
  static const uint16 TAG_FillOrder                   = 0x010A;
  static const uint16 TAG_DocumentName                = 0x010D;
  static const uint16 TAG_ImageDescription            = 0x010E;
  static const uint16 TAG_Make                        = 0x010F;
  static const uint16 TAG_Model                       = 0x0110;
  static const uint16 TAG_StripOffsets                = 0x0111;
  static const uint16 TAG_Orientation                 = 0x0112;
  static const uint16 TAG_SamplesPerPixel             = 0x0115;
  static const uint16 TAG_RowsPerStrip                = 0x0116;
  static const uint16 TAG_StripByteCounts             = 0x0117;
  static const uint16 TAG_XResolution                 = 0x011A;
  static const uint16 TAG_YResolution                 = 0x011B;
  static const uint16 TAG_PlanarConfiguration         = 0x011C;
  static const uint16 TAG_ResolutionUnit              = 0x0128;
  static const uint16 TAG_TransferFunction            = 0x012D;
  static const uint16 TAG_Software                    = 0x0131;
  static const uint16 TAG_DateTime                    = 0x0132;
  static const uint16 TAG_Artist                      = 0x013B;
  static const uint16 TAG_WhitePoint                  = 0x013E;
  static const uint16 TAG_PrimaryChromaticities       = 0x013F;
  static const uint16 TAG_TransferRange               = 0x0156;
  static const uint16 TAG_JPEGProc                    = 0x0200;
  static const uint16 TAG_ThumbnailOffset             = 0x0201;
  static const uint16 TAG_ThumbnailLength             = 0x0202;
  static const uint16 TAG_YCbCrCoefficients           = 0x0211;
  static const uint16 TAG_YCbCrSubSampling            = 0x0212;
  static const uint16 TAG_YCbCrPositioning            = 0x0213;
  static const uint16 TAG_ReferenceBlackWhite         = 0x0214;
  static const uint16 TAG_RelatedImageWidth           = 0x1001;
  static const uint16 TAG_RelatedImageLength          = 0x1002;
  static const uint16 TAG_CFARepeatPatternDim         = 0x828D;
  //static const uint16 TAG_CFAPattern                 =  0x828E;
  static const uint16 TAG_BatteryLevel                = 0x828F;
  static const uint16 TAG_Copyright                   = 0x8298;
  static const uint16 TAG_ExposureTime                = 0x829A;
  static const uint16 TAG_FNumber                     = 0x829D;
  static const uint16 TAG_IPTC_NAA                    = 0x83BB;
  static const uint16 TAG_ExifOffset                  = 0x8769;
  static const uint16 TAG_InterColorProfile           = 0x8773;
  static const uint16 TAG_ExposureProgram             = 0x8822;
  static const uint16 TAG_SpectralSensitivity         = 0x8824;
  static const uint16 TAG_GPSInfo                     = 0x8825;
  static const uint16 TAG_ISOSpeedRatings             = 0x8827;
  static const uint16 TAG_OECF                        = 0x8828;
  static const uint16 TAG_ExifVersion                 = 0x9000;
  static const uint16 TAG_DateTimeOriginal            = 0x9003;
  static const uint16 TAG_DateTimeDigitized           = 0x9004;
  static const uint16 TAG_ComponentsConfiguration     = 0x9101;
  static const uint16 TAG_CompressedBitsPerPixel      = 0x9102;
  static const uint16 TAG_ShutterSpeedValue           = 0x9201;
  static const uint16 TAG_ApertureValue               = 0x9202;
  static const uint16 TAG_BrightnessValue             = 0x9203;
  static const uint16 TAG_ExposureBiasValue           = 0x9204;
  static const uint16 TAG_MaxApertureValue            = 0x9205;
  static const uint16 TAG_SubjectDistance             = 0x9206;
  static const uint16 TAG_MeteringMode                = 0x9207;
  static const uint16 TAG_LightSource                 = 0x9208;
  static const uint16 TAG_Flash                       = 0x9209;
  static const uint16 TAG_FocalLength                 = 0x920A;
  static const uint16 TAG_MakerNote                   = 0x927C;
  static const uint16 TAG_UserComment                 = 0x9286;
  static const uint16 TAG_SubSecTime                  = 0x9290;
  static const uint16 TAG_SubSecTimeOriginal          = 0x9291;
  static const uint16 TAG_SubSecTimeDigitized         = 0x9292;
  static const uint16 TAG_FlashPixVersion             = 0xA000;
  static const uint16 TAG_ColorSpace                  = 0xA001;
  static const uint16 TAG_PixelXDimension             = 0xA002;
  static const uint16 TAG_PixelYDimension             = 0xA003;
  static const uint16 TAG_RelatedAudioFile            = 0xA004;
  static const uint16 TAG_InteroperabilityOffset      = 0xA005;
  static const uint16 TAG_FlashEnergy                 = 0xA20B;
  static const uint16 TAG_SpatialFrequencyResponse    = 0xA20C;
  static const uint16 TAG_FocalPlaneXResolution       = 0xA20E;
  static const uint16 TAG_FocalPlaneYResolution       = 0xA20F;
  static const uint16 TAG_FocalPlaneResolutionUnit    = 0xA210;
  static const uint16 TAG_SubjectLocation             = 0xA214;
  static const uint16 TAG_ExposureIndex               = 0xA215;
  static const uint16 TAG_SensingMethod               = 0xA217;
  static const uint16 TAG_FileSource                  = 0xA300;
  static const uint16 TAG_SceneType                   = 0xA301;
  static const uint16 TAG_CFAPattern                  = 0xA302;
  static const uint16 TAG_CustomRendered              = 0xA401;
  static const uint16 TAG_ExposureMode                = 0xA402;
  static const uint16 TAG_WhiteBalance                = 0xA403;
  static const uint16 TAG_DigitalZoomRatio            = 0xA404;
  static const uint16 TAG_FocalLengthIn35mmFilm       = 0xA405;
  static const uint16 TAG_SceneCaptureType            = 0xA406;
  static const uint16 TAG_GainControl                 = 0xA407;
  static const uint16 TAG_Contrast                    = 0xA408;
  static const uint16 TAG_Saturation                  = 0xA409;
  static const uint16 TAG_Sharpness                   = 0xA40a;
  static const uint16 TAG_SubjectDistanceRange        = 0xA40c;

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
    
    bool import_data(const char* filename);
    
    bool get_tag_value(unsigned int tag, int &value);
    bool get_tag_value(unsigned int tag, double &value);
    bool get_tag_value(unsigned int tag, char* &value);
    
    void print_debug();
  };

}} // namespace vw::camera

#endif // __VW_CAMERA_EXIF_DATA_H__
