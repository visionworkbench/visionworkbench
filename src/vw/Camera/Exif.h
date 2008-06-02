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
#ifndef __VW_CAMERA_EXIF_H__
#define __VW_CAMERA_EXIF_H__

#include <vw/Core/Exception.h>
#include <vw/Core/FundamentalTypes.h>
#include <vw/Camera/ExifData.h>

namespace vw { 
namespace camera {

  /// Convenience constants for use in querying by tag value. Names are
  /// of the form EXIF_<Field name>, with the original capitalization of
  /// the field name retained. Some obscure tags may not be listed, but
  /// can still be accessed (if stored by the camera) if their tag
  /// number in known.
  ///
  static const uint16 EXIF_InteropIndex                = 0x0001;
  static const uint16 EXIF_InteropVersion              = 0x0002;
  static const uint16 EXIF_ImageWidth                  = 0x0100;
  static const uint16 EXIF_ImageLength                 = 0x0101;
  static const uint16 EXIF_BitsPerSample               = 0x0102;
  static const uint16 EXIF_Compression                 = 0x0103;
  static const uint16 EXIF_PhotometricInterpretation   = 0x0106;
  static const uint16 EXIF_FillOrder                   = 0x010A;
  static const uint16 EXIF_DocumentName                = 0x010D;
  static const uint16 EXIF_ImageDescription            = 0x010E;
  static const uint16 EXIF_Make                        = 0x010F;
  static const uint16 EXIF_Model                       = 0x0110;
  static const uint16 EXIF_StripOffsets                = 0x0111;
  static const uint16 EXIF_Orientation                 = 0x0112;
  static const uint16 EXIF_SamplesPerPixel             = 0x0115;
  static const uint16 EXIF_RowsPerStrip                = 0x0116;
  static const uint16 EXIF_StripByteCounts             = 0x0117;
  static const uint16 EXIF_XResolution                 = 0x011A;
  static const uint16 EXIF_YResolution                 = 0x011B;
  static const uint16 EXIF_PlanarConfiguration         = 0x011C;
  static const uint16 EXIF_ResolutionUnit              = 0x0128;
  static const uint16 EXIF_TransferFunction            = 0x012D;
  static const uint16 EXIF_Software                    = 0x0131;
  static const uint16 EXIF_DateTime                    = 0x0132;
  static const uint16 EXIF_Artist                      = 0x013B;
  static const uint16 EXIF_WhitePoint                  = 0x013E;
  static const uint16 EXIF_PrimaryChromaticities       = 0x013F;
  static const uint16 EXIF_TransferRange               = 0x0156;
  static const uint16 EXIF_JPEGProc                    = 0x0200;
  static const uint16 EXIF_ThumbnailOffset             = 0x0201;
  static const uint16 EXIF_ThumbnailLength             = 0x0202;
  static const uint16 EXIF_YCbCrCoefficients           = 0x0211;
  static const uint16 EXIF_YCbCrSubSampling            = 0x0212;
  static const uint16 EXIF_YCbCrPositioning            = 0x0213;
  static const uint16 EXIF_ReferenceBlackWhite         = 0x0214;
  static const uint16 EXIF_RelatedImageWidth           = 0x1001;
  static const uint16 EXIF_RelatedImageLength          = 0x1002;
  static const uint16 EXIF_CFARepeatPatternDim         = 0x828D;
  //  static const uint16 EXIF_CFAPattern                 =  0x828E;
  static const uint16 EXIF_BatteryLevel                = 0x828F;
  static const uint16 EXIF_Copyright                   = 0x8298;
  static const uint16 EXIF_ExposureTime                = 0x829A;
  static const uint16 EXIF_FNumber                     = 0x829D;
  static const uint16 EXIF_IPTC_NAA                    = 0x83BB;
  static const uint16 EXIF_ExifOffset                  = 0x8769;
  static const uint16 EXIF_InterColorProfile           = 0x8773;
  static const uint16 EXIF_ExposureProgram             = 0x8822;
  static const uint16 EXIF_SpectralSensitivity         = 0x8824;
  static const uint16 EXIF_GPSInfo                     = 0x8825;
  static const uint16 EXIF_ISOSpeedRatings             = 0x8827;
  static const uint16 EXIF_OECF                        = 0x8828;
  static const uint16 EXIF_ExifVersion                 = 0x9000;
  static const uint16 EXIF_DateTimeOriginal            = 0x9003;
  static const uint16 EXIF_DateTimeDigitized           = 0x9004;
  static const uint16 EXIF_ComponentsConfiguration     = 0x9101;
  static const uint16 EXIF_CompressedBitsPerPixel      = 0x9102;
  static const uint16 EXIF_ShutterSpeedValue           = 0x9201;
  static const uint16 EXIF_ApertureValue               = 0x9202;
  static const uint16 EXIF_BrightnessValue             = 0x9203;
  static const uint16 EXIF_ExposureBiasValue           = 0x9204;
  static const uint16 EXIF_MaxApertureValue            = 0x9205;
  static const uint16 EXIF_SubjectDistance             = 0x9206;
  static const uint16 EXIF_MeteringMode                = 0x9207;
  static const uint16 EXIF_LightSource                 = 0x9208;
  static const uint16 EXIF_Flash                       = 0x9209;
  static const uint16 EXIF_FocalLength                 = 0x920A;
  static const uint16 EXIF_MakerNote                   = 0x927C;
  static const uint16 EXIF_UserComment                 = 0x9286;
  static const uint16 EXIF_SubSecTime                  = 0x9290;
  static const uint16 EXIF_SubSecTimeOriginal          = 0x9291;
  static const uint16 EXIF_SubSecTimeDigitized         = 0x9292;
  static const uint16 EXIF_FlashPixVersion             = 0xA000;
  static const uint16 EXIF_ColorSpace                  = 0xA001;
  static const uint16 EXIF_PixelXDimension             = 0xA002;
  static const uint16 EXIF_PixelYDimension             = 0xA003;
  static const uint16 EXIF_RelatedAudioFile            = 0xA004;
  static const uint16 EXIF_InteroperabilityOffset      = 0xA005;
  static const uint16 EXIF_FlashEnergy                 = 0xA20B;
  static const uint16 EXIF_SpatialFrequencyResponse    = 0xA20C;
  static const uint16 EXIF_FocalPlaneXResolution       = 0xA20E;
  static const uint16 EXIF_FocalPlaneYResolution       = 0xA20F;
  static const uint16 EXIF_FocalPlaneResolutionUnit    = 0xA210;
  static const uint16 EXIF_SubjectLocation             = 0xA214;
  static const uint16 EXIF_ExposureIndex               = 0xA215;
  static const uint16 EXIF_SensingMethod               = 0xA217;
  static const uint16 EXIF_FileSource                  = 0xA300;
  static const uint16 EXIF_SceneType                   = 0xA301;
  static const uint16 EXIF_CFAPattern                  = 0xA302;
  static const uint16 EXIF_CustomRendered              = 0xA401;
  static const uint16 EXIF_ExposureMode                = 0xA402;
  static const uint16 EXIF_WhiteBalance                = 0xA403;
  static const uint16 EXIF_DigitalZoomRatio            = 0xA404;
  static const uint16 EXIF_FocalLengthIn35mmFilm       = 0xA405;
  static const uint16 EXIF_SceneCaptureType            = 0xA406;
  static const uint16 EXIF_GainControl                 = 0xA407;
  static const uint16 EXIF_Contrast                    = 0xA408;
  static const uint16 EXIF_Saturation                  = 0xA409;
  static const uint16 EXIF_Sharpness                   = 0xA40a;
  static const uint16 EXIF_SubjectDistanceRange        = 0xA40c;

  // Module Specific Exceptions
  VW_DEFINE_EXCEPTION(ExifErr, vw::Exception);

  class ExifView {
    ExifData m_data;
    
  public:
    ExifView(std::string const& filename);

    // Query the data by tag ID (common tags are enumerated at the top if
    // Exif.h)
    void query_by_tag(const uint16 tag, int& value) const;
    void query_by_tag(const uint16 tag, double& value) const;
    void query_by_tag(const uint16 tag, std::string& value) const;

    
    // Camera info
    std::string get_make() const;
    std::string get_model() const;
    
    // Camera settings
    double get_f_number() const;
    double get_exposure_time() const;
    double get_iso() const;
    double get_focal_length_35mm_equiv() const;
    
    // APEX equivalents
    //
    // These functions report values in the Additive System for
    // Photographic Exposure.  For more information, see:
    // 
    //    http://en.wikipedia.org/wiki/APEX_system

    // This function returns the value for brightness using the linear
    // system
    double get_average_luminance() const;

    // These functions report values in the additive logarithmic
    // system.
    double get_aperture_value() const;
    double get_time_value() const;
    double get_exposure_value() const;
    double get_film_speed_value() const;
    double get_luminance_value() const;
    bool   get_auto_exposure_enabled() const;

    // Location of thumbnail image in file
    size_t get_thumbnail_location() const;
  };

}} // namespace vw::camera
  
#endif  // __VW_CAMERA_EXIF_H__
