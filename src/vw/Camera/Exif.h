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

/// \file Exif.h Contains class for reading EXIF jpeg data

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

#ifndef __VW_CAMERA_EXIF_H__
#define __VW_CAMERA_EXIF_H__

#include <vw/Core/Exception.h>
#include <vw/Core/FundamentalTypes.h>
#include <vw/Camera/ExifData.h>
#include <vw/Math/Vector.h>

#include <stddef.h>
#include <string>

namespace vw {
namespace camera {

  struct ExifDateTime {
    uint16 m_year;       // 4-digit
    uint8  m_month;      // 1-12
    uint8  m_day;        // 0-31
    uint8  m_hour;       // 0-23
    uint8  m_minute;    // 0-59
    uint8  m_second;    // 0-59
  };

  /// Convenience constants for use in querying by tag value. Names are
  /// of the form EXIF_<Field name>, with the original capitalization of
  /// the field name retained. Some obscure tags may not be listed, but
  /// can still be accessed (if stored by the camera) if their tag
  /// number in known.
  ///
  enum ExifTag {
    EXIF_InteropIndex                = 0x0001,
    EXIF_InteropVersion              = 0x0002,
    EXIF_ImageWidth                  = 0x0100,
    EXIF_ImageLength                 = 0x0101,
    EXIF_BitsPerSample               = 0x0102,
    EXIF_Compression                 = 0x0103,
    EXIF_PhotometricInterpretation   = 0x0106,
    EXIF_FillOrder                   = 0x010A,
    EXIF_DocumentName                = 0x010D,
    EXIF_ImageDescription            = 0x010E,
    EXIF_Make                        = 0x010F,
    EXIF_Model                       = 0x0110,
    EXIF_StripOffsets                = 0x0111,
    EXIF_Orientation                 = 0x0112,
    EXIF_SamplesPerPixel             = 0x0115,
    EXIF_RowsPerStrip                = 0x0116,
    EXIF_StripByteCounts             = 0x0117,
    EXIF_XResolution                 = 0x011A,
    EXIF_YResolution                 = 0x011B,
    EXIF_PlanarConfiguration         = 0x011C,
    EXIF_ResolutionUnit              = 0x0128,
    EXIF_TransferFunction            = 0x012D,
    EXIF_Software                    = 0x0131,
    EXIF_DateTime                    = 0x0132,
    EXIF_Artist                      = 0x013B,
    EXIF_WhitePoint                  = 0x013E,
    EXIF_PrimaryChromaticities       = 0x013F,
    EXIF_TransferRange               = 0x0156,
    EXIF_JPEGProc                    = 0x0200,
    EXIF_ThumbnailOffset             = 0x0201,
    EXIF_ThumbnailLength             = 0x0202,
    EXIF_YCbCrCoefficients           = 0x0211,
    EXIF_YCbCrSubSampling            = 0x0212,
    EXIF_YCbCrPositioning            = 0x0213,
    EXIF_ReferenceBlackWhite         = 0x0214,
    EXIF_RelatedImageWidth           = 0x1001,
    EXIF_RelatedImageLength          = 0x1002,
    EXIF_CFARepeatPatternDim         = 0x828D,
    EXIF_CFAPattern1                 = 0x828E,
    EXIF_BatteryLevel                = 0x828F,
    EXIF_Copyright                   = 0x8298,
    EXIF_ExposureTime                = 0x829A,
    EXIF_FNumber                     = 0x829D,
    EXIF_IPTC_NAA                    = 0x83BB,
    EXIF_ExifOffset                  = 0x8769,
    EXIF_InterColorProfile           = 0x8773,
    EXIF_ExposureProgram             = 0x8822,
    EXIF_SpectralSensitivity         = 0x8824,
    EXIF_GPSInfo                     = 0x8825,
    EXIF_ISOSpeedRatings             = 0x8827,
    EXIF_OECF                        = 0x8828,
    EXIF_ExifVersion                 = 0x9000,
    EXIF_DateTimeOriginal            = 0x9003,
    EXIF_DateTimeDigitized           = 0x9004,
    EXIF_ComponentsConfiguration     = 0x9101,
    EXIF_CompressedBitsPerPixel      = 0x9102,
    EXIF_ShutterSpeedValue           = 0x9201,
    EXIF_ApertureValue               = 0x9202,
    EXIF_BrightnessValue             = 0x9203,
    EXIF_ExposureBiasValue           = 0x9204,
    EXIF_MaxApertureValue            = 0x9205,
    EXIF_SubjectDistance             = 0x9206,
    EXIF_MeteringMode                = 0x9207,
    EXIF_LightSource                 = 0x9208,
    EXIF_Flash                       = 0x9209,
    EXIF_FocalLength                 = 0x920A,
    EXIF_SubjectArea                 = 0x9214,
    EXIF_MakerNote                   = 0x927C,
    EXIF_UserComment                 = 0x9286,
    EXIF_SubSecTime                  = 0x9290,
    EXIF_SubSecTimeOriginal          = 0x9291,
    EXIF_SubSecTimeDigitized         = 0x9292,
    EXIF_WinXPTitle                  = 0x9C9B, // windows-only - not part of exif standard.
    EXIF_WinXPComment                = 0x9C9C, // windows-only - not part of exif standard.
    EXIF_WinXPAuthor                 = 0x9C9D, // windows-only - not part of exif standard.
    EXIF_WinXPKeywords               = 0x9C9E, // windows-only - not part of exif standard.
    EXIF_WinXPSubject                = 0x9C9F, // windows-only - not part of exif standard.
    EXIF_FlashPixVersion             = 0xA000,
    EXIF_ColorSpace                  = 0xA001,
    EXIF_PixelXDimension             = 0xA002,
    EXIF_PixelYDimension             = 0xA003,
    EXIF_RelatedAudioFile            = 0xA004,
    EXIF_InteroperabilityOffset      = 0xA005,
    EXIF_FlashEnergy                 = 0xA20B,
    EXIF_SpatialFrequencyResponse    = 0xA20C,
    EXIF_FocalPlaneXResolution       = 0xA20E,
    EXIF_FocalPlaneYResolution       = 0xA20F,
    EXIF_FocalPlaneResolutionUnit    = 0xA210,
    EXIF_SubjectLocation             = 0xA214,
    EXIF_ExposureIndex               = 0xA215,
    EXIF_SensingMethod               = 0xA217,
    EXIF_FileSource                  = 0xA300,
    EXIF_SceneType                   = 0xA301,
    EXIF_CFAPattern                  = 0xA302,
    EXIF_CustomRendered              = 0xA401,
    EXIF_ExposureMode                = 0xA402,
    EXIF_WhiteBalance                = 0xA403,
    EXIF_DigitalZoomRatio            = 0xA404,
    EXIF_FocalLengthIn35mmFilm       = 0xA405,
    EXIF_SceneCaptureType            = 0xA406,
    EXIF_GainControl                 = 0xA407,
    EXIF_Contrast                    = 0xA408,
    EXIF_Saturation                  = 0xA409,
    EXIF_Sharpness                   = 0xA40A,
    EXIF_SubjectDistanceRange        = 0xA40C,
    EXIF_ImageUniqueId               = 0xA420
  };

  // Module Specific Exceptions
  VW_DEFINE_EXCEPTION(ExifErr, vw::Exception);

  class ExifView {
    ExifData m_data;

  public:
    ExifView(std::string const& filename);

    // Query the data by tag ID (common tags are enumerated at the top if
    // Exif.h)
    void query_by_tag(ExifTag tag, vw::int32& value) const;
    void query_by_tag(ExifTag tag, double& value) const;
    void query_by_tag(ExifTag tag, std::string& value) const;
    void query_by_tag(ExifTag tag, ExifDateTime& value) const;

    // Camera info
    std::string get_make() const;
    std::string get_model() const;

    // Date and time
    ExifDateTime get_modification_time() const;
    ExifDateTime get_capture_time() const;
    ExifDateTime get_digitization_time() const;

    // Camera settings
    double get_f_number() const;
    double get_exposure_time() const;
    double get_iso() const;
    double get_focal_length_35mm_equiv() const;
    Vector2 get_focal_length_pix() const;
    Vector2i get_image_size() const;

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
