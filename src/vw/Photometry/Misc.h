// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// \file Misc.h

#ifndef __VW_PHOTOMETRY_MISC_H__
#define __VW_PHOTOMETRY_MISC_H__

#include <sys/stat.h>
#include <math.h>
#include <time.h>
#include <vector>

#include <vw/Image/ImageViewBase.h>
#include <vw/Cartography/GeoReference.h>
#include <vw/Photometry/Reconstruct.h>

namespace vw {
namespace photometry {
  //upsamples a DEM by an upsampleFactor
  void upsample_image(std::string output_file, std::string input_file, int upsampleFactor);

  //subsamples a geo referenced tiff image by two
  void subsample_image(std::string output_file, std::string input_file);

  // Given two images and two georeferences, this function picks a set
  // of matching pixel samples between the two images.  It rejects
  // pixels that are not valid, and it should probably also reject
  // pixels that are near saturation (though it does not yet!).
  template <class ViewT>
  std::vector<Vector4> sample_images(ImageViewBase<ViewT> const& image1,
                                     ImageViewBase<ViewT> const& image2,
                                     cartography::GeoReference const& geo1,
                                     cartography::GeoReference const& geo2,
                                     int num_samples,
                                     std::string const& DEM_file,
                                     std::vector<Vector3> *normalArray,
                                     std::vector<Vector3> *xyzArray );


  /// Erases a file suffix if one exists and returns the base string
  static std::string prefix_from_filename(std::string const& filename) {
    std::string result = filename;
    int index = result.rfind(".");
    if (index != -1)
      result.erase(index, result.size());
    return result;
  }

  /// Erases a file suffix if one exists and returns the base string less3 characters
  static std::string prefix_less3_from_filename(std::string const& filename) {
    std::string result = filename;
    int index = result.rfind(".");
    if (index != -1)
      result.erase(index-3, result.size()+3);
    return result;
  }

  /// Erases a file suffix if one exists and returns the base string less3 characters
  static std::string sufix_from_filename(std::string const& filename) {
    std::string result = filename;
    int index = result.rfind("/");
    if (index != -1)
      result.erase(0, index);
    return result;
  }

  //reads the tiff DEM into a 3D coordinate
  //pos is a Vector2 of pixel coordinates, GR is georeference
  template <class ViewT>
  Vector3 pixel_to_cart (Vector2 pos, ImageViewBase<ViewT> const& img,
                         cartography::GeoReference GR);


  std::vector<std::string> parse_command_arguments(int argc, char *argv[] );


}} // end vw::photometry

#endif//__VW_PHOTOMETRY_MISC_H__
