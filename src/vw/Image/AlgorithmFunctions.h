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


/// \file AlgorithmFunctions.h
///
/// Basic algorithms operating on images. This is for functions that are non lazy.
///
#ifndef __VW_IMAGE_ALGORITHM_FUNCTIONS_H__
#define __VW_IMAGE_ALGORITHM_FUNCTIONS_H__

#include <vw/Image/ImageView.h>
#include <vw/Image/MaskViews.h>

namespace vw {

  // *******************************************************************
  // fill()
  // *******************************************************************

  /// Fill an image with a constant pixel value
  template <class ImageT, class ValueT>
  void fill( ImageViewBase<ImageT> const& image, ValueT value_ );

  // *******************************************************************
  // grassfire()
  // *******************************************************************

  /// Computes the 4-connected grassfire image of an image.
  /// (Specifically, computes the Manhattan distance from each pixel to
  /// the nearest pixel with zero value, assuming the borders of the
  /// image are zero.)
  /// - If ignore_borders is set, borders are not treated as zero value.
  template <class SourceT, class OutputT>
  void grassfire( ImageViewBase<SourceT> const& src, ImageView<OutputT>& dst,
                  bool ignore_borders=false);

  // Without destination given, return in a newly-created ImageView<int32>
  template <class SourceT>
  ImageView<int32> grassfire( ImageViewBase<SourceT> const& src, bool ignore_borders=false) {
    int32 cols = src.impl().cols(), rows = src.impl().rows();
    ImageView<int32> result( cols, rows );
    grassfire( src, result, ignore_borders );
    return result;
  }

  /// A weight at a given pixel, based on an image row. Return
  /// zero where image values are not valid, and positive where valid.
  /// - hCenterLine contains the center column at each row/col
  /// - hMaxDistArray contains the width of the column at each row/col
  double compute_line_weights(Vector2 const& pix, bool horizontal,
                              std::vector<double> const& centers,
                              std::vector<double> const& widths);

  // *******************************************************************
  // centerline_weights()
  // *******************************************************************

  /// Computes a weighting measure based on the distance from the vertical
  ///  and horizontal centerlines of and image.
  /// - For images which are regular with no large holes this can work
  ///   better than using grassfire weights.
  /// - fill_holes will assign a normal weight to holes in the image interior.
  /// - use_min_weight creates more of a rectangular instead of a
  ///   circular weight pattern.
  template <class ImageT>
  void centerline_weights( ImageT const& src, ImageView<double>& weights,
                           BBox2i roi=BBox2i(), bool fill_holes = false,
                           bool use_min_weight=false );


  // *******************************************************************
  // bounding_box()
  // *******************************************************************

  template <class ViewT>
  BBox2i bounding_box( ImageViewBase<ViewT> const& image_ ) {
    return BBox2i( 0, 0, image_.impl().cols(), image_.impl().rows() );
  }


  // *******************************************************************
  // nonzero_data_bounding_box()
  // *******************************************************************

  template <class ViewT>
  BBox2i nonzero_data_bounding_box( ImageViewBase<ViewT> const& image_ );


  // *******************************************************************
  // is_opaque()
  // *******************************************************************

  template <class ImageT>
  bool is_opaque_helper( ImageT const& image, true_type ) {
    for( int32 y=0; y<image.rows(); ++y )
      for( int32 x=0; x<image.cols(); ++x )
        if( ! (is_opaque( image(x,y) ) ) )
          return false;
    return true;
  }

  template <class ImageT>
  bool is_opaque_helper( ImageT const& /*image*/, false_type ) {
    return true;
  }

  /// Returns true if the given image is entirely opaque, or false if
  /// it is at least partially transparent.
  template <class ImageT>
  bool is_opaque( ImageViewBase<ImageT> const& image ) {
    return is_opaque_helper( image.impl(), typename PixelHasAlpha<typename ImageT::pixel_type>::type() );
  }


  // *******************************************************************
  // is_transparent()
  // *******************************************************************

  template <class ImageT>
  bool is_transparent_helper( ImageT const& image, true_type ) {
    for( int32 y=0; y<image.rows(); ++y )
      for( int32 x=0; x<image.cols(); ++x )
        if( ! is_transparent(image(x,y)) ) return false;
    return true;
  }

  template <class ImageT>
  bool is_transparent_helper( ImageT const& /*image*/, false_type ) {
    return false;
  }

  /// Returns true if the given image is entirely transparent, or false if
  /// it is opaque or only partially transparent.
  template <class ImageT>
  bool is_transparent( ImageViewBase<ImageT> const& image ) {
    return is_transparent_helper( image.impl(), typename PixelHasAlpha<typename ImageT::pixel_type>::type() );
  }

  // *******************************************************************
  // subdivide_bbox()
  // *******************************************************************

  // TODO: Move this to BBox.h!

  /// A utility routine that, given an image, returns a vector of
  /// bounding boxes for sub-regions of the image of the specified
  /// size.  Note that bounding boxes along the right and bottom edges
  /// of the image will not have the specified dimension unless the
  /// image width and height are perfectly divisible by the bounding
  /// box width and height, respectively. This routine is useful if you
  /// want to apply an operation to a large image one region at a time.
  /// It will operate on any object that has cols() and rows() methods.
  /// - If include_partials is set to false, only full size boxes will be 
  ///   included in the output.
  /// - If full_size is true, the boxes at the right and bottom edges will
  ///   be grown inward to have full size. This will result in them overlapping
  ///   with earlier boxes. This assumes include_partials = true.
  /// - Output tiles are in raster order, top left to bottom right.
  std::vector<BBox2i>
  subdivide_bbox(BBox2i const& object, int32 block_width, int32 block_height,
                 bool include_partials = true, bool full_size = false);

  template <class T>
  inline std::vector<BBox2i>
  subdivide_bbox(ImageViewBase<T> const& view, int32 block_width, int32 block_height,
                 bool include_partials = true) {
    return subdivide_bbox( bounding_box(view.impl()), block_width, block_height, include_partials );
  }
}

#include "AlgorithmFunctions.tcc"

#endif//__VW_IMAGE_ALGORITHMS_FUNCTIONS_H__
