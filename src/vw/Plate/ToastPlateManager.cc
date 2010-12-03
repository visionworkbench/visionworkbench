// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <vw/Plate/ToastPlateManager.h>
#include <vw/Image/Filter.h>
#include <vw/Image/Algorithms.h>
#include <vw/Cartography/GeoReference.h>
#include <vw/Cartography/ToastTransform.h>

using namespace vw::platefile;
using namespace vw;

template <class PixelT>
void ToastPlateManager<PixelT>::transform_image(
                     cartography::GeoReference const& georef,
                     ImageViewRef<PixelT>& image,
                     TransformRef& txref, int& level ) const {
  // Compute the pyramid level at which to store this image.  The
  // number of required levels is broken down so that the very top
  // of the pyramid covers the entire globe and has a size of
  // tile_size.
  int tile_size = this->m_platefile->default_tile_size();
  Vector2 p0 = georef.pixel_to_lonlat(Vector2(image.cols()/2,
                                              image.rows()/2));
  Vector2 p1 = georef.pixel_to_lonlat(Vector2(image.cols()/2+1,
                                              image.rows()/2));
  Vector2 p2 = georef.pixel_to_lonlat(Vector2(image.cols()/2,
                                              image.rows()/2+1));
  double degrees_per_pixel = sqrt(pow(p1.y()-p0.y(),2)+pow(p2.y()-p0.y(),2));
  level = static_cast<int>(round(log(360/degrees_per_pixel/(tile_size-1)) / log(2)));

  // Compute the resolution of the TOAST output space at the given
  // pyramid_level.  The formula below was carefully chosen to
  // ensure that the output space at each level is sized just
  // slightly shy of an even power of two; therby allowing for the
  // slight overlap of one line of pixels on the right and the
  // bottom of a toast tile.  This overlap is necessary for proper
  // positioning when these tiles are rendered as texture in a 3D
  // graphics environment.
  int resolution = (1<<level)*(tile_size-1)+1;

  // Set up the toast transform and compute the bounding box of this
  // image in the toast projection space.
  cartography::ToastTransform toast_tx( georef, resolution );
  BBox2i output_bbox = toast_tx.forward_bbox(bounding_box(image));

  vw_out(InfoMessage, "platefile")
    << "\t    Placing image at level " << level
    << " with bbox " << output_bbox << "\n"
    << "\t    (Total TOAST resolution at this level =  "
    << resolution << " pixels.)\n";

  ImageViewRef<PixelT> holding =
    transform( image, toast_tx, ZeroEdgeExtension(),
               BicubicInterpolation() );
  image = holding;
  txref = TransformRef(toast_tx);
}

template <class PixelT>
cartography::GeoReference
ToastPlateManager<PixelT>::georeference( int /*level*/ ) const {
  vw_throw( NoImplErr() << "Toast Plate Manager can not provide a georeference object." );
}

template <class PixelT>
ImageView<PixelT>
ToastPlateManager<PixelT>::fetch_child_tile(int x, int y, int level,
                                            TransactionOrNeg transaction_id) const {
  int32 num_tiles = 1 << level;
  if ( x==-1 ) {
    if ( y==-1 )
      return fetch_child_tile(num_tiles-1, num_tiles-1, level, transaction_id);
    if ( y==num_tiles )
      return fetch_child_tile(num_tiles-1, 0, level, transaction_id);
    ImageView<PixelT> tile =
      fetch_child_tile(0, num_tiles-1-y, level, transaction_id);
    if ( tile ) return rotate_180(tile);
    else return tile;
  }
  if ( x==num_tiles ) {
    if ( y==-1 )
      return fetch_child_tile(0, num_tiles-1, level, transaction_id);
    if ( y==num_tiles )
      return fetch_child_tile(0, 0, level, transaction_id);
    ImageView<PixelT> tile =
      fetch_child_tile(num_tiles-1, num_tiles-1-y, level, transaction_id);
    if ( tile ) return rotate_180(tile);
    else return tile;
  }
  if ( y==-1 ) {
    ImageView<PixelT> tile =
      fetch_child_tile(num_tiles-1-x, 0, level, transaction_id);
    if( tile ) return rotate_180(tile);
    else return tile;
  }
  if ( y==num_tiles ) {
    ImageView<PixelT> tile =
      fetch_child_tile(num_tiles-1-x, num_tiles-1, level, transaction_id);
    if( tile ) return rotate_180(tile);
    else return tile;
  }

  ImageView<PixelT> tile;
  try {
    this->m_platefile->read(tile, x, y, level, transaction_id, true); // exact_transaction_match == true
  } catch (TileNotFoundErr &e) { /* No Op */ }

  // Regardless of what happened above, we return the tile here.  If
  // the read failed, the tile will be an empty image.
  return tile;
}

template <class PixelT>
void ToastPlateManager<PixelT>::generate_mipmap_tile(
    ImageView<PixelT>& dest, int col, int row, int level,
    TransactionOrNeg transaction_id, bool preblur) const
{
  // Create an image large enough to store all of the child nodes
  int tile_size = this->m_platefile->default_tile_size();
  ImageView<PixelT> super(4*tile_size-3, 4*tile_size-3);

  bool found = false;

  // Iterate over the children, gathering them and (recursively)
  // regenerating them if necessary.
  for ( int j=-1; j<3; ++j ) {
    for ( int i=-1; i<3; ++i ) {
      ImageView<PixelT> child = fetch_child_tile(2*col+i, 2*row+j, level+1, transaction_id);
      if (child) {
        crop(super,(tile_size-1)*(i+1),(tile_size-1)*(j+1), tile_size,tile_size) = child;
        found = true;
      }
    }
  }

  if (!found) {
    dest.reset();
    return;
  }

  // In the WWT implementation of TOAST the pixel centers
  // (rather than the than pixel corners) are grid-aligned, so
  // we need to use an odd-sized antialiasing kernel instead of
  // the usual 2x2 box filter.  The following 5-pixel kernel was
  // optimized to avoid the extra blurring associated with using
  // a kernel wider than 2 pixels.  Math was involved.
  std::vector<float> kernel(5);
  kernel[0] = kernel[4] = -0.0344f;
  kernel[1] = kernel[3] = 0.2135f;
  kernel[2] = 0.6418f;

  if (preblur)
    dest = subsample(
             crop(
               separable_convolution_filter(super, kernel, kernel, NoEdgeExtension() ),
               tile_size-1, tile_size-1, 2*tile_size, 2*tile_size ), 2);
  else
    dest = subsample( crop( super, tile_size-1, tile_size-1, 2*tile_size, 2*tile_size ), 2 );
}


// Explicit template instatiation
namespace vw {
namespace platefile {

#define VW_INSTANTIATE_TOAST_PLATEMANAGER_TYPES(PIXELT)                 \
  template void                                                         \
  ToastPlateManager<PIXELT >::transform_image(                          \
                               cartography::GeoReference const& georef, \
                               ImageViewRef<PIXELT >& image,            \
                               TransformRef& txref, int& level ) const; \
template void                                                           \
ToastPlateManager<PIXELT >::generate_mipmap_tile(                       \
    ImageView<PIXELT>& dest, int col, int row, int level,               \
    TransactionOrNeg transaction_id, bool preblur) const;               \
                                                                        \
  template cartography::GeoReference                                    \
  ToastPlateManager<PIXELT >::georeference( int level ) const;          \
  template ImageView<PIXELT >                                           \
  ToastPlateManager<PIXELT >::fetch_child_tile(                         \
      int col, int row, int level, TransactionOrNeg transaction_id) const; \

  VW_INSTANTIATE_TOAST_PLATEMANAGER_TYPES(PixelGrayA<uint8>)
  VW_INSTANTIATE_TOAST_PLATEMANAGER_TYPES(PixelGrayA<int16>)
  VW_INSTANTIATE_TOAST_PLATEMANAGER_TYPES(PixelGrayA<float32>)
  VW_INSTANTIATE_TOAST_PLATEMANAGER_TYPES(PixelRGBA<uint8>)
}}
