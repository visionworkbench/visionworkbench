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


#ifndef __VW_TOOLS_GEOBLEND_H__
#define __VW_TOOLS_GEOBLEND_H__

#include <vw/tools/Common.h>
#include <vw/Mosaic/ImageComposite.h>
#include <vw/Cartography/GeoReference.h>
#include <vw/Cartography/GeoReferenceUtils.h>
#include <vw/FileIO/GdalWriteOptions.h>
#include <vw/Cartography/GeoTransform.h>
#include <vw/FileIO/DiskImageView.h>
#include <vw/Image/PerPixelViews.h>
#include <vw/Image/Filter.h>
#include <vw/FileIO/DiskImageResourceGDAL.h>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

// Global Variables from the command line.
extern std::vector<std::string> image_files;
extern std::string mosaic_name;
extern std::string output_file_type;
extern std::string channel_type_str;
extern bool draft;
extern unsigned int tilesize;
extern bool tile_output;
extern unsigned int patch_size, patch_overlap;
extern float nodata_value;
extern bool has_nodata_value;

namespace vw {

  // Creates an alpha channel based on pixels with a value of zero.  The
  // first version covers scalar types.  The remaining versions cover
  // compound pixel types.
  template <class PixelT> struct AlphaTypeFromPixelType { typedef PixelGrayA<PixelT> type; };
  template<class ChannelT> struct AlphaTypeFromPixelType<PixelGray<ChannelT> > { typedef PixelGrayA<ChannelT> type; };
  template<class ChannelT> struct AlphaTypeFromPixelType<PixelGrayA<ChannelT> > { typedef PixelGrayA<ChannelT> type; };
  template<class ChannelT> struct AlphaTypeFromPixelType<PixelRGB<ChannelT> > { typedef PixelRGBA<ChannelT> type; };
  template<class ChannelT> struct AlphaTypeFromPixelType<PixelRGBA<ChannelT> > { typedef PixelRGBA<ChannelT> type; };

  template <class PixelT> struct NonAlphaTypeFromPixelType { typedef PixelGrayA<PixelT> type; };
  template<class ChannelT> struct NonAlphaTypeFromPixelType<PixelGray<ChannelT> > { typedef PixelGray<ChannelT> type; };
  template<class ChannelT> struct NonAlphaTypeFromPixelType<PixelGrayA<ChannelT> > { typedef PixelGray<ChannelT> type; };
  template<class ChannelT> struct NonAlphaTypeFromPixelType<PixelRGB<ChannelT> > { typedef PixelRGB<ChannelT> type; };
  template<class ChannelT> struct NonAlphaTypeFromPixelType<PixelRGBA<ChannelT> > { typedef PixelRGB<ChannelT> type; };

  template <class PixelT>
  class NodataToMaskFunctor: public vw::UnaryReturnTemplateType<AlphaTypeFromPixelType> {
    PixelT m_nodata_value;

  public:
    NodataToMaskFunctor(typename PixelChannelType<PixelT>::type nodata_value = 0) : m_nodata_value(nodata_value) {}

    typename AlphaTypeFromPixelType<PixelT>::type operator() (PixelT const& pix) const {
      typedef typename AlphaTypeFromPixelType<PixelT>::type result_type;
      if (pix == m_nodata_value)
        return result_type();  // Mask pixel
      else
        return result_type(pix);
    }
  };

  template <class ViewT>
  vw::UnaryPerPixelView<ViewT, NodataToMaskFunctor<typename ViewT::pixel_type> >
  nodata_to_mask(vw::ImageViewBase<ViewT> const& view,
                 typename PixelChannelType<typename ViewT::pixel_type>::type const& nodata_value = 0 ) {
    return vw::per_pixel_filter(view.impl(), NodataToMaskFunctor<typename ViewT::pixel_type>(nodata_value));
  }

  template <class PixelT>
  class MaskToNodataFunctor: public vw::UnaryReturnTemplateType<NonAlphaTypeFromPixelType> {
    PixelT m_nodata_value;
    typedef typename PixelChannelType<PixelT>::type channel_type;

  public:
    MaskToNodataFunctor(float nodata_value = 0) : m_nodata_value((channel_type)nodata_value) {}

    typename NonAlphaTypeFromPixelType<PixelT>::type operator() (PixelT const& pix) const {
      typedef typename NonAlphaTypeFromPixelType<PixelT>::type result_type;
      if (is_transparent(pix))
        return result_type(m_nodata_value);
      else
        return result_type(pix);
    }
  };

  template <class ViewT>
  vw::UnaryPerPixelView<ViewT, MaskToNodataFunctor<typename ViewT::pixel_type> >
  mask_to_nodata(vw::ImageViewBase<ViewT> const& view, float nodata_value = 0 ) {
    return vw::per_pixel_filter(view.impl(), MaskToNodataFunctor<typename ViewT::pixel_type>(nodata_value));
  }

  // do_blend()
  //
  // This performs the actual work of geoblend
  template <class PixelT>
  void do_blend() {

    typedef typename AlphaTypeFromPixelType<PixelT>::type alpha_pixel_type;
    typedef typename PixelChannelCast<alpha_pixel_type,float32>::type float_pixel_type;

    TerminalProgressCallback tpc( "tools.geoblend", "");

    vw::mosaic::ImageComposite<float_pixel_type> composite;
    if( draft ) composite.set_draft_mode( true );

    double smallest_x_scale = vw::ScalarTypeLimits<float>::highest();
    double smallest_y_scale = vw::ScalarTypeLimits<float>::highest();
    double smallest_x_val = vw::ScalarTypeLimits<float>::highest();
    double largest_y_val = vw::ScalarTypeLimits<float>::lowest();

    // First pass, read georeferencing information and build an output
    // georef.
    tpc.set_progress_text( "Status (scanning):   " );
    SubProgressCallback scanning_pc( tpc, 0, 0.05 );
    for(unsigned i = 0; i < image_files.size(); ++i) {
      vw_out(vw::VerboseDebugMessage) << "Adding file " << image_files[i] << std::endl;
      scanning_pc.report_fractional_progress(i,image_files.size());

      cartography::GeoReference input_georef;
      read_georeference( input_georef, image_files[i] );
      DiskImageView<PixelT> source_disk_image( image_files[i] );
      vw_out(vw::VerboseDebugMessage) << "\tTransform: " << input_georef.transform()
                                      << "\t\tBBox: " << input_georef.bounding_box(source_disk_image) << std::endl;

      // Check to make sure the image has valid georeferencing
      // information.
      if( input_georef.transform() == identity_matrix<3>() ) {
        vw_out(InfoMessage) << "No georeferencing info found for image: \"" << image_files[i] << "\".  Aborting." << std::endl;
        exit(0);
      }

      Matrix3x3 affine = input_georef.transform();
      if (fabs(affine(0,0)) < smallest_x_scale || fabs(affine(1,1)) < smallest_y_scale) {
        smallest_x_scale = affine(0,0);
        smallest_y_scale = affine(1,1);
      }

      if (affine(0,2) < smallest_x_val)
        smallest_x_val = affine(0,2);

      // Note: since the y coordinates are typically flipped in a DEM,
      // we look here for the _largest_ y value, since it corresponds to
      // the upper left hand pixel.
      if (affine(1,2) > largest_y_val)
        largest_y_val = affine(1,2);

    }
    scanning_pc.report_finished();

    // Convert all of the images so they share the same scale factor.
    // Adopt the scale of the highest resolution being composited.
    Matrix3x3 output_affine = identity_matrix<3>();
    output_affine(0,0) = smallest_x_scale;
    output_affine(1,1) = smallest_y_scale;
    output_affine(0,2) = smallest_x_val;
    output_affine(1,2) = largest_y_val;

    vw_out(VerboseDebugMessage) << "Output affine transform: " << output_affine << std::endl
                                << int(smallest_x_val) << "  " << int(largest_y_val) << std::endl;

    // Take the georef from the first file (this ensures that the
    // projection and datum information is preserved...), but update the
    // affine transform.
    cartography::GeoReference output_georef;
    read_georeference( output_georef, image_files[0] );
    output_georef.set_transform(output_affine);

    tpc.set_progress_text( "Status (assembling): " );
    SubProgressCallback assembling_pc( tpc, 0.05, 0.1 );
    // Second pass: add files to the image composite.
    for(size_t i = 0; i < image_files.size(); ++i) {
      assembling_pc.report_fractional_progress(i, image_files.size() );
      cartography::GeoReference input_georef;
      read_georeference(input_georef, image_files[i]);
      DiskImageView<PixelT> source_disk_image( image_files[i] );

      cartography::GeoTransform trans(input_georef, output_georef);
      BBox2 output_bbox = trans.forward_bbox( BBox2(0,0,source_disk_image.cols(),source_disk_image.rows()) );
      vw_out(vw::VerboseDebugMessage) << "output_bbox = " << output_bbox << std::endl;

      // I've hardwired this to use nearest pixel interpolation for now
      // until we have a chance to sit down and develop a better
      // strategy for intepolating and filtering in the presence of
      // missing pixels in DEMs. -mbroxton
      if (has_nodata_value) {
        ImageViewRef<alpha_pixel_type> masked_source = crop( transform( nodata_to_mask(source_disk_image, (typename PixelChannelType<PixelT>::type)(nodata_value) ), trans, ZeroEdgeExtension(), NearestPixelInterpolation() ), output_bbox );
        composite.insert( channel_cast_rescale<float32>(masked_source), (int)output_bbox.min().x(), (int)output_bbox.min().y() );
      } else {
        ImageViewRef<alpha_pixel_type> masked_source = crop( transform( pixel_cast<alpha_pixel_type>(source_disk_image), trans, ZeroEdgeExtension(), NearestPixelInterpolation() ), output_bbox );
        composite.insert( channel_cast_rescale<float32>(masked_source), (int)output_bbox.min().x(), (int)output_bbox.min().y() );
      }

    }
    assembling_pc.report_finished();

    tpc.set_progress_text( "Status (preparing):  " );
    vw_out(vw::VerboseDebugMessage) << std::endl;
    composite.prepare( SubProgressCallback( tpc, 0.1, 0.5 ) );
    vw_out(vw::VerboseDebugMessage) << "Composite dimensions: " << composite.cols() << "  " << composite.rows() << std::endl;

    tpc.set_progress_text( "Status (blending):   " );
    SubProgressCallback blending_pc( tpc, 0.5, 1.0 );
    // Output the image in tiles, or one large image.
    if(tile_output) {
      const int dim = patch_size - patch_overlap;
      const int tile_width = composite.cols() / dim;
      const int tile_height = composite.rows() / dim;
      vw_out(vw::VerboseDebugMessage) << "Outputting composite in " << tile_width * tile_height << " tiles." << std::endl;

      for(int i=0; i < composite.rows(); i += dim) {
        for(int j=0; j < composite.cols(); j += dim) {
          BBox2i tile_bbox(j, i, dim, dim);
          if(tile_bbox.max().x() >= composite.cols())
            tile_bbox.max().x() = composite.cols();
          if(tile_bbox.max().y() >= composite.rows())
            tile_bbox.max().y() = composite.cols();
          ImageView<PixelT> tile_view =
            pixel_cast<PixelT>(crop(channel_cast_rescale<typename PixelChannelType<PixelT>::type>(composite), tile_bbox));
          cartography::GeoReference tile_georef = output_georef;

          // Adjust the affine transformation's offset to point to the upper
          // left of this tile.
          Vector2 upper_left = output_georef.pixel_to_point( Vector2(j, i) );
          Matrix3x3 tile_transform = tile_georef.transform();
          tile_transform(0,2) += upper_left[0];
          tile_transform(1,2) += upper_left[1];
          tile_georef.set_transform(tile_transform);

          // Filename for this tile.
          std::stringstream tile_filename;
          tile_filename << mosaic_name;
          tile_filename << '.' << j << '.' << i << '.';
          tile_filename << output_file_type;

          // Finally, write.
          write_gdal_image(tile_filename.str(), tile_view, tile_georef, GdalWriteOptions(), blending_pc);
        }
      }
    } else {
      vw_out(vw::VerboseDebugMessage) << "Output image:" << std::endl
                                      << "\tTransform: " << output_affine << std::endl
                                      << "\t\tBBox: " << output_georef.bounding_box(composite) << " [ W: " << output_georef.bounding_box(composite).width() << " H: " << output_georef.bounding_box(composite).height() << " ]" << std::endl << std::endl;

      std::string mosaic_filename = mosaic_name+".blend."+output_file_type;
      DiskImageResourceGDAL *out_resource;
      ImageViewRef<PixelT> out_image;

      // Specify the output image resource.
      if (has_nodata_value) {
        out_image = pixel_cast<PixelT>( mask_to_nodata( channel_cast_rescale<typename PixelChannelType<PixelT>::type>(composite), nodata_value ) );
      } else {
        out_image = pixel_cast<PixelT>( channel_cast_rescale<typename PixelChannelType<PixelT>::type>(composite) );
      }

      // Set up tiled TIFF output, if it's specified.
      if(tilesize > 0)
        out_resource = new DiskImageResourceGDAL( mosaic_filename, out_image.format(), Vector2i(tilesize, tilesize) );
      else
        out_resource = new DiskImageResourceGDAL( mosaic_filename, out_image.format() );

      // Finally, write.
      write_georeference(*out_resource, output_georef);
      write_image(*out_resource, out_image, blending_pc);
      delete out_resource;
    }
    blending_pc.report_finished();

    tpc.report_finished();
  }

} // end namespace vw

// These functions are defined by the geoblend_help.cc
#define PROTOTYPE_ALL_CHANNEL_TYPES( PIXELTYPE ) \
  void do_blend_##PIXELTYPE##_uint8();               \
  void do_blend_##PIXELTYPE##_int16();               \
  void do_blend_##PIXELTYPE##_uint16();              \
  void do_blend_##PIXELTYPE##_float32();

PROTOTYPE_ALL_CHANNEL_TYPES(PixelGray)
PROTOTYPE_ALL_CHANNEL_TYPES(PixelGrayA)
PROTOTYPE_ALL_CHANNEL_TYPES(PixelRGB)
PROTOTYPE_ALL_CHANNEL_TYPES(PixelRGBA)

#undef PROTOTYPE_ALL_CHANNEL_TYPES

#endif//__VW_TOOLS_GEOBLEND_H__
