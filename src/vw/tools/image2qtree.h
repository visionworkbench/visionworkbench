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


#ifndef __VW_TOOLS_IMAGE2QTREE_H__
#define __VW_TOOLS_IMAGE2QTREE_H__

#include <vw/Core/ProgressCallback.h>
#include <vw/Core/Log.h>
#include <vw/Math/BBox.h>
#include <vw/Math/Vector.h>
#include <vw/Image/Transform.h>
#include <vw/Image/MaskViews.h>
#include <vw/FileIO/DiskImageResource.h>
#include <vw/FileIO/DiskImageResourceJPEG.h>
#include <vw/FileIO/DiskImageResourcePNG.h>
#include <vw/FileIO/DiskImageView.h>
#include <vw/Cartography/GeoReference.h>
#include <vw/Cartography/GeoTransform.h>
#include <vw/Mosaic/CelestiaQuadTreeConfig.h>
#include <vw/Mosaic/GigapanQuadTreeConfig.h>
#include <vw/Mosaic/ImageComposite.h>
#include <vw/Mosaic/KMLQuadTreeConfig.h>
#include <vw/Mosaic/QuadTreeConfig.h>
#include <vw/Mosaic/QuadTreeGenerator.h>
#include <vw/Mosaic/UniviewQuadTreeConfig.h>
#include <vw/tools/Common.h>
#include <vw/FileIO/GdalWriteOptions.h>

#include <boost/algorithm/string/trim.hpp>
#include <boost/smart_ptr/shared_ptr.hpp>
#include <boost/foreach.hpp>

#include <fstream>

namespace fs = boost::filesystem;

template <class ImageT, class TransformT>
vw::TransformView<ImageT, TransformT>
transform_only( vw::ImageViewBase<ImageT> const& v,
                TransformT const& transform_func ) {
  return vw::TransformView<ImageT, TransformT>(v.impl(),transform_func);
}

struct Options: vw::GdalWriteOptions {

  Options() :
    output_file_type(""),
    module_name(""),
    nudge_x(0), nudge_y(0),
    tile_size(0),
    jpeg_quality(-9999),
    png_compression(99999),
    pixel_scale(0),
    pixel_offset(0),
    aspect_ratio(1),
    global_resolution(0),
    nodata(0),
    nodata_set(false),
    north(0), south(0),
    east(0), west(0),
    channel_type("DEFAULT"),
    multiband(false),
    help(false),
    normalize(false),
    terrain(false),
    manual(false),
    global(false)
  {}

  std::vector<std::string> input_files;

  std::string output_file_name;
  std::string output_file_type;
  std::string module_name;
  double      nudge_x, nudge_y;
  vw::uint32  tile_size;
  float       jpeg_quality;
  vw::uint32  png_compression;
  float       pixel_scale, pixel_offset;
  vw::int32   aspect_ratio;
  vw::uint32  global_resolution;
  float       nodata;
  bool        nodata_set;
  float       north, south, east, west;

  std::string channel_type;
  std::string mode; // Quadtree type

  bool multiband;
  bool help;
  bool normalize;
  bool terrain;
  bool manual;
  bool global;

  struct {
    vw::uint32 draw_order_offset;
    vw::uint32 max_lod_pixels;
  } kml;

  struct proj_{
    std::string type;
    double    lat, lon, scale /*=1*/;
    double    p1, p2;
    vw::int32 utm_zone;
    proj_() :
      type("DEFAULT"),
      scale(1),
      utm_zone(0) {}
  } proj;

  struct datum_ {
    std::string type;
    float sphere_radius;
    datum_() : type("NONE"), sphere_radius(0) {}
  } datum;

  void validate() {
    VW_ASSERT(!help, vw::tools::Usage());
    VW_ASSERT(input_files.size() > 0,
              vw::tools::Usage() << "Need at least one input image");

    if (datum.type == "SPHERE")
      VW_ASSERT(datum.sphere_radius > 0,
                vw::tools::Usage() << "Sphere datum override requires a radius");

    if(output_file_name.empty())
      output_file_name = fs::path(input_files[0]).replace_extension().string();

    // Handle options for manually specifying projection bounds
    if (global || north!=0 || south!=0 || east!=0 || west!=0) {
      VW_ASSERT(input_files.size() == 1,
                vw::tools::Usage() << "Cannot override georeference information on multiple images");
      VW_ASSERT(global || ((north!=south) && (east!=west)),
                vw::tools::Usage() << "If you provide one, you must provide all of: --north --south --east --west");
      if (global) {
        north = 90; south = -90; east = 180; west = -180;
      }
      manual = true;
    }

    if (mode == "NONE")
      VW_ASSERT(input_files.size() == 1,
                vw::tools::Usage() << "Non-georeferenced images cannot be composed");
    if (mode == "CELESTIA" || mode == "UNIVIEW")
      VW_ASSERT(!module_name.empty(),
                vw::tools::Usage() << "Uniview and Celestia require --module-name");

    if (proj.type == "NONE")
      VW_ASSERT(input_files.size() == 1,
                vw::tools::Usage() << "Non-georeferenced images cannot be composed");
    // Compare against flag values
    if (jpeg_quality > 0)
      vw::DiskImageResourceJPEG::set_default_quality( jpeg_quality );
    if (png_compression!=99999)
      vw::DiskImageResourcePNG::set_default_compression_level( png_compression );
  }
};

// For image stretching.
static float lo_value = vw::ScalarTypeLimits<float>::highest();
static float hi_value = vw::ScalarTypeLimits<float>::lowest();

namespace vw {
  
template <class PixelT>
void do_normal_mosaic(const Options& opt, const vw::ProgressCallback *progress) {
  DiskImageView<PixelT> img(opt.input_files[0]);
  mosaic::QuadTreeGenerator quadtree(img, opt.output_file_name);
  quadtree.set_tile_size( 256 );
  quadtree.set_file_type( "png" );

  if ( opt.mode != "NONE" ) {
    boost::shared_ptr<mosaic::QuadTreeConfig> config = mosaic::QuadTreeConfig::make(opt.mode);
    config->configure( quadtree );
  }

  vw_out() << "Generating overlay..." << std::endl;
  vw_out() << "Writing: " << opt.output_file_name << std::endl;

  quadtree.generate( *progress );
}

/// Set up the input georeference object from the file or user inputs
vw::cartography::GeoReference
make_input_georef(boost::shared_ptr<vw::DiskImageResource> file,
                  const Options& opt);

/// Load the georeference for each input image and compute the maximum resolution
std::vector<vw::cartography::GeoReference>
load_image_georeferences( const Options& opt, int& total_resolution );

template <class PixelT>
void do_mosaic(const Options& opt, const vw::ProgressCallback *progress) {

  typedef typename PixelChannelType<PixelT>::type ChannelT;

  // If we're not outputting any special sort of mosaic (just a regular old
  // quadtree, no georeferencing, no metadata), we use a different function.
  if(opt.mode == "NONE" || opt.proj.type == "NONE") {
    do_normal_mosaic<PixelT>(opt, progress);
    return;
  }

  // Read in georeference info and compute total resolution.
  int total_resolution = 1024;
  std::vector<vw::cartography::GeoReference> georeferences = load_image_georeferences( opt, total_resolution );

  boost::shared_ptr<mosaic::QuadTreeConfig> config = mosaic::QuadTreeConfig::make(opt.mode);

  // Now that we have the best resolution, we can get our output_georef.
  int xresolution = total_resolution / opt.aspect_ratio, yresolution = total_resolution;

  vw::cartography::GeoReference output_georef = config->output_georef(xresolution, yresolution);
  vw_out(VerboseDebugMessage, "tool") << "Output Georef:\n" << output_georef << std::endl;

  // Configure the composite.
  mosaic::ImageComposite<PixelT> composite;

  // Add the transformed image files to the composite.
  for(size_t i=0; i < opt.input_files.size(); i++) {
    // Get info for this file
    const std::string & filename     = opt.input_files[i];
    const vw::cartography::GeoReference& input_georef = georeferences[i];

    // Load the image and georef from the file
    boost::shared_ptr<DiskImageResource> file( DiskImageResource::open(filename) );
    vw::cartography::GeoTransform geotx( input_georef, output_georef );

    ImageViewRef<PixelT> source = DiskImageView<PixelT>( file );

    // Handle nodata values/mask
    if ( opt.nodata_set ) {
      vw_out() << "Using nodata value: " << opt.nodata << "\n";
      source = mask_to_alpha(create_mask(pixel_cast<typename PixelWithoutAlpha<PixelT>::type >(source),ChannelT(opt.nodata)));
    } else if ( file->has_nodata_read() ) {
      vw_out() << "Using nodata value: " << file->nodata_read() << "\n";
      source = mask_to_alpha(create_mask(pixel_cast<typename PixelWithoutAlpha<PixelT>::type >(source),ChannelT(file->nodata_read())));
    }

    bool global = ((boost::trim_copy(input_georef.proj4_str())=="+proj=longlat") &&
      (fabs(input_georef.lonlat_to_pixel(Vector2(-180,  0)).x()                ) < 1) &&
      (fabs(input_georef.lonlat_to_pixel(Vector2( 180,  0)).x() - source.cols()) < 1) &&
      (fabs(input_georef.lonlat_to_pixel(Vector2(   0, 90)).y()                ) < 1) &&
      (fabs(input_georef.lonlat_to_pixel(Vector2(   0,-90)).y() - source.rows()) < 1));

    // Do various modifications to the input image here.
    if( (opt.pixel_scale !=1.0) || (opt.pixel_offset != 0.0) ) {
      vw_out() << "Apply input scaling: " << opt.pixel_scale << " offset: "
                                          << opt.pixel_offset << "\n";
      source = channel_cast<ChannelT>( source * opt.pixel_scale + opt.pixel_offset );
    }

    // Normalize pixel intensity if desired
    if( opt.normalize ) {
      vw_out() << "Apply normalizing: [" << lo_value << ", " << hi_value << "]\n";
      typedef ChannelRange<ChannelT> range_type;
      source = normalize_retain_alpha(source, lo_value, hi_value,
                                      range_type::min(), range_type::max());
    }

    BBox2i bbox = geotx.forward_bbox( BBox2i(0,0,source.cols(),source.rows()) );
    if ( global ) {
      vw_out() << "\t--> Detected global overlay. Using cylindrical edge extension to hide the seam.\n";
      source = crop( transform( source, geotx, source.cols(), source.rows(), CylindricalEdgeExtension() ), bbox );
    } else { // not global
      if ( norm_2(geotx.reverse(geotx.forward(Vector2()))) >
           0.01*norm_2(Vector2(source.cols(),source.rows())) ) {
        // Check for a fault were the forward bbox is correct, however
        // running a reverse through the geotransform projects 360
        // degrees off. Below seems like the only fix possible, as I
        // believe the problem happens because Proj4's fwd_pj will
        // always clamp to [-180,180].

        // However this fix will break in the event that the projection
        // doesn't loop back on itself. However if the projection did
        // that, I don't believe the test condition for this section
        // would be able to throw. This fix would also break if there
        // was a rotation in the georeference transform. GDAL however
        // doesn't support that.

        // For an example, see WAC global mosiac with tiles past the 180
        BBox2i correction(-geotx.reverse(geotx.forward(Vector2()))[0],0,source.cols(),source.rows());
        source = crop(transform_only( crop( interpolate(source), correction ),
                                      geotx ), bbox );
      } else {
        source = transform( source, geotx, bbox );
      }
    } // end if not global

    // Images that wrap the date line must be added to the composite on both sides.
    if( bbox.max().x() > total_resolution ) {
      composite.insert( source, bbox.min().x()-total_resolution, bbox.min().y() );
    }
    // Images that are in the 180-360 range *only* go on the other side.
    if( bbox.min().x() < xresolution ) {
      composite.insert( source, bbox.min().x(), bbox.min().y() );
    }
  } // End loop through input files

  // This box represents the entire input data set, in pixels, in the output
  // projection space. This should NOT include the extra data used to hide seams and such.
  BBox2i total_bbox = composite.bbox();
  total_bbox.crop(BBox2i(0,0,xresolution,yresolution));

  VW_ASSERT(total_bbox.width() > 0 && total_bbox.height() > 0,
            LogicErr() << "Total bbox is empty. Georeference calculation is probably incorrect.");

  if(opt.mode == "KML") {
    BBox2i bbox = total_bbox;
    // Compute a tighter Google Earth coordinate system aligned bounding box.
    int dim = 2 << (int)(log( (double)(std::max)(bbox.width(),bbox.height()) )/log(2.));
    if( dim > total_resolution )
      dim = total_resolution;
    total_bbox = BBox2i( (bbox.min().x()/dim)*dim, (bbox.min().y()/dim)*dim, dim, dim );
    if( ! total_bbox.contains( bbox ) ) {
      if( total_bbox.max().x() == xresolution )
        total_bbox.min().x() -= dim;
      else
        total_bbox.max().x() += dim;
      if( total_bbox.max().y() == yresolution )
        total_bbox.min().y() -= dim;
      else
        total_bbox.max().y() += dim;
    }
  }

  // Prepare the composite.
  if(!opt.multiband)
    composite.set_draft_mode( true );
  composite.prepare( total_bbox, *progress );
  VW_ASSERT(composite.rows() > 0 && composite.cols() > 0,
            LogicErr() << "Composite image is empty. Georeference calculation is probably incorrect.");

  mosaic::QuadTreeGenerator quadtree( composite, opt.output_file_name );

  // This whole bit here is terrible. This functionality should be moved into the Config base class somehow.
  if( opt.mode == "KML" ) {
    mosaic::KMLQuadTreeConfig *c2 = dynamic_cast<mosaic::KMLQuadTreeConfig*>(config.get());
    BBox2 ll_bbox( -180.0 + (360.0*total_bbox.min().x())/xresolution,
                    180.0 - (360.0*total_bbox.max().y())/yresolution,
                            (360.0*total_bbox.width())  /xresolution,
                            (360.0*total_bbox.height()) /yresolution );

    c2->set_longlat_bbox( ll_bbox );
    c2->set_max_lod_pixels( opt.kml.max_lod_pixels );
    c2->set_draw_order_offset( opt.kml.draw_order_offset );
  } else if( opt.mode == "CELESTIA" ) {
    mosaic::CelestiaQuadTreeConfig *c2 = dynamic_cast<mosaic::CelestiaQuadTreeConfig*>(config.get());
    c2->set_module(opt.module_name);
  } else if( opt.mode == "UNIVIEW" ) {
    mosaic::UniviewQuadTreeConfig *c2 = dynamic_cast<mosaic::UniviewQuadTreeConfig*>(config.get());
    c2->set_terrain(opt.terrain);
    c2->set_module(opt.module_name);
  } else if ( opt.mode == "GIGAPAN" ) {
    mosaic::GigapanQuadTreeConfig *c2 = dynamic_cast<mosaic::GigapanQuadTreeConfig*>(config.get());
    BBox2 ll_bbox( -180.0 + (360.0*total_bbox.min().x())/xresolution,
                    180.0 - (360.0*total_bbox.max().y())/yresolution,
                            (360.0*total_bbox.width())  /xresolution,
                            (360.0*total_bbox.height()) /yresolution );
    c2->set_longlat_bbox( ll_bbox );
  }

  config->configure(quadtree);

  if (opt.tile_size > 0)
    quadtree.set_tile_size(opt.tile_size);
  if (!opt.output_file_type.empty())
    quadtree.set_file_type(opt.output_file_type);

  // This box represents the input data, shifted such that total_bbox.min() is
  // the origin, and cropped to the size of the output resolution.
  BBox2i data_bbox = composite.bbox();
  data_bbox.crop( BBox2i(0,0,total_bbox.width(),total_bbox.height()));

  quadtree.set_crop_bbox(data_bbox);

  // Generate the composite.
  vw_out() << "Generating overlay..." << std::endl;
  vw_out() << "Writing: " << opt.output_file_name << std::endl;

  quadtree.generate(*progress);
}

// Define all of the function instantiations here, they are defined in
//  image2qtree_help.cc.
#define PROTOTYPE_ALL_CHANNEL_TYPES( PIXELTYPE )        \
  void do_mosaic_##PIXELTYPE##_uint8(const Options&, const vw::ProgressCallback*); \
  void do_mosaic_##PIXELTYPE##_int16(const Options&, const vw::ProgressCallback*); \
  void do_mosaic_##PIXELTYPE##_uint16(const Options&, const vw::ProgressCallback*); \
  void do_mosaic_##PIXELTYPE##_float32(const Options&, const vw::ProgressCallback*);

PROTOTYPE_ALL_CHANNEL_TYPES(PixelGrayA)
PROTOTYPE_ALL_CHANNEL_TYPES(PixelRGBA)

#undef PROTOTYPE_ALL_CHANNEL_TYPES

} // end namespace vw

#endif//__VW_TOOLS_IMAGE2QTREE_H__
