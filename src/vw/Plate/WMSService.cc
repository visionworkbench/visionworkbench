// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <vw/Plate/WMSService.h>
#include <vw/Plate/PlateCarreePlateManager.h>
#include <vw/Plate/PlateView.h>
#include <vw/Plate/Exception.h>
#include <vw/Cartography/GeoTransform.h>
#include <vw/FileIO/DiskImageResource.h>
#include <vw/Math/Matrix.h>
#include <vw/Image/Filter.h>
#include <vw/Image/Algorithms.h>

#include <boost/regex.hpp>
#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;

#include <boost/foreach.hpp>

using namespace vw::platefile;

// ----------------------------------------------------------------------------------
//                                 PRIVATE METHODS
// ----------------------------------------------------------------------------------

namespace {
  std::vector<std::string> glob_plate_filenames(std::string const& root_directory)  {

    std::vector<std::string> result;

    if ( !fs::exists( root_directory ) )
      vw::vw_throw(vw::IOErr() << "Could not access the root directory: \"" << root_directory << "\"");

    // Create a regular expression for matching the pattern 'plate_<number>.blob"
    boost::regex re;
    re.assign(".*\\.plate", boost::regex_constants::icase);

    fs::directory_iterator end_itr; // default construction yields past-the-end

    // Iterate through the files in the platefile directory and return
    // any that match the regex above.
    for ( fs::directory_iterator itr( root_directory ); itr != end_itr; ++itr ) {
      if (boost::regex_match(itr->leaf(), re))
        result.push_back(itr->leaf());
    }

    return result;
  }
}

WMSServiceImpl::WMSServiceRecord WMSServiceImpl::add_platefile(std::string root_directory,
                                                               std::string plate_filename,
                                                               boost::shared_ptr<PlateFile> platefile) {
    // Build up an WMSServiceRecord
    WMSServiceRecord rec;
    rec.plate_name = plate_filename;
    rec.full_path = root_directory + "/" + plate_filename;
    rec.platefile = platefile;

    // Store the record in a std::map by platefile_id
    m_indices[platefile->index_header().platefile_id()] = rec;
    return rec;
}

/// Fetch an WMSServiceRecord for a given platefile_id, or throw an
/// exception if no record is found.  This is the first step in most
/// of the RPC handlers below, so the code is factored out here for
/// convenience.
WMSServiceImpl::WMSServiceRecord WMSServiceImpl::get_platefile(int platefile_id) {

  if (m_indices.find(platefile_id) == m_indices.end())
    vw_throw(InvalidPlatefileErr() << "No platefile matching this platefile_id could be found.");

  return m_indices[platefile_id];

}

// ----------------------------------------------------------------------------------
//                                 PUBLIC METHODS
// ----------------------------------------------------------------------------------

WMSServiceImpl::WMSServiceImpl(std::string root_directory, std::string cache_directory)
  : m_root_directory( fs::system_complete(root_directory).string() ),
    m_cache_directory( fs::system_complete(cache_directory).string() ) {

  // Search for all platefiles in the given root_directory.  A
  // platefile is any directory ending in *.plate.
  std::vector<std::string> platefiles = glob_plate_filenames(m_root_directory);
  if (platefiles.size() < 1)
    vw_out(InfoMessage, "plate:index_service") << "Warning: could not find any platefiles in the root directory.";

  BOOST_FOREACH( const std::string& name, platefiles ) {
    boost::shared_ptr<PlateFile> platefile(new PlateFile(m_root_directory + "/" + name));
    this->add_platefile(m_root_directory, name, platefile);
  }
}

void WMSServiceImpl::GetTile(::google::protobuf::RpcController* /*controller*/,
                     const ::vw::platefile::WMSTileRequest* request,
                     ::vw::platefile::WMSTileResponse* response,
                     ::google::protobuf::Closure* done) {

  response->set_filename(create_image(request));
  done->Run();
}

namespace {

  using namespace vw;
  using namespace vw::math;
  using namespace vw::cartography;

  GeoReference build_georef(const vw::BBox2& lonlat, const vw::BBox2i& pixels) {
      // Build up a GeoReference for the requested image
      // XXX: (hardwired for moon right now...)
      GeoReference request_georef;
      const double LUNAR_RADIUS = 1737400;
      request_georef.set_datum(Datum("D_MOON", "MOON", "Reference Meridian", LUNAR_RADIUS, LUNAR_RADIUS, 0.0));

      Matrix3x3 request_xform;
      request_xform.set_identity();
      request_xform(0,0) = lonlat.width() / pixels.width();
      request_xform(0,2) = lonlat.min().x();
      request_xform(1,1) = lonlat.height() / pixels.height();
      request_xform(1,2) = lonlat.min().y();
      request_georef.set_transform(request_xform);
      return request_georef;
  }

  // The max number of tiles that a request can pull in.
  const uint32 MAX_TILE_COVERAGE = 9;

  template <typename PixelT>
  std::string create_image_helper(const std::string& cache_dir, boost::shared_ptr<PlateFile> plate, const vw::BBox2& lonlat, const vw::BBox2i& pixels) {
    std::ostringstream fn;

    // XXX: PNG is hard-coded.
    fn << cache_dir
       << "/" << plate->index_header().platefile_id()
       << "/" << pixels.width() << "x" << pixels.height()
       << "/" << lonlat.min().y() << "," << lonlat.min().x() << ","
              << lonlat.max().y() << "," << lonlat.max().x() << ".png";

    const fs::path& filename(fn.str());

    if (fs::exists(filename)) {
      vw_out(DebugMessage, "plate.WMSService") << "Cache Hit: " << filename << std::endl;
      return filename.string();
    }

    GeoReference request_georef = build_georef(lonlat, pixels);

    PlateCarreePlateManager<PixelT> manager(plate);

    int32 tile_size = plate->default_tile_size();
    int32 level = plate->num_levels()-1;
    GeoTransform trans(manager.georeference(level), request_georef);
    // mapping from plate space (at given level) to request space

    // We're looking for a level that maximizes resolution without pulling more than MAX_TILE_COVERAGE tiles.
    while (kml_image_tiles(trans.reverse_bbox(pixels), tile_size).size() > MAX_TILE_COVERAGE) {
      trans = GeoTransform(manager.georeference(--level), request_georef);
    }

    create_directories(filename.branch_path());

    PlateView<PixelT> image(plate);
    image.set_level(level);
    write_image(filename.string(),
        crop(
          transform(image, trans, ConstantEdgeExtension(), BilinearInterpolation()),
          pixels));

    return filename.string();
  }

  BBox2 adapt_bbox(const BBoxContainer& in_box) {
    return BBox2(in_box.origin_x(), in_box.origin_y(), in_box.width(), in_box.height());
  }
}

std::string WMSServiceImpl::create_image(const WMSTileRequest* req) {

  //vw_out(DebugMessage) << "Request: " << req->DebugString() << std::endl;

  boost::shared_ptr<PlateFile> plate = get_platefile(req->platefile_id()).platefile;
  const std::string description = plate->index_header().description();

  BBox2 lonlat(adapt_bbox(req->lonlat())),
        pixels(adapt_bbox(req->pixels()));

  if (description.find("DEM") != std::string::npos) {
    // If the description has DEM in it, return it has a height map
    return create_image_helper<PixelGray<int16> >(m_cache_directory, plate, lonlat, pixels);
  } else {
    // Otherwise, treat it as imagery
    return create_image_helper<PixelRGBA<uint8> >(m_cache_directory, plate, lonlat, pixels);
  }
}
