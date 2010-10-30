// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <vw/gui/TileGenerator.h>
#include <vw/Plate/HTTPUtils.h>
#include <vw/Image/ViewImageResource.h>
#include <vw/Image/ImageResourceView.h>
#include <vw/Image/Statistics.h>
#include <vw/Image/MaskViews.h>
#include <vw/Core/Debugging.h>

#include <QBuffer>
#include <QImage>
#include <QHttp>
#include <QByteArray>

using namespace vw;
using namespace vw::gui;
using namespace vw::platefile;

#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;

// --------------------------------------------------------------------------------
//                                Utility Functions
// --------------------------------------------------------------------------------

BBox2i vw::gui::tile_to_bbox(Vector2i tile_size, int col, int row, int level, int max_level) {
  if (col < 0 || row < 0 || col >= (1 << max_level) || row >= (1 << max_level) ) {
    return BBox2i();
  } else {
    BBox2i result(tile_size[0]*col, tile_size[1]*row, tile_size[0], tile_size[1]);
    return result * (1 << (max_level - level));
  }
}

std::list<TileLocator> vw::gui::bbox_to_tiles(Vector2i tile_size, BBox2i bbox, int level, int max_level, int transaction_id, bool exact_transaction_id_match) {
  std::list<TileLocator> results;

  // Compute the bounding box at the current level.
  BBox2i level_bbox = bbox / (1 << (max_level - level));

  // Grow that bounding box to align with tile boundaries
  BBox2i aligned_level_bbox = level_bbox;
  aligned_level_bbox.min().x() = ( (level_bbox.min().x() / tile_size[0]) * tile_size[0] );
  aligned_level_bbox.min().y() = ( (level_bbox.min().y() / tile_size[1]) * tile_size[1] );
  aligned_level_bbox.max().x() = ( int(ceilf( float(level_bbox.max().x()) / float(tile_size[0]) ))
                                   * tile_size[0] );
  aligned_level_bbox.max().y() = ( int(ceilf( float(level_bbox.max().y()) / float(tile_size[1]) ))
                                   * tile_size[1] );

  int tile_y = aligned_level_bbox.min().y() / tile_size[1];
  int dest_row = 0;
  while ( tile_y < aligned_level_bbox.max().y() / tile_size[1] ) {

    int tile_x = aligned_level_bbox.min().x() / tile_size[0];
    int dest_col = 0;
    while ( tile_x < aligned_level_bbox.max().x() / tile_size[0] ) {
      BBox2i tile_bbox(dest_col, dest_row, tile_size[0], tile_size[1]);
      TileLocator loc;
      loc.col = tile_x;
      loc.row = tile_y;
      loc.level = level;
      loc.transaction_id = transaction_id;
      loc.exact_transaction_id_match = exact_transaction_id_match;
      results.push_back(loc);

      ++tile_x;
      dest_col += tile_size[0];
    }
    ++tile_y;
    dest_row += tile_size[1];
  }
  return results;
}



// --------------------------------------------------------------------------
//                              TILE GENERATOR
// --------------------------------------------------------------------------

boost::shared_ptr<TileGenerator> TileGenerator::create(std::string filename_) {

  // Remove trailing /
  boost::trim_right_if(filename_, boost::is_any_of("/"));
  Url u(filename_);

  try {
    if (u.scheme() == "http") {
      return boost::shared_ptr<TileGenerator>( new WebTileGenerator(u.string(),17));
    } else if (u.scheme() == "file") {
      if (fs::extension(u.path()) == ".plate")
        return boost::shared_ptr<TileGenerator>( new PlatefileTileGenerator(u.path()) );
      else if (u.path() == "testpattern")
        return boost::shared_ptr<TileGenerator>( new TestPatternTileGenerator(256) );
      else
        return boost::shared_ptr<TileGenerator>( new ImageTileGenerator(u.path()) );
    } else {
      std::cerr << "Could not open " << u << ":\n\t" << "No handler for url scheme " << u.scheme() << std::endl;
    }
  } catch (const vw::Exception& e) {
    std::cerr << "Could not open " << u << ":\n\t" << e.what() << std::endl;
  }
  exit(EXIT_FAILURE);
}


// --------------------------------------------------------------------------
//                         TEST PATTERN TILE GENERATOR
// --------------------------------------------------------------------------

boost::shared_ptr<ViewImageResource> TestPatternTileGenerator::generate_tile(TileLocator const& /*tile_info*/) {
  ImageView<PixelRGBA<uint8> > tile(m_tile_size, m_tile_size);
  for (int j = 0; j < m_tile_size; ++j){
    for (int i = 0; i < m_tile_size; ++i){
      if (abs(i - j) < 10 || abs(i - (m_tile_size - j)) < 10)
        tile(i,j) = PixelRGBA<uint8>(255,0,0,255);
      else
        tile(i,j) = PixelRGBA<uint8>(0,0,0,255);
    }
  }
  boost::shared_ptr<ViewImageResource> result( new ViewImageResource(tile) );
  return result;
}

Vector2 TestPatternTileGenerator::minmax() { return Vector2(0.0, 1.0); }

PixelRGBA<float32> TestPatternTileGenerator::sample(int /*x*/, int /*y*/, int /*level*/, int /*transaction_id*/) {
  PixelRGBA<float32> result;
  return result;
}

int TestPatternTileGenerator::cols() const { return 2048; }
int TestPatternTileGenerator::rows() const { return 2048; }
PixelFormatEnum TestPatternTileGenerator::pixel_format() const { return VW_PIXEL_RGBA; }
ChannelTypeEnum TestPatternTileGenerator::channel_type() const { return VW_CHANNEL_UINT8; }
Vector2i TestPatternTileGenerator::tile_size() const {
  return Vector2i(m_tile_size, m_tile_size);
}
int32 TestPatternTileGenerator::num_levels() const {
  return 4;
}

// --------------------------------------------------------------------------
//                            WEB TILE GENERATOR
// --------------------------------------------------------------------------

HttpDownloadThread::HttpDownloadThread(const Url& u)
  : m_base_url(u), m_http(new QHttp(QString::fromStdString(m_base_url.hostname()), m_base_url.port() != 0 ? m_base_url.port() : 80, NULL))
{
  connect(m_http.get(), SIGNAL(requestFinished(int, bool)),this, SLOT(request_finished(int, bool)));
}

void HttpDownloadThread::run() {
  QThread::exec();
}

HttpDownloadThread::~HttpDownloadThread() { }

int HttpDownloadThread::get(const std::vector<std::string>& path, int transaction_id,
                            bool exact_transaction_id_match) {
  Url req_url(m_base_url);
  Url::split_t req_path = req_url.path_split();
  req_path.insert(req_path.end(), path.begin(), path.end());
  req_url.path_join(req_path);

  std::cout << "URL STRING: " << req_url << std::endl;

  // Set up request buffer
  RequestBufferPtr buf(new RequestBuffer(req_url));

  int request_id;
  req_url.query().set("nocache", "1");
  req_url.query().set("transaction_id", vw::stringify(transaction_id));

  if (exact_transaction_id_match)
    req_url.query().set("exact", "1");

  vw_out() << "\t --> Fetching " << req_url << "\n";

  request_id = m_http->get(QString::fromStdString(req_url.string()), &buf->buffer);

  Mutex::Lock lock(m_mutex);
  m_requests[request_id] = buf;

  return request_id;
}

bool HttpDownloadThread::result_available(int request_id) {
  Mutex::Lock lock(m_mutex);
  return m_requests[request_id]->ready;
}

vw::ImageView<vw::PixelRGBA<float> > HttpDownloadThread::pop_result(int request_id) {
  Mutex::Lock lock(m_mutex);
  map_t::iterator i = m_requests.find(request_id);
  VW_ASSERT(i != m_requests.end(), LogicErr() << "Asked for a request that doesn't exist!");
  RequestBufferPtr r = i->second;
  m_requests.erase(i);
  return r->result;
}

void HttpDownloadThread::request_finished(int request_id, bool error) {
  Mutex::Lock lock(m_mutex);

  map_t::iterator i = m_requests.find(request_id);
  // setHost and things count as requests
  if (i == m_requests.end())
    return;

  QHttpResponseHeader hdr = m_http->lastResponse();
  RequestBufferPtr buf = i->second;

  typedef PixelRGBA<float> pixel_t;

  switch (hdr.statusCode()) {
    case 200:
      // handled below
      break;
    case 404:
      buf->result.set_size(1,1);
      buf->result(0,0) = pixel_t(0.0f,0.1f,0.0f,1.0f);
      buf->ready = true;
      return;
    default:
      error = true;
      vw_out(WarningMessage, "console") <<
        "Request " << request_id << " failed with status " << hdr.statusCode() << " for URL: " << buf->url << "\n";
      break;
    }

    std::string filetype;
    if (hdr.contentType() == "image/png")
      filetype = "png";
    else if (hdr.contentType() == "image/jpeg" || hdr.contentType() == "image/jpg")
      filetype = "jpg";
    else if (hdr.contentType() == "image/tiff" || hdr.contentType() == "image/tif")
      filetype = "tif";
    else {
      vw_out(WarningMessage, "console") << "unrecognized content-type: " << hdr.contentType().toStdString() << "\n";
      error = true;
    }

    if (error) {
      buf->result.set_size(1,1);
      buf->result(0,0) = pixel_t(1.0f,0.0f,0.0f,1.0f);
      buf->ready = true;
      return;
    }

    const QByteArray& bytes = buf->buffer.data();

    boost::scoped_ptr<SrcImageResource> rsrc(
        SrcMemoryImageResource::open(
          filetype, reinterpret_cast<const uint8*>(bytes.constData()), bytes.size()));

    try {
      if (rsrc->channel_type() == VW_CHANNEL_UINT8) {
        buf->result = channel_cast_rescale<float>(ImageView<PixelRGBA<uint8> >(*rsrc));
      } else if (rsrc->channel_type() == VW_CHANNEL_UINT16) {
        buf->result = channel_cast_rescale<float>(ImageView<PixelRGBA<int16> >(*rsrc));
      } else {
        vw_out() << "WARNING: Image contains unsupported channel type: " << rsrc->channel_type() << "\n";
      }
      buf->ready = true;
    } catch (IOErr &e) {
      vw_out(WarningMessage) << "Could not parse network tile: " << buf->url << std::endl;
      buf->ready = true;
    }
}


WebTileGenerator::WebTileGenerator(const Url& url, int levels) :
  m_tile_size(256), m_levels(levels), m_download_thread(url) {
  m_download_thread.start();
}

boost::shared_ptr<ViewImageResource> WebTileGenerator::generate_tile(TileLocator const& tile_info) {

  std::vector<std::string> path(3);
  path.push_back(vw::stringify(tile_info.level));
  path.push_back(vw::stringify(tile_info.col));
  path.push_back(vw::stringify(tile_info.row) + ".png");

  int request_id = m_download_thread.get(path,
                                         tile_info.transaction_id,
                                         tile_info.exact_transaction_id_match);
  // TODO: busywait :(
  while(!m_download_thread.result_available(request_id))
    usleep(100);

  vw::ImageView<vw::PixelRGBA<float> > result = m_download_thread.pop_result(request_id);
  return boost::shared_ptr<ViewImageResource>( new ViewImageResource( result ) );
}

Vector2 WebTileGenerator::minmax() { return Vector2(0.0, 1.0); }

PixelRGBA<float32> WebTileGenerator::sample(int /*x*/, int /*y*/, int /*level*/, int /*transaction_id*/) {
  PixelRGBA<float32> result;
  return result;
}

int WebTileGenerator::cols() const { return m_tile_size * (1 << (m_levels-1)); }
int WebTileGenerator::rows() const { return m_tile_size * (1 << (m_levels-1)); }
PixelFormatEnum WebTileGenerator::pixel_format() const { return VW_PIXEL_RGBA; }
ChannelTypeEnum WebTileGenerator::channel_type() const { return VW_CHANNEL_UINT8; }
Vector2i WebTileGenerator::tile_size() const {
  return Vector2i(m_tile_size, m_tile_size);
}
int32 WebTileGenerator::num_levels() const {
  return m_levels;
}

// --------------------------------------------------------------------------
//                         PLATE FILE TILE GENERATOR
// --------------------------------------------------------------------------

PlatefileTileGenerator::PlatefileTileGenerator(std::string platefile_name) :
  m_platefile(new vw::platefile::PlateFile(platefile_name)) {
  m_num_levels = m_platefile->num_levels();
  std::cout << "\t--> Loading platefile \"" << platefile_name << "\" with "
            << m_num_levels << " levels.\n";
}


#define VW_DELEGATE_BY_PIXEL_TYPE(func, arg1, arg2)                                  \
  switch (this->pixel_format()) {                                                    \
  case VW_PIXEL_GRAY:                                                                \
  case VW_PIXEL_GRAYA:                                                               \
      if (this->channel_type() == VW_CHANNEL_UINT8) {                                \
        return func<PixelGrayA<uint8> >(arg1, arg2);                                 \
      } else if (this->channel_type() == VW_CHANNEL_INT16) {                         \
        return func<PixelGrayA<int16> >(arg1, arg2);                                 \
      } else if (this->channel_type() == VW_CHANNEL_FLOAT32) {                       \
        return func<PixelGrayA<float> >(arg1, arg2);                                 \
      } else {                                                                       \
        std::cout << "This platefile has a channel type that is not yet support by vwv.\n"; \
        std::cout << "Exiting...\n\n";                                               \
        exit(0);                                                                     \
      }                                                                              \
      break;                                                                         \
    case VW_PIXEL_RGB:                                                               \
    case VW_PIXEL_RGBA:                                                              \
      if (this->channel_type() == VW_CHANNEL_UINT8) {                                \
        return func<PixelRGBA<uint8> >(arg1, arg2);                                  \
      } else if (this->channel_type() == VW_CHANNEL_UINT16) {                        \
        return func<PixelRGBA<uint16> >(arg1, arg2);                                 \
      } else {                                                                       \
        std::cout << "This platefile has a channel type that is not yet support by vwv.\n"; \
        std::cout << "Exiting...\n\n";                                               \
        exit(0);                                                                     \
      }                                                                              \
      break;                                                                         \
    default:                                                                         \
      std::cout << "This platefile has a pixel format that is not yet support by vwv.\n"; \
      std::cout << "Exiting...\n\n";                                                 \
      exit(0);                                                                       \
    }                                                                                \

template<class PixelT>
Vector2 minmax_impl(TileLocator const& tile_info,
                    boost::shared_ptr<vw::platefile::PlateFile> platefile) {
  ImageView<PixelT> tile;
  platefile->read(tile, tile_info.col, tile_info.row, tile_info.level, tile_info.transaction_id);
  typename PixelChannelType<PixelT>::type min, max;
  min_max_channel_values(alpha_to_mask(tile), min, max);
  Vector2 result(min, max);
  std::cout << "Here is the original answer: " << result << "\n";
  result /= ChannelRange<typename PixelChannelType<PixelT>::type>::max();
  std::cout << "NEW MIN AND MAX: " << result << "\n";
  return result;
}

Vector2 PlatefileTileGenerator::minmax() {
  try {
    TileLocator loc;
    loc.col = 0;
    loc.row = 0;
    loc.level = 0;
    VW_DELEGATE_BY_PIXEL_TYPE(minmax_impl, loc, m_platefile)
  } catch (platefile::TileNotFoundErr &e) {
    return Vector2(0,1);
  }
}

// template<class PixelT>
// std::string sample_impl(TileLocator const& tile_info,
//                         Vector2 const& px_loc) {
//   ImageView<PixelT> tile;
//   platefile->read(tile, tile_info.col, tile_info.row, tile_info.level, tile_info.transaction_id);
//   VW_ASSERT(px_loc[0] >= 0 && px_loc[0] < tile.cols() &&
//             px_loc[1] >= 0 && px_loc[1] < tile.rows(),
//             ArgumentErr() << "sample_impl() invalid pixel location");
//   return tile(px_loc[0], px_loc[1]);
// }

PixelRGBA<float32> PlatefileTileGenerator::sample(int /*x*/, int /*y*/, int /*level*/, int /*transaction_id*/) {
  // TileLocator tile_loc;
  // tile_loc.col = floor(x/this->tile_size[0]);
  // tile_loc.row = floor(y/this->tile_size[1]);
  // tile_loc.level = this->num_levels();
  // px_loc = Vector2(x % this->tile_size[0],
  //                  y % this->tile_size[1]);

  try {
    return PixelRGBA<float32>(1.0, 0.0, 0.0, 1.0);
    //    VW_DELEGATE_BY_PIXEL_TYPE(sample_tile_impl, tile_loc, px_loc)
  } catch (platefile::TileNotFoundErr &e) {
    ImageView<PixelGrayA<uint8> > blank_tile(1,1);
    return PixelRGBA<float32>();
  }
}

template <class PixelT>
boost::shared_ptr<ViewImageResource> generate_tile_impl(TileLocator const& tile_info,
                                      boost::shared_ptr<vw::platefile::PlateFile> platefile) {
  ImageView<PixelT> tile(1,1);
  try {

    platefile->read(tile, tile_info.col, tile_info.row,
                    tile_info.level, tile_info.transaction_id,
                    tile_info.exact_transaction_id_match);

  } catch (platefile::TileNotFoundErr &e) {

    ImageView<PixelRGBA<uint8> > blank_tile(1,1);
    blank_tile(0,0) = PixelRGBA<uint8>(0, 20, 0, 255);
    return boost::shared_ptr<ViewImageResource>( new ViewImageResource(blank_tile) );

  } catch (vw::IOErr &e) {

    std::cout << "WARNING: AMQP ERROR -- " << e.what() << "\n";

  }
  return boost::shared_ptr<ViewImageResource>( new ViewImageResource(tile) );
}

boost::shared_ptr<ViewImageResource> PlatefileTileGenerator::generate_tile(TileLocator const& tile_info) {

  vw_out(DebugMessage, "gui") << "Request to generate platefile tile "
                              << tile_info.col << " " << tile_info.row
                              << " @ " << tile_info.level << "\n";

  VW_DELEGATE_BY_PIXEL_TYPE(generate_tile_impl, tile_info, m_platefile)

  // If we get to here, then there was no support for the pixel format.
  vw_throw(NoImplErr() << "Unsupported pixel format or channel type in TileGenerator.\n");
}

int PlatefileTileGenerator::cols() const {
  return this->tile_size()[0] * (1 << (m_num_levels-1));
}

int PlatefileTileGenerator::rows() const {
  return this->tile_size()[1] * (1 << (m_num_levels-1));
}

PixelFormatEnum PlatefileTileGenerator::pixel_format() const {
  return m_platefile->pixel_format();
}

ChannelTypeEnum PlatefileTileGenerator::channel_type() const {
  return m_platefile->channel_type();
}

Vector2i PlatefileTileGenerator::tile_size() const {
  return Vector2i(m_platefile->default_tile_size(),
                  m_platefile->default_tile_size());
}

int32 PlatefileTileGenerator::num_levels() const {
  return m_num_levels;
}


// --------------------------------------------------------------------------
//                             IMAGE TILE GENERATOR
// --------------------------------------------------------------------------

ImageTileGenerator::ImageTileGenerator(std::string filename) :
  m_filename(filename), m_rsrc( DiskImageResource::open(filename) ) {
  vw_out() << "\t--> Loading image: " << filename << ".\n";
}


// This little template makes the code below much cleaner.
template <class PixelT>
boost::shared_ptr<ViewImageResource> do_image_tilegen(boost::shared_ptr<SrcImageResource> rsrc,
                                                      BBox2i tile_bbox,
                                                      int level, int num_levels) {
  ImageView<PixelT> tile(tile_bbox.width(), tile_bbox.height());
  rsrc->read(tile.buffer(), tile_bbox);
  ImageView<PixelT> reduced_tile = subsample(tile, (1 << ((num_levels-1) - level)));
  return boost::shared_ptr<ViewImageResource>( new ViewImageResource(reduced_tile) );
}

boost::shared_ptr<ViewImageResource> ImageTileGenerator::generate_tile(TileLocator const& tile_info) {

  // Compute the bounding box of the image and the tile that is being
  // requested.  The bounding box of the tile depends on the pyramid
  // level we are looking at.
  BBox2i image_bbox(0,0,m_rsrc->cols(),m_rsrc->rows());
  BBox2i tile_bbox = tile_to_bbox(this->tile_size(), tile_info.col,
                                  tile_info.row, tile_info.level, this->num_levels());

  // Check to make sure the image intersects the bounding box.  Print
  // an error to screen and return an empty tile if it does not.
  if (!image_bbox.intersects(tile_bbox)) {
    vw_out() << "WARNING in ImageTileGenerator: a tile was requested that doesn't exist.";
    ImageView<PixelGray<uint8> > blank_tile(this->tile_size()[0], this->tile_size()[1]);
    return boost::shared_ptr<ViewImageResource>( new ViewImageResource(blank_tile) );
  }

  // Make sure we don't access any pixels outside the image boundary
  // by cropping the tile to the image dimensions.
  tile_bbox.crop(image_bbox);

  switch (this->pixel_format()) {
  case VW_PIXEL_GRAY:
    if (this->channel_type() == VW_CHANNEL_UINT8) {
      return do_image_tilegen<PixelGray<uint8> >(m_rsrc, tile_bbox,
                                                 tile_info.level, this->num_levels());
    } else if (this->channel_type() == VW_CHANNEL_INT16) {
      return do_image_tilegen<PixelGray<int16> >(m_rsrc, tile_bbox,
                                                 tile_info.level, this->num_levels());
    } else if (this->channel_type() == VW_CHANNEL_UINT16) {
      return do_image_tilegen<PixelGray<uint16> >(m_rsrc, tile_bbox,
                                                  tile_info.level, this->num_levels());
    } else if (this->channel_type() == VW_CHANNEL_FLOAT32) {
      return do_image_tilegen<PixelGray<float> >(m_rsrc, tile_bbox,
                                                  tile_info.level, this->num_levels());
    } else {
      std::cout << "This platefile has a channel type that is not yet support by vwv.\n";
      std::cout << "Exiting...\n\n";
      exit(0);
  }
  break;

  case VW_PIXEL_GRAYA:
    if (this->channel_type() == VW_CHANNEL_UINT8) {
      return do_image_tilegen<PixelGrayA<uint8> >(m_rsrc, tile_bbox,
                                                  tile_info.level, this->num_levels());
    } else if (this->channel_type() == VW_CHANNEL_INT16) {
      return do_image_tilegen<PixelGrayA<int16> >(m_rsrc, tile_bbox,
                                                  tile_info.level, this->num_levels());
    } else if (this->channel_type() == VW_CHANNEL_UINT16) {
      return do_image_tilegen<PixelGrayA<uint16> >(m_rsrc, tile_bbox,
                                                   tile_info.level, this->num_levels());
    } else if (this->channel_type() == VW_CHANNEL_FLOAT32) {
      return do_image_tilegen<PixelGrayA<float> >(m_rsrc, tile_bbox,
                                                  tile_info.level, this->num_levels());
    } else {
      std::cout << "This image has a channel type that is not yet support by vwv.\n";
      std::cout << "Exiting...\n\n";
      exit(0);
      }

    break;

  case VW_PIXEL_RGB:
    if (this->channel_type() == VW_CHANNEL_UINT8) {
      return do_image_tilegen<PixelRGB<uint8> >(m_rsrc, tile_bbox,
                                                tile_info.level, this->num_levels());
    } else if (this->channel_type() == VW_CHANNEL_UINT16) {
      return do_image_tilegen<PixelRGB<uint16> >(m_rsrc, tile_bbox,
                                                tile_info.level, this->num_levels());
    } else {
      std::cout << "This image has a channel type that is not yet support by vwv.\n";
      std::cout << "Exiting...\n\n";
      exit(0);
    }

    break;

  case VW_PIXEL_RGBA:
    if (this->channel_type() == VW_CHANNEL_UINT8) {
      return do_image_tilegen<PixelRGBA<uint8> >(m_rsrc, tile_bbox,
                                                 tile_info.level, this->num_levels());
    } else if (this->channel_type() == VW_CHANNEL_UINT16) {
      return do_image_tilegen<PixelRGBA<uint16> >(m_rsrc, tile_bbox,
                                                  tile_info.level, this->num_levels());
    } else {
      std::cout << "This image has a channel type that is not yet support by vwv.\n";
      std::cout << "Exiting...\n\n";
      exit(0);
    }

    break;

  default:
    std::cout << "This image has a pixel format that is not yet support by vwv.\n";
    std::cout << "Exiting...\n\n";
    exit(0);
  }

  vw_throw(NoImplErr() << "Unsupported pixel format or channel type in TileGenerator.\n");

}

Vector2 ImageTileGenerator::minmax() {
  vw_throw(NoImplErr() << VW_CURRENT_FUNCTION << " not implemented.");
}

PixelRGBA<float32> ImageTileGenerator::sample(int /*x*/, int /*y*/, int /*level*/, int /*transaction_id*/) {
  vw_throw(NoImplErr() << VW_CURRENT_FUNCTION << " not implemented.");
}


int ImageTileGenerator::cols() const {
  return m_rsrc->cols();
}

int ImageTileGenerator::rows() const {
  return m_rsrc->rows();
}

PixelFormatEnum ImageTileGenerator::pixel_format() const {
  return m_rsrc->pixel_format();
}

ChannelTypeEnum ImageTileGenerator::channel_type() const {
  return m_rsrc->channel_type();
}

Vector2i ImageTileGenerator::tile_size() const {
  return m_rsrc->block_read_size();
}

int32 ImageTileGenerator::num_levels() const {
  int32 max_dimension = std::max(this->cols(), this->rows());
  int32 max_tilesize = std::max(this->tile_size()[0], this->tile_size()[1]);
  return 1 + boost::numeric_cast<int32>(ceil(log(float(max_dimension) / float(max_tilesize)) / log(2)));
}
