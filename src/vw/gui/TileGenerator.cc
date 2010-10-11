// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <vw/gui/TileGenerator.h>
#include <vw/Image/ViewImageResource.h>
#include <vw/Image/ImageResourceView.h>
#include <vw/Image/Statistics.h>
#include <vw/Image/MaskViews.h>

#include <QUrl>
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
  if (col < 0 || row < 0 || col >= pow(2, max_level) || row >= pow(2, max_level) ) {
    return BBox2i();
  } else {
    BBox2i result(tile_size[0]*col, tile_size[1]*row, tile_size[0], tile_size[1]);
    return result * pow(2,max_level - level);
  }
}

std::list<TileLocator> vw::gui::bbox_to_tiles(Vector2i tile_size, BBox2i bbox, int level, int max_level, int transaction_id, bool exact_transaction_id_match) {
  std::list<TileLocator> results;

  // Compute the bounding box at the current level.
  BBox2i level_bbox = bbox / pow(2,max_level - level);

  // Grow that bounding box to align with tile boundaries
  BBox2i aligned_level_bbox = level_bbox;
  aligned_level_bbox.min().x() = ( (level_bbox.min().x() / tile_size[0]) * tile_size[0] );
  aligned_level_bbox.min().y() = ( (level_bbox.min().y() / tile_size[1]) * tile_size[1] );
  aligned_level_bbox.max().x() = ( int(ceilf( float(level_bbox.max().x()) / tile_size[0] ))
                                   * tile_size[0] );
  aligned_level_bbox.max().y() = ( int(ceilf( float(level_bbox.max().y()) / tile_size[1] ))
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

boost::shared_ptr<TileGenerator> TileGenerator::create(std::string filename) {

  // Strip off any trailing slashes to make sure we aren't
  // accidentally misparsing a platefile name.
  if (filename[filename.size()-1] == '/')
    filename.erase(filename.size()-1, 1);

  try {

    // If ends in .plate, then assume platefile.
    if ( fs::extension(filename) == ".plate") {

      return boost::shared_ptr<TileGenerator>( new PlatefileTileGenerator(filename) );

    // If begins with http://, then assume web tiles.
    } else if ( filename.find("http://") == 0) {

      return boost::shared_ptr<TileGenerator>( new WebTileGenerator(filename,17));

    // If testpattern, then we use the testpattern tile generator
    } else if (filename == "testpattern") {

      vw_out() << "\t--> Starting vwv in testpattern mode.\n";
      return boost::shared_ptr<TileGenerator>( new TestPatternTileGenerator(256) );

    // Otherwise, assume an image.
    } else {
      return boost::shared_ptr<TileGenerator>( new ImageTileGenerator(filename) );
    }
  } catch (vw::IOErr &e) {
    std::cout << "An error occurred opening \"" << filename << "\":\n\t" << e.what() << "\n";
    exit(0);
  }
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

HttpDownloadThread::HttpDownloadThread() {
  m_http = new QHttp(NULL);
  connect(m_http, SIGNAL(requestStarted(int)),this, SLOT(request_started(int)));
  connect(m_http, SIGNAL(requestFinished(int, bool)),this, SLOT(request_finished(int, bool)));
  connect(m_http, SIGNAL(responseHeaderReceived(QHttpResponseHeader)),this, SLOT(response_header_received(QHttpResponseHeader)));
}

void HttpDownloadThread::run() {
  QThread::exec();
}

HttpDownloadThread::~HttpDownloadThread() {
  if (m_http)
    delete m_http;
}

int HttpDownloadThread::get(std::string url_string, int transaction_id,
                            bool exact_transaction_id_match) {
  Mutex::Lock lock(m_mutex);
  QUrl url(url_string.c_str());

  std::cout << "URL STRING: " << url.toString().toStdString() << "\n";
  std::cout << "URL FRAG: " << url.fragment().toStdString() << "\n";

  /// XXX: Hard coding file_type as PNG for now. FIXME!!!!
  std::string file_type = "png";

  int request_id;
  if (m_http) {
    m_http->setHost(url.host());

    // Set up request buffer
    RequestBuffer buf;
    buf.file_type = file_type;
    buf.url = url_string;

    std::ostringstream path_with_opts;
    path_with_opts << url.path().toStdString() << "?nocache=1&transaction_id=" << transaction_id;
    if (exact_transaction_id_match)
      path_with_opts << "&exact=1";
    vw_out() << "\t --> Fetching " << path_with_opts.str() << "\n";
    QString final_path_str(path_with_opts.str().c_str());
    request_id = m_http->get (final_path_str, buf.buffer.get());
    m_requests[request_id] = buf;
    m_current_request = request_id;
  }
  return request_id;
}

bool HttpDownloadThread::result_available(int request_id) {
  Mutex::Lock lock(m_mutex);
  return m_requests[request_id].finished;
}

vw::ImageView<vw::PixelRGBA<float> > HttpDownloadThread::pop_result(int request_id) {
  Mutex::Lock lock(m_mutex);
  vw::ImageView<vw::PixelRGBA<float> > result = m_requests[request_id].result;
  m_requests.erase(request_id);
  return result;
}

void HttpDownloadThread::request_started(int /*request_id*/) {}

void HttpDownloadThread::response_header_received( const QHttpResponseHeader & resp ) {
  Mutex::Lock lock(m_mutex);

  std::map<int, RequestBuffer>::iterator request_iter = m_requests.find(m_current_request);
  if (request_iter != m_requests.end()) {
    RequestBuffer &buf = request_iter->second;
    buf.status = resp.statusCode();
    if (resp.contentType() == "image/png")
      buf.file_type = "png";
    else if (resp.contentType() == "image/jpeg" || resp.contentType() == "image/jpg")
      buf.file_type = "jpg";
    else if (resp.contentType() == "image/tiff" || resp.contentType() == "image/tif")
      buf.file_type = "tif";
    else if (resp.contentType() == "text/html") {
      /* do nothing... this is likely a 404 */
    } else
      vw_out() << "WARNING: unrecognized content-type: " << resp.contentType().toStdString() << "\n";
  }
}

void HttpDownloadThread::request_finished(int request_id, bool error) {
  Mutex::Lock lock(m_mutex);
  std::map<int, RequestBuffer>::iterator request_iter = m_requests.find(request_id);
  if (request_iter != m_requests.end()) {
    RequestBuffer &buf = request_iter->second;

    if (buf.status != 404 && buf.status != 200) {
      std::cout << "WARNING: Request " << request_id << " failed with status " << buf.status
                << " for URL: " << buf.url << "\n";
      ImageView<PixelRGBA<float> > vw_image(1,1);
      vw_image(0,0) = PixelRGBA<float>(1.0,0.0,0.0,1.0);
      buf.result = vw_image;
    } else if (buf.status == 404) {
      ImageView<PixelRGBA<float> > vw_image(1,1);
      vw_image(0,0) = PixelRGBA<float>(0.0,0.1,0.0,1.0);
      buf.result = vw_image;
    }

    // Handle the error case.
    if (error || buf.buffer->buffer().length() == 0 || buf.status != 200) {
      buf.finished = true;
      return;
    }

    // Save the data to a temporary file, and then read the image file
    // using the Vision Workbench FileIO subsystem.
    std::string temp_filename = platefile::TemporaryTileFile::unique_tempfile_name(buf.file_type);
    std::ofstream of(temp_filename.c_str());
    if ( !(of.good()) )
      vw_throw(IOErr() << "Could not open temporary tile file for writing: " << temp_filename);
    of.write(buf.buffer->buffer().data(), buf.buffer->buffer().length());
    of.close();
    TemporaryTileFile temp_tile_file(temp_filename);

    // Now read the image out of the tempfile and save the decoded
    // pixels as the result.
    try {
      boost::shared_ptr<DiskImageResource> rsrc(DiskImageResource::open(temp_filename));
      if (rsrc->channel_type() == VW_CHANNEL_UINT8) {
        ImageView<PixelRGBA<uint8> > vw_image = temp_tile_file.read<PixelRGBA<uint8> >();
        buf.result = channel_cast_rescale<float>(vw_image);
      } else if (rsrc->channel_type() == VW_CHANNEL_UINT16) {
        ImageView<PixelRGBA<int16> > vw_image = temp_tile_file.read<PixelRGBA<int16> >();
        buf.result = channel_cast_rescale<float>(vw_image);
      } else {
        vw_out() << "WARNING: Image contains unsupported channel type: "
                 << rsrc->channel_type() << "\n";
      }
      buf.finished = true;
    } catch (IOErr &e) {
      vw_out(WarningMessage) << "Could not read data from temporary file: "
                             << temp_filename << "\n";
      buf.finished = true;
    }
  }
}


WebTileGenerator::WebTileGenerator(std::string url, int levels) :
  m_tile_size(256), m_levels(levels), m_url(url) {
  m_download_thread.start();
}

boost::shared_ptr<ViewImageResource> WebTileGenerator::generate_tile(TileLocator const& tile_info) {

  std::ostringstream full_url;
  full_url << m_url << "/" << tile_info.level
           << "/" << tile_info.col
           << "/" << tile_info.row << ".png";
  int request_id = m_download_thread.get(full_url.str(),
                                         tile_info.transaction_id,
                                         tile_info.exact_transaction_id_match);
  while(!m_download_thread.result_available(request_id));

  vw::ImageView<vw::PixelRGBA<float> > result = m_download_thread.pop_result(request_id);
  return boost::shared_ptr<ViewImageResource>( new ViewImageResource( result ) );
}

Vector2 WebTileGenerator::minmax() { return Vector2(0.0, 1.0); }

PixelRGBA<float32> WebTileGenerator::sample(int /*x*/, int /*y*/, int /*level*/, int /*transaction_id*/) {
  PixelRGBA<float32> result;
  return result;
}

int WebTileGenerator::cols() const { return m_tile_size * pow(2,m_levels-1); }
int WebTileGenerator::rows() const { return m_tile_size * pow(2,m_levels-1); }
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
  return this->tile_size()[0] * pow(2, m_num_levels-1);
}

int PlatefileTileGenerator::rows() const {
  return this->tile_size()[1] * pow(2, m_num_levels-1);
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
boost::shared_ptr<ViewImageResource> do_image_tilegen(boost::shared_ptr<ImageResource> rsrc,
                                                      BBox2i tile_bbox,
                                                      int level, int num_levels) {
  ImageView<PixelT> tile(tile_bbox.width(), tile_bbox.height());
  rsrc->read(tile.buffer(), tile_bbox);
  ImageView<PixelT> reduced_tile = subsample(tile, pow(2,(num_levels-1) - level));
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
  return Vector2(0,1.0); // TODO: Implement this properly
}

PixelRGBA<float32> ImageTileGenerator::sample(int /*x*/, int /*y*/, int /*level*/, int /*transaction_id*/) {
  PixelRGBA<float32> result; // TODO: Implement this properly
  return result;
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
  return m_rsrc->block_size();
}

int32 ImageTileGenerator::num_levels() const {
  int32 max_dimension = std::max(this->cols(), this->rows());
  int32 max_tilesize = std::max(this->tile_size()[0], this->tile_size()[1]);
  return ceil(log(float(max_dimension) / max_tilesize) / log(2));
}
