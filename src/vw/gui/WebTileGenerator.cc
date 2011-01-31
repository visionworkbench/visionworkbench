// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <vw/gui/WebTileGenerator.h>
#include <vw/FileIO/MemoryImageResource.h>

using namespace vw::platefile;

namespace vw { namespace gui {

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

  std::string filetype;
  if (error) {
    vw_out(WarningMessage, "console") << "Connection Error: " << m_http->errorString().toStdString() << "\n";
  } else {
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
  path[0] = vw::stringify(tile_info.level);
  path[1] = vw::stringify(tile_info.col);
  path[2] = vw::stringify(tile_info.row) + ".png";

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

}}
