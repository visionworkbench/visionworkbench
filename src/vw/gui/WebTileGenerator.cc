// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <vw/gui/WebTileGenerator.h>
#include <vw/FileIO/MemoryImageResource.h>
#include <vw/Image/PixelTypes.h>

using namespace vw::platefile;

namespace {
  size_t DEFAULT_BUFFER_SIZE = 32768;
}

namespace vw { namespace gui {

BlockingDownloader::BlockingDownloader() {
  m_http  = new QHttp(this);

  bool ret;
  ret = connect(m_http, SIGNAL(responseHeaderReceived(const QHttpResponseHeader&)), this, SLOT(onResponseHeaderReceived(const QHttpResponseHeader&)));
  VW_ASSERT(ret, NotFoundErr() << "Failed to connect to responseHeaderReceived");
  ret = connect(m_http, SIGNAL(readyRead(const QHttpResponseHeader&)), this, SLOT(onReadyRead(const QHttpResponseHeader&)));
  VW_ASSERT(ret, NotFoundErr() << "Failed to connect to readyRead");
  ret = connect(m_http, SIGNAL(requestFinished(int, bool)), this, SLOT(onRequestFinished(int, bool)));
  VW_ASSERT(ret, NotFoundErr() << "Failed to connect to onRequestFinished");
}

BlockingDownloader::~BlockingDownloader() {
  m_http->disconnect();
  delete m_http;
}

void BlockingDownloader::onResponseHeaderReceived(const QHttpResponseHeader& resp) {
  if (resp.hasContentLength())
    m_alloc = resp.contentLength();
  else
    m_alloc = DEFAULT_BUFFER_SIZE;
  m_result->data.reset(new uint8[m_alloc]);
}

void BlockingDownloader::onReadyRead(const QHttpResponseHeader& /*resp*/) {
  while (m_http->bytesAvailable() > ssize_t(m_alloc - m_result->size)) {
    vw_out(WarningMessage) << "Resizing read buffer from " << m_alloc << " to " << m_alloc*2 << "!\n";
    m_alloc *= 2;
    boost::shared_array<const uint8> newdata(new uint8[m_alloc]);
    std::copy(m_result->data.get(), m_result->data.get()+m_alloc/2, const_cast<uint8*>(newdata.get()));
    m_result->data = newdata;
  }
  int64 r = m_http->read(reinterpret_cast<char*>(const_cast<uint8*>(m_result->data.get())) + m_result->size, m_alloc - m_result->size);
  VW_ASSERT(r >= 0, IOErr() << "Failed to read bytes from wire");
  m_result->size += r;
}

BlockingDownloader::Result* BlockingDownloader::get(const Url& u_, int transaction, bool exact) {
  VW_ASSERT(!u_.hostname().empty(), ArgumentErr() << "Need a hostname");

  Url url(u_);

  url.query().set("nocache", "1");
  url.query().set("transaction_id", vw::stringify(transaction));
  if (exact)
    url.query().set("exact", "1");

  vw_out() << "  --> Fetching " << url << " ";

  m_http->setHost(url.hostname().c_str(), url.port() ? url.port() : 80);

  m_result.reset(new Result());
  m_result->size = m_alloc = 0;

  QEventLoop wait;
  connect(m_http, SIGNAL(done(bool)), &wait, SLOT(quit()));
  m_request = m_http->get(url.string().c_str(), 0);
  wait.exec();

  vw_out() << "[" << m_result->status << ", " << m_result->size << "]";

  if (m_result->status > 0) {
    vw_out() << '\n';
    return m_result.release();
  }
  else {
    if (m_result->status == 0)
      m_result->msg += "Timed out?";
    vw_out() << " [" << m_result->msg << "]\n";
    return 0;
  }
}

void BlockingDownloader::onRequestFinished(int request_id, bool error) {
  if (request_id != m_request)
    return;

  if (error) {
    m_result->msg = std::string("Connection Error: ") + m_http->errorString().toStdString();
    m_result->status = -1;
    return;
  }

  QHttpResponseHeader hdr = m_http->lastResponse();
  m_result->status   = hdr.statusCode();
  m_result->mimetype = std::string(hdr.contentType().toStdString());
  // result data is set already, behind the scenes
}

WebTileGenerator::WebTileGenerator(const Url& base_url, int levels)
  : m_tile_size(256), m_levels(levels), m_base_url(base_url) {}

boost::shared_ptr<SrcImageResource> WebTileGenerator::generate_tile(TileLocator const& tile_info) {

  Url url(m_base_url);
  Url::split_t path = url.path_split();
  path.reserve(path.size()+3);
  path.push_back(vw::stringify(tile_info.level));
  path.push_back(vw::stringify(tile_info.col));
  path.push_back(vw::stringify(tile_info.row) + ".png");
  url.path_join(path);

  BlockingDownloader d;
  boost::scoped_ptr<BlockingDownloader::Result> r(d.get(url, tile_info.transaction_id, tile_info.exact_transaction_id_match));

  typedef PixelRGBA<float> pixel_t;
  typedef boost::shared_ptr<SrcImageResource> Ptr;
  if (!r)
    return Ptr(make_point_src(pixel_t(1.f, 0.f, 0.f, 1.f)));
  else if (r->status == 404)
    return Ptr(make_point_src(pixel_t(0.f, 0.1f, 0.f, 1.f)));
  else if (r->status != 200)
    return Ptr(make_point_src(pixel_t(1.f, 0.f, 0.f, 1.f)));

  return Ptr(SrcMemoryImageResource::open( r->mimetype, r->data, r->size));
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
