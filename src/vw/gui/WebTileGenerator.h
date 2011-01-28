// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifndef __VW_GUI_WEBTILEGENERATOR_H__
#define __VW_GUI_WEBTILEGENERATOR_H__

#include <vw/gui/TileGenerator.h>
#include <vw/Plate/HTTPUtils.h>

#include <QThread>
#include <QHttp>
#include <QBuffer>

namespace vw {
namespace gui {

  class HttpDownloadThread : public QThread {
    Q_OBJECT

    struct RequestBuffer {
        platefile::Url url;
        bool ready;
        QBuffer buffer;
        vw::ImageView<vw::PixelRGBA<float> > result;
      public:
        RequestBuffer(const platefile::Url& url_)
          : url(url_), ready(false) { buffer.open(QIODevice::WriteOnly); }
    };
    typedef boost::shared_ptr<RequestBuffer> RequestBufferPtr;

    typedef std::map<int, RequestBufferPtr> map_t;

    const platefile::Url m_base_url;
    const boost::shared_ptr<QHttp> m_http;
    map_t m_requests;
    vw::Mutex m_mutex;

  protected:
    void run();
  public:
    HttpDownloadThread(const platefile::Url& u);
    virtual ~HttpDownloadThread();
    int get(const std::vector<std::string>& path, int transaction_id, bool exact_transaction_id_match);
    bool result_available(int request_id);
    vw::ImageView<vw::PixelRGBA<float> > pop_result(int request_id);
  public slots:
    void request_finished(int id, bool error);
  };

  class WebTileGenerator : public TileGenerator {
    int m_tile_size;
    int m_levels;
    HttpDownloadThread m_download_thread;

  public:
    WebTileGenerator(const platefile::Url& url, int levels);
    virtual ~WebTileGenerator() {}

    virtual boost::shared_ptr<ViewImageResource> generate_tile(TileLocator const& tile_info);
    virtual Vector2 minmax();
    virtual PixelRGBA<float> sample(int x, int y, int level, int transaction_id);

    virtual int cols() const;
    virtual int rows() const;
    virtual PixelFormatEnum pixel_format() const;
    virtual ChannelTypeEnum channel_type() const;
    virtual Vector2i tile_size() const;
    virtual int32 num_levels() const;
  };


}} // namespace vw::gui

#endif

