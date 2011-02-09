// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifndef __VW_GUI_WEBTILEGENERATOR_H__
#define __VW_GUI_WEBTILEGENERATOR_H__

#include <vw/gui/TileGenerator.h>
#include <vw/Plate/HTTPUtils.h>

#include <QEventLoop>
#include <QHttp>
#include <QBuffer>

namespace vw {
namespace gui {

  class BlockingDownloader : public QObject {
      Q_OBJECT
    public:
      struct Result {
        int status;
        std::string msg;
        std::string mimetype;
        boost::shared_array<const uint8> data;
        size_t size;
      };

      BlockingDownloader();
      ~BlockingDownloader();
      Result* get(const platefile::Url& u_, int transaction, bool exact);

    private:
      QHttp* m_http;
      std::auto_ptr<Result> m_result;
      int m_request;
      bool m_done;
      size_t m_alloc;
    private slots:
      void onResponseHeaderReceived(const QHttpResponseHeader& resp);
      void onReadyRead(const QHttpResponseHeader& resp);
      void onRequestFinished(int request_id, bool error);
  };

  class WebTileGenerator : public TileGenerator {
    int m_tile_size;
    int m_levels;
    const platefile::Url m_base_url;

  public:
    WebTileGenerator(const platefile::Url& url, int levels);
    virtual ~WebTileGenerator() {}

    virtual boost::shared_ptr<SrcImageResource> generate_tile(TileLocator const& tile_info);
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

