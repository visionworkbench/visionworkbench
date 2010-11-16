// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifndef __VW_GUI_TILEGENERATOR_H__
#define __VW_GUI_TILEGENERATOR_H__

#include <vw/Core/FundamentalTypes.h>
#include <vw/Math/BBox.h>
#include <vw/Image/ImageView.h>
#include <vw/FileIO/DiskImageResource.h>
#include <vw/Image/PixelTypes.h>
#include <vw/Image/ViewImageResource.h>
#include <vw/Plate/PlateFile.h>

// Qt
#include <QObject>
#include <QThread>
#include <QHttp>
#include <QBuffer>


namespace vw {
namespace gui {

  struct TileLocator {
    int col;
    int row;
    int level;
    int transaction_id;
    bool exact_transaction_id_match;

    bool is_valid() const {
      return (col >= 0 && row >= 0 && col < pow(2, level) && row < pow(2, level));
    }
  };

  // Given a tile index, return the bounding box of that tile coverage
  // in the bottom (i.e. highest resolution) level of the source image
  // pyramid.
  BBox2i tile_to_bbox(Vector2i tile_size, int col, int row, int level, int max_level);

  std::list<TileLocator> bbox_to_tiles(Vector2i tile_size, BBox2i bbox, int level,
                                       int max_level, int transaction_id,
                                       bool exact_transaction_id_match);

  // --------------------------------------------------------------------------
  //                              TILE GENERATOR
  // --------------------------------------------------------------------------

  class TileGenerator : public QObject {
  public:
    virtual ~TileGenerator() {}
    virtual boost::shared_ptr<ViewImageResource> generate_tile(TileLocator const& tile_info) = 0;
    virtual Vector2 minmax() = 0;
    virtual PixelRGBA<float> sample(int x, int y, int level, int transaction_id) = 0;

    virtual int cols() const = 0;
    virtual int rows() const = 0;
    virtual PixelFormatEnum pixel_format() const = 0;
    virtual ChannelTypeEnum channel_type() const = 0;
    virtual Vector2i tile_size() const = 0;
    virtual int32 num_levels() const = 0;

    // Use this method to generate the correct type of TileGenerator
    // for a given filename.
    static boost::shared_ptr<TileGenerator> create(std::string filename);
  };

  // --------------------------------------------------------------------------
  //                         TEST PATTERN TILE GENERATOR
  // --------------------------------------------------------------------------

  class TestPatternTileGenerator : public TileGenerator {
    int m_tile_size;

  public:
    TestPatternTileGenerator(int tile_size) : m_tile_size(tile_size) {}
    virtual ~TestPatternTileGenerator() {}

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

  // --------------------------------------------------------------------------
  //                          WEB TILE GENERATOR
  // --------------------------------------------------------------------------
  class HttpDownloadThread : public QThread {
    Q_OBJECT

    QHttp *m_http;
    int m_current_request;

    struct RequestBuffer {
      boost::shared_ptr<QByteArray> bytes;
      boost::shared_ptr<QBuffer> buffer;
      bool finished;
      std::string file_type;
      int status;
      std::string url;
      vw::ImageView<vw::PixelRGBA<float> > result;

      RequestBuffer() {
        bytes.reset(new QByteArray());
        buffer.reset(new QBuffer(bytes.get()));
        buffer->open(QIODevice::WriteOnly);
        finished = false;
        status = 0;
      }
    };
    std::map<int, RequestBuffer> m_requests;
    vw::Mutex m_mutex;

  protected:
    void run();
  public:
    HttpDownloadThread();
    virtual ~HttpDownloadThread();
    int get(std::string url, int transaction_id, bool exact_transaction_id_match);
    bool result_available(int request_id);
    vw::ImageView<vw::PixelRGBA<float> > pop_result(int request_id);
  public slots:
    void request_started(int id);
    void response_header_received( const QHttpResponseHeader & resp );
    void request_finished(int id, bool error);
    void state_changed(int state);
  };

  class WebTileGenerator : public TileGenerator {
    Q_OBJECT

    int m_tile_size;
    int m_levels;
    std::string m_url;
    bool m_download_complete;
    int m_request_id;
    HttpDownloadThread m_download_thread;

  public:
    WebTileGenerator(std::string url, int levels);
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

  // --------------------------------------------------------------------------
  //                         PLATE FILE TILE GENERATOR
  // --------------------------------------------------------------------------

  class PlatefileTileGenerator : public TileGenerator {
    boost::shared_ptr<vw::platefile::PlateFile> m_platefile;
    int m_num_levels;

  public:
    PlatefileTileGenerator(std::string platefile_name);
    virtual ~PlatefileTileGenerator() {}

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

  // --------------------------------------------------------------------------
  //                             IMAGE TILE GENERATOR
  // --------------------------------------------------------------------------

  class ImageTileGenerator : public TileGenerator {
    std::string m_filename;
    boost::shared_ptr<SrcImageResource> m_rsrc;

  public:
    ImageTileGenerator(std::string filename);
    virtual ~ImageTileGenerator() {}
    virtual PixelRGBA<float> sample(int x, int y, int level, int transaction_id);

    virtual boost::shared_ptr<ViewImageResource> generate_tile(TileLocator const& tile_info);
    virtual Vector2 minmax();

    virtual int cols() const;
    virtual int rows() const;
    virtual PixelFormatEnum pixel_format() const;
    virtual ChannelTypeEnum channel_type() const;
    virtual Vector2i tile_size() const;
    virtual int32 num_levels() const;
  };


}} // namespace vw::gui

#endif // __VW_GUI_TILEGENERATOR_H__

