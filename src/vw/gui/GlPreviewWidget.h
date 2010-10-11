// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// \file vwv_GlPreviewWidget.h
///
/// The Vision Workbench image viewer.
///
#ifndef __VW_GUI__PREVIEW_GL_WIDGET_H__
#define __VW_GUI_PREVIEW_GL_WIDGET_H__

// Qt
#include <QGLWidget>
#include <QGLFormat>
#include <QPoint>

// Vision Workbench
#include <vw/Core/Thread.h>
#include <vw/Core/Log.h>
#include <vw/Image/ImageResource.h>
#include <vw/Image/ImageViewRef.h>
#include <vw/Image/ImageView.h>
#include <vw/Image/Manipulation.h>
#include <vw/Image/Statistics.h>
#include <vw/FileIO/DiskImageResource.h>
#include <vw/FileIO/DiskImageView.h>
#include <vw/Math/BBox.h>
#include <vw/Math/Vector.h>

// STL
#include <string>
#include <list>

#include <vw/gui/TextureCache.h>

class QMouseEvent;
class QWheelEvent;
class QPoint;

namespace vw {
namespace gui {

  // A simple class for keeping track of crosshair locations and colors.
  class PointList {
    std::list<vw::Vector2> m_points;
    vw::Vector3 m_color;
  public:
    PointList(vw::Vector3 const& color) : m_color(color) {}
    PointList(std::list<vw::Vector2> const& points, vw::Vector3 const& color) :
      m_color(color) {
      this->push_back(points);
    }

    std::list<vw::Vector2> const& points() const { return m_points; }
    vw::Vector3 color() const { return m_color; }

    void push_back(vw::Vector2 pt) { m_points.push_back(pt); }
    void push_back(std::list<vw::Vector2> pts) {
      std::list<vw::Vector2>::iterator iter  = pts.begin();
      while (iter != pts.end()) {
        m_points.push_back(*iter);
        ++iter;
      }
    }
  };


  class GlPreviewWidget : public QGLWidget, public CachedTextureRenderer {
    Q_OBJECT

    public:

    // Constructors/Destructor
    GlPreviewWidget(QWidget *parent, std::string filename, QGLFormat const& frmt, int transaction_id);
    virtual ~GlPreviewWidget();

    virtual GLuint allocate_texture(boost::shared_ptr<ViewImageResource> tile);
    virtual void deallocate_texture(GLuint texture_id);

    // Set a default size for this widget.  This is usually overridden
    // by parent views.
    virtual QSize sizeHint () const { return QSize(500,500); }

    // Image Manipulation Methods
    void zoom(float scale);
    void size_to_fit();

  public slots:

    // This timer is used by the Texture Fetch thread to inform the
    // GlPreviewWidget when new textures are available for drawing.
    // Timer callback is called 30 times per second.
    void timer_callback() {
      if (m_needs_redraw) {
        update();
        m_needs_redraw = false;
      }
    }

    void set_nodata_value(float nodata_value) {
      m_nodata_value = nodata_value;
      m_use_nodata = 1;
    }

    void normalize() {
      std::cout << "CALLING NORMALIZE!\n";
      Vector2 result = m_tile_generator->minmax();
      std::cout << "Got result = " << result << "\n";
      m_image_min = result[0];
      m_image_max = result[1];
      m_offset = -m_image_min;
      m_gain = 1/(m_image_max-m_image_min);
      update();
    }

  protected:

    // Setup
    void initializeGL();
    void resizeGL(int width, int height);

    // Event handlers
    void paintEvent(QPaintEvent *event);
    void mousePressEvent(QMouseEvent *event);
    void mouseMoveEvent(QMouseEvent *event);
    void mouseDoubleClickEvent(QMouseEvent *event);
    void wheelEvent(QWheelEvent *event);
    void enterEvent(QEvent *event);
    void leaveEvent(QEvent *event);
    void keyPressEvent(QKeyEvent *event);

  private:
    // Drawing is driven by QPaintEvents, which call out to drawImage()
    // and drawLegend()
    void drawImage();
    void drawLegend(QPainter *painter);
    void updateCurrentMousePosition();

    // Image & OpenGL
    GLuint m_glsl_program;
    bool m_show_legend;
    bool m_bilinear_filter;
    bool m_use_colormap;
    bool m_show_tile_boundaries;

    // Timers and updates
    QTimer *m_timer;

    // Image tiles and the texture cache
    boost::shared_ptr<TileGenerator> m_tile_generator;
    boost::shared_ptr<GlTextureCache> m_gl_texture_cache;
    PixelRGBA<float> m_last_pixel_sample;

    // Adjustment mode
    enum AdjustmentMode { TransformAdjustment, GainAdjustment,
                          OffsetAdjustment, GammaAdjustment };
    AdjustmentMode m_adjust_mode;

    // Mouse positions and legend information
    QPoint lastPos;
    QPoint currentImagePos;
    std::string m_legend_status;

    // Dimensions & stats
    int m_viewport_width;
    int m_viewport_height;
    vw::float32 m_image_min;
    vw::float32 m_image_max;
    vw::float32 m_nodata_value;
    int m_use_nodata;

    // Image Parameters
    vw::BBox2 m_current_viewport;
    float m_gain;
    float m_offset;
    float m_gamma;
    int m_current_transaction_id;
    bool m_exact_transaction_id_match;
    int m_current_level;

    enum DisplayChannel { DisplayRGBA = 0, DisplayR, DisplayG, DisplayB, DisplayA };
    int m_display_channel;
    int m_colorize_display;
    int m_hillshade_display;
  };

}} // namespace vw::gui

#endif  // __VW_GUI_PREVIEW_GL_WIDGET_H__
