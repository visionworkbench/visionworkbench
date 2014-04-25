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


/// \file vwv_GlPreviewWidget.h
///
/// The Vision Workbench image viewer.
///
#ifndef __VW_GUI__PREVIEW_GL_WIDGET_H__
#define __VW_GUI_PREVIEW_GL_WIDGET_H__

// Qt
#include <QWidget>
//#include <QGLFormat>
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
#include <vector>
#include <list>

#include <vw/gui/TextureCache.h>

class QMouseEvent;
class QWheelEvent;
class QPoint;
class QResizeEvent;

namespace vw {
namespace gui {

  // A class to keep all data associated with an image file
  struct imageData{
    std::string name;
    ImageView<float> img;
    void read(std::string const& image);
  };
  
  vw::Vector2 QPoint2Vec(QPoint const& qpt);
  QPoint Vec2QPoint(vw::Vector2 const& V);

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


  class GlPreviewWidget : public QWidget, public CachedTextureRenderer {
    Q_OBJECT

    public:

    // Constructors/Destructor
    GlPreviewWidget(QWidget *parent, std::vector<std::string> const& images,
                    //QGLFormat const& frmt,
                    int transaction_id);
    virtual ~GlPreviewWidget();

    virtual GLuint allocate_texture(boost::shared_ptr<SrcImageResource> tile);
    virtual void deallocate_texture(GLuint texture_id);

    // Set a default size for this widget.  This is usually overridden
    // by parent views.
    virtual QSize sizeHint () const { return QSize(500,500); }

    // Image Manipulation Methods
    void zoom(double scale);
    void size_to_fit();

  public slots:

    // This timer is used by the Texture Fetch thread to inform the
    // GlPreviewWidget when new textures are available for drawing.
    // Timer callback is called 30 times per second.
    void timer_callback() {
      update();
    }

    void set_nodata_value(double nodata_value) {
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
    }

  protected:

    // Setup
    void resizeEvent(QResizeEvent*);

    // Event handlers
    void paintEvent(QPaintEvent *event);
    void mousePressEvent(QMouseEvent *event);
    void mouseReleaseEvent(QMouseEvent *event);
    void mouseMoveEvent(QMouseEvent *event);
    void mouseDoubleClickEvent(QMouseEvent *event);
    void wheelEvent(QWheelEvent *event);
    void enterEvent(QEvent *event);
    void leaveEvent(QEvent *event);
    void keyPressEvent(QKeyEvent *event);

  private:
    // Drawing is driven by QPaintEvents, which call out to drawImage()
    // and drawLegend()
    void drawImage(QPainter* paint);
    void drawLegend(QPainter *paint);
    vw::Vector2 world2pixel(vw::Vector2 const& p);
    vw::Vector2 pixel2world(vw::Vector2 const& pix);
    void updateCurrentMousePosition();

    // Image & OpenGL
    GLuint m_glsl_program;
    bool m_show_legend;
    bool m_bilinear_filter;
    bool m_use_colormap;
    bool m_show_tile_boundaries;

    // Image tiles and the texture cache
    std::vector<imageData> m_images;
    BBox2i m_images_box;
    
    boost::shared_ptr<TileGenerator> m_tile_generator;
    boost::shared_ptr<GlTextureCache> m_gl_texture_cache;
    PixelRGBA<float> m_last_pixel_sample;

    // Adjustment mode
    enum AdjustmentMode { NoAdjustment,
                          TransformAdjustment, GainAdjustment,
                          OffsetAdjustment, GammaAdjustment };
    AdjustmentMode m_adjust_mode;

    // Mouse positions and legend information
    vw::Vector2 m_curr_pixel_pos, m_curr_world_pos;
    QPoint m_last_viewport_min, m_mouse_press_pos;
    std::string m_legend_status;

    // Dimensions and stats
    int m_window_width;  // the width  of the plotting window in screen pixels
    int m_window_height; // the height of the plotting window in screen pixels
    vw::float32 m_image_min;
    vw::float32 m_image_max;
    vw::float32 m_nodata_value;
    int m_use_nodata;

    // Image Parameters
    vw::BBox2 m_current_viewport;
    double m_gain, m_last_gain;
    double m_offset, m_last_offset;
    double m_gamma, m_last_gamma;
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
