// __BEGIN_LICENSE__
// 
// Copyright (C) 2008 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration
// (NASA).  All Rights Reserved.
// 
// Copyright 2008 Carnegie Mellon University. All rights reserved.
// 
// This software is distributed under the NASA Open Source Agreement
// (NOSA), version 1.3.  The NOSA has been approved by the Open Source
// Initiative.  See the file COPYING at the top of the distribution
// directory tree for the complete NOSA document.
// 
// THE SUBJECT SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY WARRANTY OF ANY
// KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT
// LIMITED TO, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL CONFORM TO
// SPECIFICATIONS, ANY IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR
// A PARTICULAR PURPOSE, OR FREEDOM FROM INFRINGEMENT, ANY WARRANTY THAT
// THE SUBJECT SOFTWARE WILL BE ERROR FREE, OR ANY WARRANTY THAT
// DOCUMENTATION, IF PROVIDED, WILL CONFORM TO THE SUBJECT SOFTWARE.
//
// __END_LICENSE__

/// \file vwv_GlPreviewWidget.h
///
/// The Vision Workbench image viewer. 
///

#ifndef __VWV_PREVIEW_GL_WIDGET_H__
#define __VWV_PREVIEW_GL_WIDGET_H__

// Qt
#include <QGLWidget>
#include <QGLFormat>
#include <QPoint>

// Vision Workbench
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

class QMouseEvent;
class QWheelEvent;
class QPoint;

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


class GlPreviewWidget : public QGLWidget {
  Q_OBJECT

public:

  // Constructors/Destructor
  GlPreviewWidget(QWidget *parent) : QGLWidget(parent) {
    setup();
  }

  GlPreviewWidget(QWidget *parent, std::string filename) : QGLWidget(parent) {
    setup();
    set_image_from_file(filename);
  }

  template <class ViewT>
  GlPreviewWidget(QWidget *parent, vw::ImageViewBase<ViewT> const& view) : QGLWidget(parent) {
    setup();
    set_image(view);
  }

  virtual ~GlPreviewWidget();


  // Set a default size for this widget.  This is usually overridden
  // by parent views.
  virtual QSize sizeHint () const { return QSize(500,500); }

  // Image Manipulation Methods
  void zoom(float scale);
  void normalize_image();
  void add_crosshairs(std::list<vw::Vector2> const& points, vw::Vector3 const& color);
  void clear_crosshairs(); 
  void size_to_fit();

  /// Replace the current image in the widget with the supplied image view.
  template <class ViewT>
  void set_image(vw::ImageViewBase<ViewT> const& view) {
    m_true_image_format = view.format();
    vw::ImageView<vw::PixelRGBA<vw::float32> > im = vw::channel_cast<vw::float32>(view.impl());
    m_image = im;
    rebind_textures();
    update();
    size_to_fit();
  }

public slots:

  void set_nodata_value(float nodata_value) { 
    m_nodata_value = nodata_value; 
    m_use_nodata = 1; 
  }

  void set_data_range(float lo, float hi) { 
    m_image_min = lo; 
    m_image_max = hi; 
    if (m_image_max > 1.0 || m_image_min < 0.0)
      normalize_image();
  }

  /// Replace the current image in the widget with the image in the
  /// supplied file.
  void set_image_from_file(std::string const& filename) {
    
    // First, we will determine the native format of the data in the
    // file.
    vw::DiskImageResource *rsrc = vw::DiskImageResource::open(filename);
    m_true_image_format = rsrc->format();
    delete rsrc;

    // Then we open it, and set it as the current image for rendering.
    vw::DiskImageView<vw::PixelRGBA<vw::float32> > input_image(filename);
    m_image = input_image;
    rebind_textures();
    update();
    size_to_fit();
  }

protected:

  // Setup
  void setup();
  void initializeGL();
  void resizeGL(int width, int height);
  void rebind_textures();

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
  void deallocate_textures();
  
  // Image & OpenGL
  vw::ImageViewRef<vw::PixelRGBA<vw::float32> > m_image;
  std::vector<vw::BBox2i> m_bboxes;
  vw::ImageFormat m_true_image_format;
  std::vector<GLuint> m_textures;
  GLuint m_glsl_program;
  bool m_draw_texture;
  bool m_show_legend;
  bool m_bilinear_filter;
  bool m_use_colormap;
  
  // Adjustment mode
  enum AdjustmentMode { TransformAdjustment, GainAdjustment, OffsetAdjustment, GammaAdjustment };
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

  // Crosshair overlays
  std::vector<PointList> m_crosshairs;

  enum DisplayChannel { DisplayRGBA = 0, DisplayR, DisplayG, DisplayB, DisplayA };
  int m_display_channel;
};

#endif  // __VWV_PREVIEW_GL_WIDGET_H__
