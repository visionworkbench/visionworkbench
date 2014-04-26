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


/// \file vwv_GlPreviewWidget.cc
///
/// The Vision Workbench image viewer.
///

// Qt
#include <QtGui>

// Vision Workbench
#include <vw/Image.h>
#include <vw/FileIO.h>
#include <vw/gui/GlPreviewWidget.h>
using namespace vw;
using namespace vw::gui;

void imageData::read(std::string const& image){
  name = image;
  img = DiskImageView<float>(name);

  boost::shared_ptr<DiskImageResource> rsrc(DiskImageResource::open(name));
  ChannelTypeEnum channel_type = rsrc->channel_type();
  double nodata_val = -32768;
  if ( rsrc->has_nodata_read() ) {
    nodata_val = rsrc->nodata_read();
  }
  
  // Must scale the image values to uint8
  if(channel_type == VW_CHANNEL_UINT8){
    
    // Set no-data pixels to 0    
    for (int col = 0; col < img.cols(); col++){
      for (int row = 0; row < img.rows(); row++){
        img(col, row) = std::max((int)img(col, row), 0);
      }
    }
    
  }else{

    // Normalize to 0 - 255
    double mn = DBL_MAX, mx = -DBL_MAX;
    for (int col = 0; col < img.cols(); col++){
      for (int row = 0; row < img.rows(); row++){
        if (img(col, row) <= nodata_val) continue;
        if (img(col, row) < mn) mn = img(col, row);
        if (img(col, row) > mx) mx = img(col, row);
      }
    }
    if (mn >= mx){
      for (int col = 0; col < img.cols(); col++){
        for (int row = 0; row < img.rows(); row++){
          img(col, row) = 0.0;
        }
      }
    }else{
      for (int col = 0; col < img.cols(); col++){
        for (int row = 0; row < img.rows(); row++){
          img(col, row) = round(255*(std::max(double(img(col, row)), mn)
                                       - mn)/(mx-mn));
        }
      }
    }
  }
  
}

vw::Vector2 vw::gui::QPoint2Vec(QPoint const& qpt) {
  return vw::Vector2(qpt.x(), qpt.y());
}
QPoint vw::gui::Vec2QPoint(vw::Vector2 const& V) {
  return QPoint(round(V.x()), round(V.y()));
}

// --------------------------------------------------------------
//               GlPreviewWidget Public Methods
// --------------------------------------------------------------

GlPreviewWidget::GlPreviewWidget(QWidget *parent,
                                 std::vector<std::string> const& images,
                                 int transaction_id)
  : QWidget(parent)
{

  // Set default values
  m_nodata_value = 0;
  m_use_nodata = 0;
  m_image_min = 0;
  m_image_max = 1.0;

  // Set some reasonable defaults
  m_bilinear_filter = true;
  m_use_colormap = false;
  m_adjust_mode = NoAdjustment;
  m_display_channel = DisplayRGBA;
  m_colorize_display = false;
  m_hillshade_display = false;
  m_show_tile_boundaries = false;

  // Set up shader parameters
  m_gain = 1.0;
  m_offset = 0.0;
  m_gamma = 1.0;
  m_current_transaction_id = transaction_id;
  m_exact_transaction_id_match = false;

  // Set mouse tracking
  this->setMouseTracking(true);

  // Set the size policy that the widget can grow or shrink and still
  // be useful.
  this->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
  this->setFocusPolicy(Qt::ClickFocus);


  int num_images = images.size();
  m_images.resize(num_images);
  for (int i = 0; i < num_images; i++){
    m_images[i].read(images[i]);
    m_images_box.grow(bounding_box(m_images[i].img));
  }
  
  size_to_fit();
}


GlPreviewWidget::~GlPreviewWidget() {
}

void GlPreviewWidget::size_to_fit() {
  double aspect = double(m_window_width) / m_window_height;
  int maxdim = std::max(m_images_box.width(),m_images_box.height());
  if (m_images_box.width() > m_images_box.height()) {
    double width = maxdim;
    double height = maxdim/aspect;
    double extra = height - m_images_box.height();
    m_current_viewport = BBox2(Vector2(0.0, -extra/2),
                                Vector2(width, height-extra/2));
  } else {
    double width = maxdim*aspect;
    double height = maxdim;
    double extra = width - m_images_box.width();
    m_current_viewport = BBox2(Vector2(-extra/2, 0.0),
                                Vector2(width-extra/2, height));
  }
  update();
}

void GlPreviewWidget::zoom(double scale) {
  // Check to make sure we haven't hit our zoom limits...
  if (m_current_viewport.width()/scale > 1.0 &&
      m_current_viewport.height()/scale > 1.0 &&
      m_current_viewport.width()/scale < 20*m_images_box.width() &&
      m_current_viewport.height()/scale < 20*m_images_box.height()) {
    m_current_viewport = (m_current_viewport - m_curr_world_pos) / scale + m_curr_world_pos;
  }
  update(); // will call paintEvent()
}

void GlPreviewWidget::resizeEvent(QResizeEvent*){
  QRect v       = this->geometry();
  m_window_width = v.width();
  m_window_height = v.height();
  size_to_fit();
  return;
}


// --------------------------------------------------------------
//             GlPreviewWidget Private Methods
// --------------------------------------------------------------

void GlPreviewWidget::drawImage(QPainter* paint) {

  // The portion of the image to draw
  for (int i = 0; i < (int)m_images.size(); i++){
    
    BBox2i image_box = m_current_viewport;
    image_box.crop(bounding_box(m_images[i].img));
    
    // See where it fits on the screen
    BBox2i pixel_box;
    pixel_box.grow(world2pixel(image_box.min()));
    pixel_box.grow(world2pixel(image_box.max()));
    
    ImageView<float> img = crop(m_images[i].img, image_box);
    QImage qimg(img.cols(), img.rows(), QImage::Format_RGB888);
    for (int x = 0; x < img.cols(); ++x) {
      for (int y = 0; y < img.rows(); ++y) {
        qimg.setPixel(x, y, qRgb(img(x, y), img(x, y), img(x, y)));
      }
    }
    
    QRect rect(pixel_box.min().x(), pixel_box.min().y(),
               pixel_box.width(), pixel_box.height());
    paint->drawImage (rect, qimg);
  }
  
  return;
    
}

vw::Vector2 GlPreviewWidget::world2pixel(vw::Vector2 const& p){
  // Convert a position in the world coordinate system to a pixel value,
  // relative to the image seen on screen (the origin is an image corner).
  double x = m_window_width*((p.x() - m_current_viewport.min().x())
                             /m_current_viewport.width());
  double y = m_window_height*((p.y() - m_current_viewport.min().y())
                              /m_current_viewport.height());
  return vw::Vector2(round(x), round(y));
}

vw::Vector2 GlPreviewWidget::pixel2world(vw::Vector2 const& pix){
  // Convert a pixel on the screen (the origin is an image corner),
  // to global world coordinates.
  double x = m_current_viewport.min().x()
    + m_current_viewport.width() * double(pix.x()) / m_window_width;
  double y = m_current_viewport.min().y()
    + m_current_viewport.height() * double(pix.y()) / m_window_height;
  return vw::Vector2(x, y);
}

void GlPreviewWidget::updateCurrentMousePosition() {
  m_curr_world_pos = pixel2world(m_curr_pixel_pos);
}

// --------------------------------------------------------------
//             GlPreviewWidget Event Handlers
// --------------------------------------------------------------

void GlPreviewWidget::paintEvent(QPaintEvent * /* event */) {
  QPainter paint(this);
  drawImage(&paint);
}

void GlPreviewWidget::mousePressEvent(QMouseEvent *event) {
  m_mouse_press_pos = event->pos();

  m_curr_pixel_pos = QPoint2Vec(m_mouse_press_pos);
  m_last_gain = m_gain;     // Store this so the user can do linear
  m_last_offset = m_offset; // and nonlinear steps.
  m_last_gamma = m_gamma;
  m_last_viewport_min = QPoint( m_current_viewport.min().x(),
                                m_current_viewport.min().y() );
  updateCurrentMousePosition();
}

void GlPreviewWidget::mouseReleaseEvent ( QMouseEvent *event ){

  QPoint mouse_rel_pos = event->pos();

  int tol = 5;
  if (std::abs(mouse_rel_pos.x() - m_mouse_press_pos.x()) < tol &&
      std::abs(mouse_rel_pos.y() - m_mouse_press_pos.y()) < tol
      ){
    // If the mouse was released too close to where it was clicked,
    // do nothing.
    return;
  }

  m_current_viewport -= (pixel2world(QPoint2Vec(event->pos())) -
                         pixel2world(QPoint2Vec(m_mouse_press_pos)));

  update(); // will call paintEvent()
  
  return;
}

void GlPreviewWidget::mouseMoveEvent(QMouseEvent *event) {
  // Diff variables are just the movement of the mouse normalized to
  // 0.0-1.0;
  double x_diff = double(event->x() - m_curr_pixel_pos.x()) / m_window_width;
  double y_diff = double(event->y() - m_curr_pixel_pos.y()) / m_window_height;
  double width = m_current_viewport.width();
  double height = m_current_viewport.height();

  // Right mouse button just kicks up the gain of all mouse actions
  if (event->buttons() & Qt::LeftButton || event->buttons() & Qt::RightButton) {

    if ( event->buttons() & Qt::RightButton ) {
      x_diff *= 2;
      y_diff *= 2;
    }

    std::ostringstream s;
    switch (m_adjust_mode) {

    case NoAdjustment:
      break;
    case TransformAdjustment:
      m_current_viewport.min().x() =
        m_last_viewport_min.x() - x_diff * width;
      m_current_viewport.min().y() =
        m_last_viewport_min.y() - y_diff * height;
      m_current_viewport.max().x() =
        m_last_viewport_min.x() - x_diff * width + width;
      m_current_viewport.max().y() =
        m_last_viewport_min.y() - y_diff * height + height;
      break;

    case GainAdjustment:
      m_gain = m_last_gain * pow(2.0,x_diff);
      s << "Gain: " << (log(m_gain)/log(2)) << " f-stops\n";
      break;

    case OffsetAdjustment:
      m_offset = m_last_offset +
        (pow(100,fabs(x_diff))-1.0)*(x_diff > 0 ? 0.1 : -0.1);
      s << "Offset: " << m_offset << "\n";
      break;

    case GammaAdjustment:
      m_gamma = m_last_gamma * pow(2.0,x_diff);
      s << "Gamma: " << m_gamma << "\n";
      break;
    }
  }

  updateCurrentMousePosition();
}

void GlPreviewWidget::mouseDoubleClickEvent(QMouseEvent *event) {
  m_curr_pixel_pos = QPoint2Vec(event->pos());
  updateCurrentMousePosition();
}

void GlPreviewWidget::wheelEvent(QWheelEvent *event) {
  int num_degrees = event->delta();
  double num_ticks = double(num_degrees) / 360;

  // 2.0 chosen arbitrarily here as a reasonable scale factor giving good
  // sensitivity of the mousewheel. Shift zooms 50 times slower.
  double scale_factor = 2;
  if (event->modifiers() & Qt::ShiftModifier)
    scale_factor *= 50;

  double mag = fabs(num_ticks/scale_factor);
  double scale = 1;
  if (num_ticks > 0)
    scale = 1+mag;
  else if (num_ticks < 0)
    scale = 1-mag;

  zoom(scale);

  m_curr_pixel_pos = QPoint2Vec(event->pos());
  updateCurrentMousePosition();
}


void GlPreviewWidget::enterEvent(QEvent */*event*/) {
}

void GlPreviewWidget::leaveEvent(QEvent */*event*/) {
}

void GlPreviewWidget::keyPressEvent(QKeyEvent *event) {

  std::ostringstream s;

  switch (event->key()) {
  case Qt::Key_Plus:   // Increase transaction id
  case Qt::Key_Minus:  // Decrease transaction id
  case Qt::Key_F:  // Size to fit
    size_to_fit();
    break;
  case Qt::Key_I:  // Toggle bilinear/nearest neighbor interp
    m_bilinear_filter = !m_bilinear_filter;
    break;
  case Qt::Key_C:  // Activate colormap
    m_use_colormap = !m_use_colormap;
    break;
  case Qt::Key_H:  // Activate hillshade
    if ( m_hillshade_display == 0 ) {
      m_hillshade_display = 1;
    } else {
      m_hillshade_display *= 3;
    }
    if ( m_hillshade_display > 100 || m_hillshade_display < 0 )
      m_hillshade_display = 0;
    break;
  case Qt::Key_T:  // Activate tile boundaries
    m_show_tile_boundaries = !m_show_tile_boundaries;
    break;
  case Qt::Key_G:  // Gain adjustment mode
    if (m_adjust_mode == GainAdjustment) {
      m_adjust_mode = TransformAdjustment;
    } else {
      m_adjust_mode = GainAdjustment;
      s << "Gain: " << log(m_gain)/log(2) << " f-stops\n";
    }
    break;
  case Qt::Key_O:  // Offset adjustment mode
    if (m_adjust_mode == OffsetAdjustment) {
      m_adjust_mode = TransformAdjustment;
    } else {
      m_adjust_mode = OffsetAdjustment;
      s << "Offset: " << m_offset;
    }
    break;
  case Qt::Key_V:  // Gamma adjustment mode
    if (m_adjust_mode == GammaAdjustment) {
      m_adjust_mode = TransformAdjustment;
    } else {
      m_adjust_mode = GammaAdjustment;
      s << "Gamma: " << m_gamma;
    }
    break;
  case Qt::Key_1:  // Display red channel only
    m_display_channel = DisplayR;
    break;
  case Qt::Key_2:  // Display green channel only
    m_display_channel = DisplayG;
    break;
  case Qt::Key_3:  // Display blue channel only
    m_display_channel = DisplayB;
    break;
  case Qt::Key_4:  // Display alpha channel only
    m_display_channel = DisplayA;
    break;
  case Qt::Key_0:  // Display all color channels
    m_display_channel = DisplayRGBA;
    break;
  default:
    QWidget::keyPressEvent(event);
  }
}

