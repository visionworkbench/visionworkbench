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

/// \file vwv_GlPreviewWidget.cc
///
/// The Vision Workbench image viewer. 
///

// Qt
#include <QtGui>

// Vision Workbench
#include <vw/Image.h>
#include <vw/FileIO.h>
using namespace vw;

#include <vw/gui/vwv_GlPreviewWidget.h>

const std::string g_FRAGMENT_PROGRAM =  
"uniform sampler2D tex;                             \n"
"                                                   \n"
"uniform float gain;\n"
"uniform float offset;\n"
"uniform float gamma;\n"
"uniform float nodata;\n"
"uniform int use_nodata;\n"
"uniform int display_channel;  // 0 = RGBA, 1 = R, 2 = G, 3 = B, 4 = A\n"
"\n"
"void main() { \n"
"  vec4 g = vec4(gain, gain, gain, 1.0);\n"
"  vec4 o = vec4(offset, offset, offset, 0.0);\n"
"  vec4 v = vec4(gamma, gamma, gamma, 1.0);\n"
"\n"
"  vec4 src_tex = texture2D(tex,gl_TexCoord[0].st);\n"
"  vec4 selected_tex = src_tex;\n"
"  if (display_channel == 1) {\n"
"    selected_tex = vec4(src_tex[0],src_tex[0],src_tex[0],src_tex[3]);\n"
"  } else if (display_channel == 2) {\n"
"    selected_tex = vec4(src_tex[1],src_tex[1],src_tex[1],src_tex[3]);\n"
"  } else if (display_channel == 3) {\n"
"    selected_tex = vec4(src_tex[2],src_tex[2],src_tex[2],src_tex[3]);\n"
"  } else if (display_channel == 4) {\n"
"    selected_tex = vec4(src_tex[3],src_tex[3],src_tex[3],1.0);\n"
"  } \n"
"  vec4 final_tex = pow(g*(selected_tex + o), v);\n"
"  if (use_nodata == 1) { \n"
"    if (src_tex[0] == nodata) { \n"
"      final_tex = vec4(0.0,0,0,0.0);\n"
"    } \n"
"  }\n"
"  gl_FragColor = final_tex; \n"
"}\n";


// --------------------------------------------------------------
//                       GLSL DEBUGGING
// --------------------------------------------------------------

void print_shader_info_log(GLuint obj)
{
  int infologLength = 0;
  int charsWritten  = 0;
  char *infoLog;

  glGetShaderiv(obj, GL_INFO_LOG_LENGTH,&infologLength);

  if (infologLength > 0) {
    infoLog = (char *)malloc(infologLength);
    glGetShaderInfoLog(obj, infologLength, &charsWritten, infoLog);
    std::ostringstream err;
    err << "<h4>An error occured while compiling the GLSL shader:</h4><p><h5><tt>" << infoLog << "</tt></h5>";
    QMessageBox::critical(0, "GLSL Shader Error", 
                            err.str().c_str());
    free(infoLog);
  }
}

void print_program_info_log(GLuint obj)
{
  int infologLength = 0;
  int charsWritten  = 0;
  char *infoLog;
  
  glGetProgramiv(obj, GL_INFO_LOG_LENGTH,&infologLength);
  
  if (infologLength > 0) {
    infoLog = (char *)malloc(infologLength);
    glGetProgramInfoLog(obj, infologLength, &charsWritten, infoLog);
    std::ostringstream err;
    err << "<h4>An error occured while linking the GLSL program:</h4><p><h5><tt>" << infoLog << "</tt></h5>";
    QMessageBox::critical(0, "GLSL Program Error", 
                          err.str().c_str());
    printf("%s\n",infoLog);
    free(infoLog);
  }
}

void check_gl_errors( void ) 
{ 
  GLenum errCode; 
  const GLubyte *errStr; 
  if ( ( errCode = glGetError( ) ) != GL_NO_ERROR ) { 
    errStr = gluErrorString( errCode ); 
    std::cout << "OpenGL ERROR (" << int(errCode) << "): " << errStr << "\n";
  } 
} 
// --------------------------------------------------------------
//               GlPreviewWidget Public Methods
// --------------------------------------------------------------

GlPreviewWidget::~GlPreviewWidget() {
  deallocate_textures();
}

void GlPreviewWidget::size_to_fit() {
  float aspect = float(m_viewport_width) / m_viewport_height;
  int maxdim = std::max(m_image.cols(),m_image.rows());
  if (m_image.cols() > m_image.rows()) {
    float width = maxdim;
    float height = maxdim/aspect;
    float extra = height - m_image.rows();
    m_current_viewport = BBox2(Vector2(0.0, -extra/2), 
                                Vector2(width, height-extra/2));
  } else {
    float width = maxdim*aspect;
    float height = maxdim;
    float extra = width - m_image.cols();
    m_current_viewport = BBox2(Vector2(-extra/2, 0.0), 
                                Vector2(width-extra/2, height));
  }
  update();
}

void GlPreviewWidget::zoom(float scale) {
  float mid_x = m_current_viewport.min().x() + m_current_viewport.width()/2;
  float mid_y = m_current_viewport.min().y() + m_current_viewport.height()/2;
  
  // Check to make sure we haven't hit our zoom limits...
  if (m_current_viewport.width()/scale > 1.0 && 
      m_current_viewport.height()/scale > 1.0 &&
      m_current_viewport.width()/scale < 4*m_image.cols() && 
      m_current_viewport.height()/scale < 4*m_image.rows()) {
    m_current_viewport.min().x() = (m_current_viewport.min().x() - mid_x) / scale + mid_x;
    m_current_viewport.max().x() = (m_current_viewport.max().x() - mid_x) / scale + mid_x;
    m_current_viewport.min().y() = (m_current_viewport.min().y() - mid_y) / scale + mid_y;
    m_current_viewport.max().y() = (m_current_viewport.max().y() - mid_y) / scale + mid_y;
    update();
  }
  m_show_legend = false;
}

void GlPreviewWidget::normalize_image() {
  m_offset = -m_image_min;
  m_gain = 1/(m_image_max-m_image_min);
  update();
}

// --------------------------------------------------------------
//             GlPreviewWidget Private Methods
// --------------------------------------------------------------

void GlPreviewWidget::drawImage() {

  // Make this context current, and store the current OpenGL state
  // before we start to modify it.
  makeCurrent();
  glPushAttrib(GL_ALL_ATTRIB_BITS);
  
  // Activate our GLSL fragment program and set up the uniform
  // variables in the shader
  glUseProgram(m_glsl_program);
  GLint gain_loc = glGetUniformLocation(m_glsl_program,"gain");
  glUniform1f(gain_loc,m_gain);
  GLint offset_loc = glGetUniformLocation(m_glsl_program,"offset");
  glUniform1f(offset_loc,m_offset);
  GLint gamma_loc = glGetUniformLocation(m_glsl_program,"gamma");
  glUniform1f(gamma_loc,m_gamma);
  GLint nodata_loc = glGetUniformLocation(m_glsl_program,"nodata");
  glUniform1f(nodata_loc,m_nodata_value);
  GLint use_nodata_loc = glGetUniformLocation(m_glsl_program,"use_nodata");
  glUniform1i(use_nodata_loc,m_use_nodata);
  GLint display_channel_loc = glGetUniformLocation(m_glsl_program,"display_channel");
  glUniform1i(display_channel_loc,m_display_channel);

  // Set the background color and viewport.
  qglClearColor(QColor(0, 25, 50)); // Bluish-green background
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glViewport(0,0,m_viewport_width,m_viewport_height);

  // Set up the orthographic view of the scene.  The exact extent of
  // the view onto the scene depends on the current panning and zoom
  // in the UI.
  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();
  glOrtho(m_current_viewport.min().x(), m_current_viewport.max().x(), 
          -m_current_viewport.max().y(), -m_current_viewport.min().y(),
          -1.0, 1.0);

  // Set up the modelview matrix, and bind the image as the texture we
  // are about to use.
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();

  glEnable( GL_TEXTURE_2D );
  for (unsigned int i=0; i < m_bboxes.size(); ++i) {

    // Bind the current texture.
    glBindTexture( GL_TEXTURE_2D, m_textures[i] );
    
    // Set up bilinear or nearest neighbor filtering.
    if (m_bilinear_filter) {
      // When the texture area is small, bilinear filter the closest mipmap
      glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST_MIPMAP_LINEAR );
      glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR );
    } else {
      // When the texture area is small, pick the nearest neighbor in the closest mipmap
      glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST_MIPMAP_NEAREST );
      glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST );
    }
    
    // We will draw the image as a texture on this quad.
    qglColor(Qt::white);
    glBegin(GL_QUADS);
    glTexCoord2d( 0.0 , 0.0); 
    glVertex2d( m_bboxes[i].min().x() , -(m_bboxes[i].min().y()) );
    glTexCoord2d( 0.0 , 1 ); 
    glVertex2d( m_bboxes[i].min().x() , -(m_bboxes[i].max().y()) );
    glTexCoord2d( 1.0 , 1.0 );
    glVertex2d( m_bboxes[i].max().x() , -(m_bboxes[i].max().y()) );
    glTexCoord2d( 1.0 , 0.0 ); 
    glVertex2d( m_bboxes[i].max().x() , -(m_bboxes[i].min().y()) );
    glEnd();
  }

  // Disable texture mapping and GLSL shaders
  glDisable( GL_TEXTURE_2D );
  glUseProgram(0);

  // Draw crosshairs
  glLineWidth(1.0);
  for (unsigned i = 0; i < m_crosshairs.size(); ++i) {
    Vector3 color = m_crosshairs[i].color();
    glColor3f(color[0], color[1], color[2]);
    glBegin(GL_LINES);
    std::list<Vector2>::const_iterator iter = m_crosshairs[i].points().begin();
    while (iter != m_crosshairs[i].points().end() ) {
      Vector2 point = *iter;
      glVertex2d( point[0]-3 , -point[1]);
      glVertex2d( point[0]+3 , -point[1]);
      glVertex2d( point[0], -point[1]-3);
      glVertex2d( point[0], -point[1]+3);
      ++iter;
    }
    glEnd();
  }    
  
  // Restore the previous OpenGL state so that we don't trample on the
  // QPainter elements of the window.
  glMatrixMode(GL_MODELVIEW);
  glPopMatrix();
  glMatrixMode(GL_PROJECTION);
  glPopMatrix();
  glPopAttrib();
}

void GlPreviewWidget::drawLegend(QPainter* painter) {

  // Extract the value for the pixel currently under the mouse
  PixelRGBA<float32> pix_value;
  if (currentImagePos.x() >= 0 && currentImagePos.x() < m_image.cols() &&
      currentImagePos.y() >= 0 && currentImagePos.y() < m_image.rows()) {
    pix_value = m_image(currentImagePos.x(), currentImagePos.y());
    
    const char* pixel_name = vw::pixel_format_name(m_true_image_format.pixel_format);
    const char* channel_name = vw::channel_type_name(m_true_image_format.channel_type);
    int num_channels = vw::num_channels(m_true_image_format.pixel_format);
    std::ostringstream pix_value_ostr;
    pix_value_ostr << "Pos: ( " << currentImagePos.x() << " " << currentImagePos.y() << " ) --> Val: [ ";

    // The following code is very messy and should be replaced and/or
    // simplified once we have better control over whether automatic
    // conversions happen when reading in images using the FileIO
    // code.
    bool round;
    float scale_factor;
    switch (m_true_image_format.channel_type) {
    case VW_CHANNEL_INT8:
      scale_factor = ScalarTypeLimits<int8>::highest();
      round = true;
      break;
    case VW_CHANNEL_BOOL:
    case VW_CHANNEL_CHAR:
    case VW_CHANNEL_UINT8:
    case VW_CHANNEL_GENERIC_1_BYTE:
      scale_factor = ScalarTypeLimits<uint8>::highest();
      round = true;
      break;
    case VW_CHANNEL_INT16:
      scale_factor = ScalarTypeLimits<int16>::highest();
      round = true;
      break;
    case VW_CHANNEL_UINT16:
    case VW_CHANNEL_GENERIC_2_BYTE:
      scale_factor = ScalarTypeLimits<uint16>::highest();
      round = true;
      break;
    case VW_CHANNEL_FLOAT16:
      scale_factor = 1.0;
      round = false;
      break;
    case VW_CHANNEL_INT32:
      scale_factor = ScalarTypeLimits<int32>::highest();
      round = true;
      break;
    case VW_CHANNEL_UINT32:
    case VW_CHANNEL_GENERIC_4_BYTE:
      scale_factor = ScalarTypeLimits<uint32>::highest();
      round = true;
      break;
    case VW_CHANNEL_FLOAT32:
      scale_factor = 1.0;
      round = false;
      break;
    case VW_CHANNEL_INT64:
      scale_factor = ScalarTypeLimits<int64>::highest();
      round = true;
      break;
    case VW_CHANNEL_UINT64:
      scale_factor = ScalarTypeLimits<uint64>::highest();
      round = true;
      break;
    case VW_CHANNEL_FLOAT64:
      scale_factor = 1.0;
      round = false;
      break;
    default:
      vw_throw( ArgumentErr() << "vwv_GlPreviewWidget() : Unrecognized or unsupported channel type (" << m_true_image_format.channel_type << ")." );
      return; // never reached
    }
    
    // We don't print out the channel value for the alpha channel for
    // those pixel types that have alpha.  This is messy too, and
    // should also be replaced by a better method...
    if (m_true_image_format.pixel_format == vw::VW_PIXEL_GRAYA ||
        m_true_image_format.pixel_format == vw::VW_PIXEL_RGBA ||
        m_true_image_format.pixel_format == vw::VW_PIXEL_GRAY_MASKED ||
        m_true_image_format.pixel_format == vw::VW_PIXEL_RGB_MASKED ||
        m_true_image_format.pixel_format == vw::VW_PIXEL_XYZ_MASKED ||
        m_true_image_format.pixel_format == vw::VW_PIXEL_LUV_MASKED ||
        m_true_image_format.pixel_format == vw::VW_PIXEL_LAB_MASKED) {
      num_channels -= 1;
    }

    for (unsigned i=0; i < num_channels; ++i) {
      if (round)
        pix_value_ostr << llroundf(pix_value[i] * scale_factor) << " ";
      else
        pix_value_ostr << (pix_value[i] * scale_factor) << " ";
    }

    pix_value_ostr << "]";


    const int Margin = 11;
    const int Padding = 6;
    QTextDocument textDocument;
    textDocument.setDefaultStyleSheet("* { color: #00FF00; font-family: courier, serif; font-size: 12 }");
    std::ostringstream legend_text;
    legend_text << "<p align=\"right\">" << m_legend_status << "<br>"
                << "[ " << m_image.cols() << " x " << m_image.rows() << " ] <br>"
                << "Pixel Format: " << pixel_name << "  Channel Type: " << channel_name << "<br>"
                << "Current Pixel Range: [ " << -m_offset << " " 
                << ( -m_offset+(1/m_gain) ) << " ]<br>"
                << pix_value_ostr.str() << "<br>"
                << "</p>";
    textDocument.setHtml(legend_text.str().c_str());
    textDocument.setTextWidth(textDocument.size().width());
  
    QRect rect(QPoint(0,0), textDocument.size().toSize() + QSize(2 * Padding, 2 * Padding));
    painter->translate(width() - rect.width() - Margin, height() - rect.height() - Margin);
    //     painter->setPen(QColor(255, 239, 239));
    //     painter->drawRect(rect);
    painter->translate(Padding, Padding);
    textDocument.drawContents(painter);
  }
}

void GlPreviewWidget::updateCurrentMousePosition() {
  float x_loc = m_current_viewport.min().x() + m_current_viewport.width() * float(lastPos.x()) / m_viewport_width;
  float y_loc = m_current_viewport.min().y() + m_current_viewport.height() * float(lastPos.y()) / m_viewport_height;
  currentImagePos = QPoint(x_loc,y_loc);
}


// --------------------------------------------------------------
//             GlPreviewWidget Setup Methods
// --------------------------------------------------------------
void GlPreviewWidget::setup() {
  if (!QGLFormat::hasOpenGL()) {
    vw::vw_out(0) << "This system has no OpenGL support.\nExiting\n\n";
    exit(1);
  }

  m_nodata_value = 0;
  m_use_nodata = 0;
  m_image_min = 0;
  m_image_max = 1.0;

  // Set some reasonable defaults
  m_draw_texture = true;
  m_show_legend = false;
  m_bilinear_filter = true;
  m_use_colormap = false;
  m_adjust_mode = TransformAdjustment;
  m_display_channel = DisplayRGBA;
  
  // Set up shader parameters
  m_gain = 1.0;
  m_offset = 0.0;
  m_gamma = 1.0;
  
  // Set mouse tracking
  this->setMouseTracking(true);

  // Set the size policy that the widget can grow or shrink and still
  // be useful.
  this->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
}

void GlPreviewWidget::deallocate_textures() {
  makeCurrent();
  for (unsigned i = 0; i < m_textures.size(); ++i) {
    glDeleteTextures(1,&(m_textures[i]));
  }
  m_textures.clear();
}

void GlPreviewWidget::rebind_textures() {
  makeCurrent();

  // Deallocate any previously allocated textures
  deallocate_textures();

  // The GPU has limited VRAM, so any one texture must be limited in
  // size.  We query the GPU for its max texture size here, and choose
  // a final texture size that is 1/4 this size.  (This leaves ample
  // room for use to store GL_RGBA16F_ARB (16-bit float) images on the
  // GPU.)
  GLint texSize; 
  glGetIntegerv(GL_MAX_TEXTURE_SIZE, &texSize);
  m_bboxes = vw::image_blocks(m_image, 512, 512);
  
  // Allocate the OpenGL texture patches
  m_textures.resize(m_bboxes.size());
  std::cout << "\t--> Caching textures... " << std::flush;
  for (unsigned i=0; i<m_bboxes.size(); ++i) {
    glGenTextures(1, &(m_textures[i]) );
    glBindTexture(GL_TEXTURE_2D, m_textures[i]);
    
    glTexEnvf( GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE );
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_GENERATE_MIPMAP, GL_TRUE); 

    // For debugging: 
    //     glTexImage2D(GL_PROXY_TEXTURE_2D, 0, GL_RGBA16F_ARB, 
    //                  m_image.cols(), m_image.rows(), 0, 
    //                  GL_RGBA, GL_FLOAT, NULL);
    //     GLint width; 
    //     glGetTexLevelParameteriv(GL_PROXY_TEXTURE_2D, 0, GL_TEXTURE_WIDTH, &width); 
    //     std::cout << "WIDTH: " << width << "\n";
    //     if (width==0) { /* Can't use that texture */ }

    // Grab the appropriate pixel data from our image.
    vw::ImageView<vw::PixelRGBA<vw::float32> > block = crop(m_image, m_bboxes[i]);

    // Copy the texture data over into 16-bit floating point texture
    // memory.  This prevents the data from being clamped in the range
    // [0,1].
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA16F_ARB, 
                 block.cols(), block.rows(), 0, 
                 GL_RGBA, GL_FLOAT, &(block(0,0)) );
  }
  std::cout << "done.\n";
}


void GlPreviewWidget::initializeGL() {  

  // Set up the texture mode to replace (rather than blend...)
  glShadeModel(GL_FLAT);
  
  // This is required for supporting Alpha in textures.
  glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glEnable (GL_BLEND); 

  // Set up the fragment shader.
  const char* fragment_prog_ptr = g_FRAGMENT_PROGRAM.c_str();

  // For debugging:
  //  std::cout << "***\n" << std::string(fragment_prog_ptr) << "***\n";
  
  GLuint m_fragment_shader = glCreateShader(GL_FRAGMENT_SHADER);
  glShaderSource(m_fragment_shader, 1, &fragment_prog_ptr, NULL);
  glCompileShader(m_fragment_shader);
  print_shader_info_log(m_fragment_shader);
  
  m_glsl_program = glCreateProgram();
  glAttachShader(m_glsl_program, m_fragment_shader);
  glLinkProgram(m_glsl_program);
  print_program_info_log(m_glsl_program);
}

void GlPreviewWidget::resizeGL(int width, int height) {
  m_viewport_width = width;
  m_viewport_height = height;
  size_to_fit();
}


// --------------------------------------------------------------
//             GlPreviewWidget Event Handlers
// --------------------------------------------------------------

void GlPreviewWidget::paintEvent(QPaintEvent * /* event */) { 
  QPainter painter(this);
  drawImage();
  if (m_show_legend)
    drawLegend(&painter);
}

void GlPreviewWidget::mousePressEvent(QMouseEvent *event) { 
  m_show_legend = true;
  grabKeyboard();
  lastPos = event->pos();
  updateCurrentMousePosition();
}

void GlPreviewWidget::mouseMoveEvent(QMouseEvent *event) {
  // Left mouse button moves the image around
  if (event->buttons() & Qt::LeftButton) {
    float x_diff = float(event->x() - lastPos.x()) / m_viewport_width;
    float y_diff = float(event->y() - lastPos.y()) / m_viewport_height;
    float ticks;

    std::ostringstream s; 
    switch (m_adjust_mode) {

    case TransformAdjustment:
      m_current_viewport.min().x() -= x_diff * m_current_viewport.width();
      m_current_viewport.min().y() -= y_diff * m_current_viewport.height();
      m_current_viewport.max().x() -= x_diff * m_current_viewport.width();
      m_current_viewport.max().y() -= y_diff * m_current_viewport.height();
      break;

    case GainAdjustment:
      // The number '5' below adjusts the sensitivity.
      ticks = pow(2, 5 * x_diff);
      if (m_gain * ticks > 1e-8 && m_gain * ticks < 1e8)
        m_gain *= ticks;
      s << "Gain: " << (log(m_gain)/log(2)) << " f-stops\n";
      m_legend_status = s.str();
      break;

    case OffsetAdjustment:
      m_offset += x_diff * (m_image_max - m_image_min);
      s << "Offset: " << m_offset << "\n";
      m_legend_status = s.str();
      break;

    case GammaAdjustment:
      // The number '5.0' below adjust the sensitivity.
      ticks = pow(2, x_diff * 5.0);
      if (m_gamma * ticks > 0.01 && m_gamma * ticks < 10.0)
        m_gamma *= ticks;
      s << "Gamma: " << m_gamma << "\n";
      m_legend_status = s.str();
      break;
    }

  } else if (event->buttons() & Qt::RightButton) {
    m_gain += GLfloat(event->x() - lastPos.x()) / m_viewport_width *10;
  } 

  // Regardless, we store the current position for the text legend.
  lastPos = event->pos();
  updateCurrentMousePosition();
  update();
}

void GlPreviewWidget::mouseDoubleClickEvent(QMouseEvent * /*event*/) {
  m_draw_texture = !m_draw_texture;
  update();
}

void GlPreviewWidget::wheelEvent(QWheelEvent *event) {
  int num_degrees = event->delta() / 8;
  float num_ticks = num_degrees / 15;

  float scale = pow(2,num_ticks/5);
  zoom(scale);

  m_show_legend = true;
  grabKeyboard();
  lastPos = event->pos();
  updateCurrentMousePosition();
}


void GlPreviewWidget::enterEvent(QEvent */*event*/) {
  m_show_legend = true;
  grabKeyboard();
  update();
}

void GlPreviewWidget::leaveEvent(QEvent */*event*/) {
  m_show_legend = false;
  releaseKeyboard();
  update();
}

void GlPreviewWidget::keyPressEvent(QKeyEvent *event) {

  std::ostringstream s; 
  
  switch (event->key()) {
  case Qt::Key_Plus:   // Zoom inxo
    zoom(2.0);
    break;
  case Qt::Key_Minus:  // Zoom out
    zoom(0.5);
    break;
  case Qt::Key_F:  // Size to fit
    size_to_fit();
    break;
  case Qt::Key_N:  // Toggle bilinear/nearest neighbor interp
    m_bilinear_filter = !m_bilinear_filter;
    update();
    break;
  case Qt::Key_C:  // Activate colormap
    m_use_colormap = !m_use_colormap;
    update();
    break;
  case Qt::Key_R:  // Normalize the image
    normalize_image();
    update();
    break;
  case Qt::Key_G:  // Gain adjustment mode
    if (m_adjust_mode == GainAdjustment) {
      m_adjust_mode = TransformAdjustment;
      m_legend_status = "";
    } else {
      m_adjust_mode = GainAdjustment;
      s << "Gain: " << log(m_gain)/log(2) << " f-stops\n";
      m_legend_status = s.str();
    }
    update();
    break;
  case Qt::Key_O:  // Offset adjustment mode
    if (m_adjust_mode == OffsetAdjustment) {
      m_adjust_mode = TransformAdjustment;
      m_legend_status = "";
    } else {
      m_adjust_mode = OffsetAdjustment;
      s << "Offset: " << m_offset;
      m_legend_status = s.str();
    }
    update();
    break;
  case Qt::Key_V:  // Gamma adjustment mode
    if (m_adjust_mode == GammaAdjustment) {
      m_adjust_mode = TransformAdjustment;
      m_legend_status = "";
    } else {
      m_adjust_mode = GammaAdjustment;
      s << "Gamma: " << m_gamma;
      m_legend_status = s.str();
    }
    update();
    break;
  case Qt::Key_1:  // Display red channel only
    m_display_channel = DisplayR;
    update();
    break;
  case Qt::Key_2:  // Display green channel only
    m_display_channel = DisplayG;
    update();
    break;
  case Qt::Key_3:  // Display blue channel only
    m_display_channel = DisplayB;
    update();
    break;
  case Qt::Key_4:  // Display alpha channel only
    m_display_channel = DisplayA;
    update();
    break;
  case Qt::Key_0:  // Display all color channels
    m_display_channel = DisplayRGBA;
    update();
    break;
  default: 
    QWidget::keyPressEvent(event);
  }
}

void GlPreviewWidget::add_crosshairs(std::list<Vector2> const& points, Vector3 const& color) {
  m_crosshairs.push_back(PointList(points, color));
  update();
}

void GlPreviewWidget::clear_crosshairs() {
  m_crosshairs.clear(); 
  update();
}
