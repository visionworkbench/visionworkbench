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

#ifdef __linux__
// This is required to get prototypes, according to the opengl linux abi
#define GL_GLEXT_PROTOTYPES 1
#endif

#ifdef __APPLE__
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#else // Linux
#include <GL/gl.h>
#include <GL/glu.h>
#endif

// Qt
#include <QtGui>

// Vision Workbench
#include <vw/Image.h>
#include <vw/FileIO.h>
#include <vw/gui/TileGenerator.h>
#include <vw/gui/GlPreviewWidget.h>
using namespace vw;
using namespace vw::gui;

const std::string g_FRAGMENT_PROGRAM =
"uniform sampler2D tex;                             \n"
"                                                   \n"
"uniform double gain;\n"
"uniform double offset;\n"
"uniform double gamma;\n"
"uniform double nodata;\n"
"uniform int use_nodata;\n"
"uniform int display_channel;  // 0 = RGBA, 1 = R, 2 = G, 3 = B, 4 = A\n"
"uniform int colorize_display;\n"
"uniform int hillshade_display;\n"
"uniform int current_level;\n"
"\n"
"void main() { \n"
"  bool alpha_selected = false; \n"
"  vec4 g = vec4(gain, gain, gain, 1.0);\n"
"  vec4 o = vec4(offset, offset, offset, 0.0);\n"
"  vec4 v = vec4(gamma, gamma, gamma, 1.0);\n"
"\n"
"  vec4 src_tex = texture2D(tex,gl_TexCoord[0].st);\n"
"  vec4 selected_tex = src_tex;\n"
"  if ( hillshade_display > 0 ) {\n"
"    alpha_selected = true; \n"
"    vec4 right = texture2D(tex,gl_TexCoord[0].st+vec2(0.01,0.0));\n"
"    vec4 bot   = texture2D(tex,gl_TexCoord[0].st+vec2(0.0,0.01));\n"
"    vec3 normal = vec3(right.y - src_tex.y, bot.y - src_tex.y, (double(hillshade_display) * double(hillshade_display)) / double(current_level*current_level) );\n"
"    vec3 light = vec3(0.5774,0.5774,0.5774);\n"
"    double intense = dot(light,normalize(normal));\n"
"    selected_tex = vec4(intense,intense,intense,src_tex[3]);\n"
"  } else if (display_channel == 1) {\n"
"    selected_tex = vec4(src_tex[0],src_tex[0],src_tex[0],1.0);\n"
"  } else if (display_channel == 2) {\n"
"    selected_tex = vec4(src_tex[1],src_tex[1],src_tex[1],1.0);\n"
"  } else if (display_channel == 3) {\n"
"    selected_tex = vec4(src_tex[2],src_tex[2],src_tex[2],1.0);\n"
"  } else if (display_channel == 4) {\n"
"    alpha_selected = true;            \n"
"    selected_tex = vec4(src_tex[3],src_tex[3],src_tex[3],1.0);\n"
"  } \n"
"  \n"
"  if (!alpha_selected) {\n"
"    vec4 final_tex = pow((g*selected_tex)+o, v);\n"
"    if (use_nodata == 1) { \n"
"      if (src_tex[0] == nodata) { \n"
"        final_tex = vec4(0.0,0,0,0.0);\n"
"      } \n"
"    }\n"
"    gl_FragColor = final_tex; \n"
"  } else {\n"
"    gl_FragColor = selected_tex;\n"
"  }\n"
"}\n";

namespace {
  vw::Vector2 QPoint2Vec(QPoint const& qpt) {
    return vw::Vector2(qpt.x(), qpt.y());
  }
  QPoint Vec2QPoint(vw::Vector2 const& V) {
    return QPoint(round(V.x()), round(V.y()));
  }
}

// --------------------------------------------------------------
//                       GLSL DEBUGGING
// --------------------------------------------------------------
void print_shader_info_log(GLuint obj)
{
  GLint infologLength = 0;
  GLint charsWritten  = 0;
  char *infoLog;

  //glGetShaderiv(obj, GL_INFO_LOG_LENGTH,&infologLength);

  if (infologLength > 1) {
    infoLog = (char *)malloc(infologLength);
//     glGetShaderInfoLog(obj, infologLength, &charsWritten, infoLog);
    std::ostringstream err;
    err << "<h4>An error occured while compiling the GLSL shader:</h4><p><h5><tt>" << infoLog << "</tt></h5>";
    QMessageBox::critical(0, "GLSL Shader Error",
                            err.str().c_str());
    free(infoLog);
  }
}

void print_program_info_log(GLuint obj)
{
  GLint infologLength = 0;
  GLint charsWritten  = 0;
  char *infoLog;

//   glGetProgramiv(obj, GL_INFO_LOG_LENGTH,&infologLength);

  if (infologLength > 1) {
    infoLog = (char *)malloc(infologLength);
//     glGetProgramInfoLog(obj, infologLength, &charsWritten, infoLog);
    std::ostringstream err;
    err << "<h4>An error occured while linking the GLSL program:</h4><p><h5><tt>" << infoLog << "</tt></h5>";
    QMessageBox::critical(0, "GLSL Program Error",
                          err.str().c_str());
    printf("%s\n",infoLog);
    free(infoLog);
  }
}

void vw::gui::check_gl_errors()
{
  GLenum errCode;
  const GLubyte *errStr;
//   if ( ( errCode = glGetError( ) ) != GL_NO_ERROR ) {
//     errStr = gluErrorString( errCode );
//     vw_out(ErrorMessage) << "OpenGL ERROR (" << int(errCode) << "): " << errStr << "\n";
//   }
}

// --------------------------------------------------------------
//               GlPreviewWidget Public Methods
// --------------------------------------------------------------

GlPreviewWidget::GlPreviewWidget(QWidget *parent,
                                 std::vector<std::string> const& images,
                                 //QGLFormat const& frmt,
                                 int transaction_id)
  : QWidget(parent)
{

  // Verify that our OpenGL formatting options stuck
//   if (!QGLFormat::hasOpenGL()) {
//     vw::vw_out() << "This system has no OpenGL support.\nExiting\n\n";
//     exit(1);
//   }

//   check_gl_errors();

//   if (!format().sampleBuffers())
//     std::cout << "\n\nCould not activate FSAA; results will be suboptimal\n\n";
//   if (!format().doubleBuffer())
//     std::cout << "\n\nCould not set double buffering; results will be suboptimal\n\n";

  // Set default values
  m_nodata_value = 0;
  m_use_nodata = 0;
  m_image_min = 0;
  m_image_max = 1.0;

  // Set some reasonable defaults
  m_show_legend = false;
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
  
  m_tile_generator = TileGenerator::create(images[0]);
  m_gl_texture_cache.reset( new GlTextureCache(m_tile_generator) );
  size_to_fit();
}


GlPreviewWidget::~GlPreviewWidget() {
  m_gl_texture_cache.reset();
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
  m_show_legend = false;

  // Check to make sure we haven't hit our zoom limits...
  if (m_current_viewport.width()/scale > 1.0 &&
      m_current_viewport.height()/scale > 1.0 &&
      m_current_viewport.width()/scale < 20*m_images_box.width() &&
      m_current_viewport.height()/scale < 20*m_images_box.height()) {
    m_current_viewport = (m_current_viewport - m_curr_world_pos) / scale + m_curr_world_pos;
  }
  update(); // will call paintEvent()
}

// --------------------------------------------------------------
//             GlPreviewWidget Setup Methods
// --------------------------------------------------------------

GLuint GlPreviewWidget::allocate_texture(boost::shared_ptr<SrcImageResource> tile) {
  GLuint texture_id;

  //makeCurrent();
//   glEnable( GL_TEXTURE_2D );
//   glGenTextures(1, &(texture_id) );
//   glBindTexture(GL_TEXTURE_2D, texture_id);
//   glTexEnvf( GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE );
//   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
//   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
//   glTexParameteri(GL_TEXTURE_2D, GL_GENERATE_MIPMAP, GL_TRUE);

  vw_out(VerboseDebugMessage, "gui") << "Allocating texture for "
                                     << tile->cols() << "x" << tile->rows() << " sized tile.\n";

  // std::cout << "This image has " << tile->channels() << " channels\n";
  // std::cout << "           ond " << tile->channel_type() << " channel type\n";

  // We save VRAM by copying the image data over in its native
  // number of channels.
  GLuint texture_pixel_type = GL_RGBA32F_ARB;
  GLuint source_pixel_type = GL_RGBA;

  if (tile->channels() == 1) {
    texture_pixel_type = GL_LUMINANCE32F_ARB;
    source_pixel_type = GL_LUMINANCE;
  } else if (tile->channels() == 2) {
    texture_pixel_type = GL_LUMINANCE_ALPHA32F_ARB;
    source_pixel_type = GL_LUMINANCE_ALPHA;
  } else if (tile->channels() == 3) {
    texture_pixel_type = GL_RGB32F_ARB;
    source_pixel_type = GL_RGB;
  } else if (tile->channels() == 4) {
    texture_pixel_type = GL_RGBA32F_ARB;
    source_pixel_type = GL_RGBA;
  } else {
    vw_throw(vw::ArgumentErr() << "GlPreviewWidget: allocate_texture() failed."
             << " Unsupported number of channels (" << tile->channels() << ").");
  }

  // Set the GL channel type for source data.
  GLuint source_channel_type = GL_FLOAT;

  switch (tile->channel_type()) {
  case vw::VW_CHANNEL_UINT8:
    source_channel_type = GL_UNSIGNED_BYTE;
    break;
  case vw::VW_CHANNEL_INT8:
    source_channel_type = GL_BYTE;
    break;
  case vw::VW_CHANNEL_UINT16:
    source_channel_type = GL_UNSIGNED_SHORT;
    break;
  case vw::VW_CHANNEL_INT16:
    source_channel_type = GL_SHORT;
    break;
  case vw::VW_CHANNEL_UINT32:
    source_channel_type = GL_UNSIGNED_INT;
    break;
  case vw::VW_CHANNEL_INT32:
    source_channel_type = GL_INT;
    break;
  case vw::VW_CHANNEL_FLOAT32:
    source_channel_type = GL_FLOAT;
    break;
  default:
    vw_throw(vw::ArgumentErr() << "GlPreviewWidget: allocate_texture() failed."
             << " Unsupported channel type (" << tile->channel_type() << ").");
  }

//   glTexImage2D(GL_TEXTURE_2D, 0, texture_pixel_type,
//                tile->cols(), tile->rows(), 0,
//                source_pixel_type, source_channel_type, tile->native_ptr().get() );

//   glBindTexture(GL_TEXTURE_2D, 0);
//   glDisable( GL_TEXTURE_2D );
  return texture_id;
}

void GlPreviewWidget::deallocate_texture(GLuint texture_id) {
  //makeCurrent();
//   glEnable( GL_TEXTURE_2D );
//   glDeleteTextures(1,&(texture_id));
//   glDisable( GL_TEXTURE_2D );
}

// ---------------------------------------------------------

void GlPreviewWidget::initializeGL() {

  // Set up the texture mode to replace (rather than blend...)
//   glShadeModel(GL_FLAT);

  // Set up the fragment shader.
  const char* fragment_prog_ptr = g_FRAGMENT_PROGRAM.c_str();
  GLuint m_fragment_shader;// = glCreateShader(GL_FRAGMENT_SHADER);
//   glShaderSource(m_fragment_shader, 1, &fragment_prog_ptr, NULL);
//   glCompileShader(m_fragment_shader);
//   print_shader_info_log(m_fragment_shader);

//   m_glsl_program = glCreateProgram();
//   glAttachShader(m_glsl_program, m_fragment_shader);
//   glLinkProgram(m_glsl_program);
//   print_program_info_log(m_glsl_program);

  // Now that GL is setup, we can start the Qt Timer
  m_timer = new QTimer(this);
  connect(m_timer, SIGNAL(timeout()), this, SLOT(timer_callback()));
  m_timer->start(33);
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
    
  // Make this context current, and store the current OpenGL state
  // before we start to modify it.
  //makeCurrent();
//   glPushAttrib(GL_ALL_ATTRIB_BITS);

//   // Before we draw this frame, we will check to see whether there are
//   // any new texture to upload or delete from the texture cache.  If
//   // there are, we perform at least one of these operations.
//   this->process_allocation_requests();

//   // Activate our GLSL fragment program and set up the uniform
//   // variables in the shader
//   glUseProgram(m_glsl_program);
//   GLint gain_loc = glGetUniformLocation(m_glsl_program,"gain");
//   glUniform1f(gain_loc,m_gain);
//   GLint offset_loc = glGetUniformLocation(m_glsl_program,"offset");
//   glUniform1f(offset_loc,m_offset);
//   GLint gamma_loc = glGetUniformLocation(m_glsl_program,"gamma");
//   glUniform1f(gamma_loc,m_gamma);
//   GLint nodata_loc = glGetUniformLocation(m_glsl_program,"nodata");
//   glUniform1f(nodata_loc,m_nodata_value);
//   GLint use_nodata_loc = glGetUniformLocation(m_glsl_program,"use_nodata");
//   glUniform1i(use_nodata_loc,m_use_nodata);
//   GLint display_channel_loc = glGetUniformLocation(m_glsl_program,"display_channel");
//   glUniform1i(display_channel_loc,m_display_channel);
//   GLint colorize_display_loc = glGetUniformLocation(m_glsl_program,"colorize_display");
//   glUniform1i(colorize_display_loc,m_colorize_display);
//   GLint hillshade_display_loc = glGetUniformLocation(m_glsl_program,"hillshade_display");
//   glUniform1i(hillshade_display_loc,m_hillshade_display);
//   GLint current_level_loc = glGetUniformLocation(m_glsl_program,"current_level");
//   glUniform1i(current_level_loc,m_current_level);
//   glUseProgram(0);

//   // Set the background color and viewport.
//   //qglClearColor(QColor(0, 25, 50)); // Bluish-green background
//   //qglClearColor(QColor(0, 0, 0)); // Black background
//   glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
//   glViewport(0,0,m_window_width,m_window_height);

//   // Set up the orthographic view of the scene.  The exact extent of
//   // the view onto the scene depends on the current panning and zoom
//   // in the UI.
//   glMatrixMode(GL_PROJECTION);
//   glPushMatrix();
//   glLoadIdentity();
//   glOrtho(m_current_viewport.min().x(), m_current_viewport.max().x(),
//           -m_current_viewport.max().y(), -m_current_viewport.min().y(),
//           -1.0, 1.0);

//   // Set up the modelview matrix, and bind the image as the texture we
//   // are about to use.
//   glMatrixMode(GL_MODELVIEW);
//   glPushMatrix();
//   glLoadIdentity();

//   // This is required for supporting Alpha in textures.
//   glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
//   glEnable (GL_BLEND);

  // Compute the current level of detail (limit to a minimum tile
  // level of 0 and maximum size that depends on how many tiles the
  // file contains.)
//   Vector2i tile_size = m_tile_generator->tile_size();
  BBox2i image_bbox(0,0,m_images_box.width(),m_images_box.height());

//   int max_level = m_tile_generator->num_levels()-1;
//   int level = max_level - log(double(m_current_viewport.width()) / m_window_width) / log(2.0)+1;
//   if (level < 0) level = 0;
//   if (level > max_level) level = max_level;
//   m_current_level = level;

//   std::list<TileLocator> tiles = bbox_to_tiles(tile_size, m_current_viewport, level, max_level, m_current_transaction_id, m_exact_transaction_id_match);
//   std::list<TileLocator>::iterator tile_iter = tiles.begin();

//   while (tile_iter != tiles.end()) {
//     BBox2i texture_bbox = tile_to_bbox(tile_size, tile_iter->col, tile_iter->row, tile_iter->level, max_level);
//     texture_bbox.crop(image_bbox);
//     if (tile_iter->is_valid() && image_bbox.width() > 0 && image_bbox.height() > 0) {

//       // Fetch the texture out of the cache.  If the texture is not
//       // currently in the cache, a request for this texture will be
//       // generated so that it is available at some point in the
//       // future. will be generated if necessary. Note that this
//       // happens outside the m_gl_mutex to avoid deadlock.
//       GLuint texture_id = m_gl_texture_cache->get_texture_id(*tile_iter, this);

//       if (texture_id) {
// //         glUseProgram(m_glsl_program);

//         // Enable texturing and bind the texture ID
// //         glEnable( GL_TEXTURE_2D );
// //         glBindTexture( GL_TEXTURE_2D, texture_id );

//         // Set up bilinear or nearest neighbor filtering.
//         if (m_bilinear_filter) {
//           // When the texture area is small, bilinear filter the closest mipmap
// //           glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR );
// //           glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR );
//         } else {
//           // When the texture area is small, pick the nearest neighbor in the closest mipmap
// //           glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST );
// //           glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST );
//         }

//         // Draw the texture onto a quad.
//         //qglColor(Qt::white);
// //         glBegin(GL_QUADS);
// //         glTexCoord2d( 0.0 , 0.0);
// //         glVertex2d( texture_bbox.min().x() , -(texture_bbox.min().y()) );
// //         glTexCoord2d( 0.0 , 1 );
// //         glVertex2d( texture_bbox.min().x() , -(texture_bbox.max().y()) );
// //         glTexCoord2d( 1.0 , 1.0 );
// //         glVertex2d( texture_bbox.max().x() , -(texture_bbox.max().y()) );
// //         glTexCoord2d( 1.0 , 0.0 );
// //         glVertex2d( texture_bbox.max().x() , -(texture_bbox.min().y()) );
// //         glEnd();

// //         // Clean up
// //         glDisable( GL_TEXTURE_2D );
// //         glUseProgram(0);

// //         // Optional: draw a border around the texture
// //         if (m_show_tile_boundaries) {
// //           //qglColor(Qt::blue);
// //           glBegin(GL_LINES);
// //           glVertex2d( texture_bbox.min().x() , -(texture_bbox.min().y()) );
// //           glVertex2d( texture_bbox.min().x() , -(texture_bbox.max().y()) );

// //           glVertex2d( texture_bbox.min().x() , -(texture_bbox.max().y()) );
// //           glVertex2d( texture_bbox.max().x() , -(texture_bbox.max().y()) );

// //           glVertex2d( texture_bbox.max().x() , -(texture_bbox.max().y()) );
// //           glVertex2d( texture_bbox.max().x() , -(texture_bbox.min().y()) );

// //           glVertex2d( texture_bbox.max().x() , -(texture_bbox.min().y()) );
// //           glVertex2d( texture_bbox.min().x() , -(texture_bbox.min().y()) );
// //           glEnd();
// //         }


// //       } else {
// //         // If no texture is (yet) available, we draw a dark blue quad.
// //         glBegin(GL_QUADS);
// //         glColor3f(0.0,0.0,0.1);
// //         glVertex2d( texture_bbox.min().x() , -(texture_bbox.min().y()) );
// //         glVertex2d( texture_bbox.min().x() , -(texture_bbox.max().y()) );
// //         glVertex2d( texture_bbox.max().x() , -(texture_bbox.max().y()) );
// //         glVertex2d( texture_bbox.max().x() , -(texture_bbox.min().y()) );
// //         glEnd();

// //         // Optional: draw a border around the texture
// //         if (m_show_tile_boundaries) {
// //           //qglColor(Qt::red);
// //           glBegin(GL_LINES);
// //           glVertex2d( texture_bbox.min().x() , -(texture_bbox.min().y()) );
// //           glVertex2d( texture_bbox.min().x() , -(texture_bbox.max().y()) );

// //           glVertex2d( texture_bbox.min().x() , -(texture_bbox.max().y()) );
// //           glVertex2d( texture_bbox.max().x() , -(texture_bbox.max().y()) );

// //           glVertex2d( texture_bbox.max().x() , -(texture_bbox.max().y()) );
// //           glVertex2d( texture_bbox.max().x() , -(texture_bbox.min().y()) );

// //           glVertex2d( texture_bbox.max().x() , -(texture_bbox.min().y()) );
// //           glVertex2d( texture_bbox.min().x() , -(texture_bbox.min().y()) );
// //           glEnd();
// //       }
      
//       }
//     }
    
//     // Move onto the next tile
//     ++tile_iter;
//   }
//   //  std::cout << "\n";

//   // Restore the previous OpenGL state so that we don't trample on the
//   // QPainter elements of the window.
// //   glMatrixMode(GL_MODELVIEW);
// //   glPopMatrix();
// //   glMatrixMode(GL_PROJECTION);
// //   glPopMatrix();
// //   glPopAttrib();
}

void GlPreviewWidget::drawLegend(QPainter* paint) {

  // Extract the value for the pixel currently under the mouse
  if (m_curr_world_pos.x() >= 0 && m_curr_world_pos.x() < m_images_box.width() &&
      m_curr_world_pos.y() >= 0 && m_curr_world_pos.y() < m_images_box.height()) {


    // Note: we sample directly from the OpenGL texture buffer for
    // now, but eventually we would prefer to sample the actual value
    // from the tile_generator itself.
    double raw_pixels[4];
//     glReadPixels(m_curr_pixel_pos.x(),m_window_height-m_curr_pixel_pos.y(),
//                  1,1,GL_RGBA,GL_FLOAT,&raw_pixels);
    PixelRGBA<float32> pix_value(raw_pixels[0], raw_pixels[1], raw_pixels[2], raw_pixels[3]);
    //    PixelRGBA<float32> pix_value = m_last_pixel_sample;

    const char* pixel_name = vw::pixel_format_name(m_tile_generator->pixel_format());
    const char* channel_name = vw::channel_type_name(m_tile_generator->channel_type());
    int num_channels = vw::num_channels(m_tile_generator->pixel_format());

//     // Compute the tile location
//     int tile_x = m_curr_world_pos.x() /
//       pow(2,(m_tile_generator->num_levels()-1) - m_current_level) /
//       m_tile_generator->tile_size()[0];
//     int tile_y = m_curr_world_pos.y() /
//       pow(2,(m_tile_generator->num_levels()-1) - m_current_level) /
//       m_tile_generator->tile_size()[1];

    std::ostringstream pix_value_ostr;
//     pix_value_ostr << "Pos: ( " << m_curr_world_pos.x() << " " << m_curr_world_pos.y() << " )"
//                    << "[ " << tile_x << " " << tile_y << " ] --> Val: [ ";

    // The following code is very messy and should be replaced and/or
    // simplified once we have better control over whether automatic
    // conversions happen when reading in images using the FileIO
    // code.
    bool round;
    double scale_factor;
    switch (m_tile_generator->channel_type()) {
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
      vw_throw( ArgumentErr()
                << "vwv_GlPreviewWidget() : Unrecognized or unsupported channel type ("
                << m_tile_generator->channel_type() << ")." );
      return; // never reached
    }

    // We don't print out the channel value for the alpha channel for
    // those pixel types that have alpha.  This is messy too, and
    // should also be replaced by a better method...
    bool has_alpha = false;
    if (m_tile_generator->pixel_format() == vw::VW_PIXEL_GRAYA ||
        m_tile_generator->pixel_format() == vw::VW_PIXEL_RGBA ||
        m_tile_generator->pixel_format() == vw::VW_PIXEL_GRAY_MASKED ||
        m_tile_generator->pixel_format() == vw::VW_PIXEL_RGB_MASKED ||
        m_tile_generator->pixel_format() == vw::VW_PIXEL_XYZ_MASKED ||
        m_tile_generator->pixel_format() == vw::VW_PIXEL_LUV_MASKED ||
        m_tile_generator->pixel_format() == vw::VW_PIXEL_LAB_MASKED) {
      has_alpha = true;
      num_channels -= 1;
    }

    for (int i=0; i < num_channels; ++i) {
      if (round)
        pix_value_ostr << llroundf(pix_value[i] * scale_factor) << " ";
      else
        pix_value_ostr << (pix_value[i] * scale_factor) << " ";
    }


    // We don't print out the channel value for the alpha channel for
    // those pixel types that have alpha.  This is messy too, and
    // should also be replaced by a better method...
    if (has_alpha) {
      pix_value_ostr << "   " << (pix_value[3] * scale_factor);
    }
    pix_value_ostr << "]";

    const int Margin = 11;
    const int Padding = 6;
    QTextDocument textDocument;
    textDocument.setDefaultStyleSheet("* { color: #00FF00; font-family: courier, serif; font-size: 12 }");
    std::ostringstream legend_text;
    legend_text << "<p align=\"right\">" << m_legend_status << "<br>"
                << "[ " << m_images_box.width() << " x "
                << m_images_box.height() << " @ " << m_current_level << " ] (t_id = "
                << m_current_transaction_id << ")<br>"
                << "Pixel Format: " << pixel_name << "  Channel Type: " << channel_name << "<br>"
                << "Current Pixel Range: [ " << -m_offset/m_gain << " "
                << ( pow(1,1/m_gamma) - m_offset )/m_gain << " ]<br>"
                << pix_value_ostr.str() << "<br>"
                << "</p>";
    textDocument.setHtml(legend_text.str().c_str());
    textDocument.setTextWidth(textDocument.size().width());

    QRect rect(QPoint(0,0), textDocument.size().toSize() + QSize(2 * Padding, 2 * Padding));
    paint->translate(width() - rect.width() - Margin, height() - rect.height() - Margin);
    //     paint->setPen(QColor(255, 239, 239));
    //     paint->drawRect(rect);
    paint->translate(Padding, Padding);
    textDocument.drawContents(paint);
  }
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
  if (m_show_legend)
    drawLegend(&paint);
  //swapBuffers();
}

void GlPreviewWidget::mousePressEvent(QMouseEvent *event) {
  m_show_legend = true;
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
      m_legend_status = s.str();
      break;

    case OffsetAdjustment:
      m_offset = m_last_offset +
        (pow(100,fabs(x_diff))-1.0)*(x_diff > 0 ? 0.1 : -0.1);
      s << "Offset: " << m_offset << "\n";
      m_legend_status = s.str();
      break;

    case GammaAdjustment:
      m_gamma = m_last_gamma * pow(2.0,x_diff);
      s << "Gamma: " << m_gamma << "\n";
      m_legend_status = s.str();
      break;
    }
  }

  // Regardless, we store the current position for the text legend.
  updateCurrentMousePosition();
}

void GlPreviewWidget::mouseDoubleClickEvent(QMouseEvent *event) {
  m_curr_pixel_pos = QPoint2Vec(event->pos());
  updateCurrentMousePosition();
//   m_last_pixel_sample = m_tile_generator->sample(m_curr_world_pos.x(),
//                                                  m_curr_world_pos.y(),
//                                                  m_current_level,
//                                                  m_current_transaction_id);
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

  m_show_legend = true;
  m_curr_pixel_pos = QPoint2Vec(event->pos());
  updateCurrentMousePosition();
}


void GlPreviewWidget::enterEvent(QEvent */*event*/) {
  m_show_legend = true;
}

void GlPreviewWidget::leaveEvent(QEvent */*event*/) {
  m_show_legend = false;
}

void GlPreviewWidget::keyPressEvent(QKeyEvent *event) {

  std::ostringstream s;

  switch (event->key()) {
  case Qt::Key_Plus:   // Increase transaction id
  case Qt::Key_Equal:
    m_current_transaction_id++;
    m_gl_texture_cache->clear();
    break;
  case Qt::Key_Minus:  // Decrease transaction id
  case Qt::Key_Underscore:
    m_current_transaction_id--;
    m_gl_texture_cache->clear();
    break;
  case Qt::Key_Greater:   // Increase transaction id
    m_current_transaction_id+=100;
    m_gl_texture_cache->clear();
    break;
  case Qt::Key_Less:  // Decrease transaction id
    m_current_transaction_id-=100;
    m_gl_texture_cache->clear();
    break;
  case Qt::Key_E:  // Toggle exact transaction id match
    m_exact_transaction_id_match = !m_exact_transaction_id_match;
    m_gl_texture_cache->clear();
    break;
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
  case Qt::Key_R: // Reload the image
    m_gl_texture_cache->clear();
    break;
  case Qt::Key_N:  // Normalize the image
    this->normalize();
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

