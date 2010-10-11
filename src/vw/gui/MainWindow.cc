// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// \file vwv_MainWindow.cc
///
/// The Vision Workbench image viewer main window class.
///
#include <QtGui>
#include <vw/config.h>
#include <vw/gui/MainWindow.h>
#include <vw/gui/GlPreviewWidget.h>
using namespace vw::gui;

#include <vw/Image/MaskViews.h>
#include <vw/FileIO/DiskImageView.h>
#include <vw/FileIO/DiskImageResource.h>
#include <vw/Image/Statistics.h>
#include <vw/Image/PixelMask.h>

#include <sstream>
namespace po = boost::program_options;

MainWindow::MainWindow(std::string input_filename,
                       float nodata_value,
                       int transaction_id,
                       bool /*do_normalize*/,
                       po::variables_map const& vm) :
  m_filename(input_filename), m_nodata_value(nodata_value), m_vm(vm) {

  // Set up the basic layout of the window and its menus
  create_actions();
  create_menus();
  create_status_bar();

  // Set the window title and add tabs
  std::string window_title = "Vision Workbench Viewer : " + m_filename;
  this->setWindowTitle(window_title.c_str());

  // Set up OpenGL context parameters
  QGLFormat gl_frmt = QGLFormat::defaultFormat();
  gl_frmt.setSampleBuffers(true);
  gl_frmt.setDoubleBuffer(true);
  gl_frmt.setSwapInterval(1);

  // Set up GlPreviewWidget
  m_preview_widget = new GlPreviewWidget(this, input_filename, gl_frmt, transaction_id);
  setCentralWidget(m_preview_widget);

  // // Set the nodata value
  // vw::DiskImageResource *rsrc = vw::DiskImageResource::open(input_filename);
  // if (m_vm.count("nodata-value")) {
  //   std::cout << "\t--> User specified nodata value: " << m_nodata_value << ".\n";
  //   m_preview_widget->set_nodata_value(m_nodata_value);
  // } else if ( rsrc->has_nodata_value() ) {
  //   m_nodata_value = rsrc->nodata_value();
  //   std::cout << "\t--> Extracted nodata value from file: " << m_nodata_value << ".\n";
  //   m_preview_widget->set_nodata_value(m_nodata_value);
  // }
  // delete rsrc;

  // if (do_normalize) {
  //   m_preview_widget->normalize();
  // }

  // Maximize the main window
  this->showMaximized();
}

//********************************************************************
//                      MAIN WINDOW SETUP
//********************************************************************

void MainWindow::create_actions() {

  // The About Box
  about_action = new QAction(tr("About VWV"), this);
  about_action->setStatusTip(tr("Show the Vision Workbench Viewere about box."));
  connect(about_action, SIGNAL(triggered()), this, SLOT(about()));

  // Exit or Quit
  exit_action = new QAction(tr("E&xit"), this);
  exit_action->setShortcut(tr("Ctrl+Q"));
  exit_action->setStatusTip(tr("Exit the application"));
  connect(exit_action, SIGNAL(triggered()), this, SLOT(close()));
}

void MainWindow::create_menus() {

  // File Menu
  file_menu = menuBar()->addMenu(tr("&File"));
  file_menu->addAction(exit_action);

  // Edit Menu
  edit_menu = menuBar()->addMenu(tr("&Edit"));

  // Help menu
  menuBar()->addSeparator();
  help_menu = menuBar()->addMenu(tr("&Help"));
  help_menu->addAction(about_action);
}

void MainWindow::create_status_bar() {
  status_label = new QLabel("");
  status_label->setAlignment(Qt::AlignHCenter);
  statusBar()->addWidget(status_label);

  // WARNING: Memory leak as currently written.  Fix this somewhow...
  //   GuiProgressCallback *clbk = new GuiProgressCallback(this, "Updating: ");
  //   statusBar()->addWidget(clbk->widget());
}

void MainWindow::update_status_bar(std::string const& s) {
  status_label->setText(QString(s.c_str()));
}

void MainWindow::about() {
  std::ostringstream about_text;
  about_text << "<h3>Vision Workbench Image Viewer (vwv)</h3>"
             << "<p>Version " << VW_PACKAGE_VERSION << "</p>"
             << "<p>Copyright &copy; 2008 NASA Ames Research Center</p>";
  QMessageBox::about(this, tr("About Vision Workbench Viewer"),
                     tr(about_text.str().c_str()));

}

void MainWindow::keyPressEvent(QKeyEvent *event) {

  std::ostringstream s;

  switch (event->key()) {
  case Qt::Key_Escape:  // Quit
    close();
    break;
  }
}
