// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// \file vwv_MainWindow.cc
///
/// The Vision Workbench image viewer main window class. 
///
#include <QtGui>
#include <vw/config.h>
#include <vw/gui/vwv_MainWindow.h>
#include <vw/gui/vwv_GlPreviewWidget.h>

#include <vw/Image/MaskViews.h>
#include <vw/FileIO/DiskImageView.h>
#include <vw/Image/Statistics.h>
#include <vw/Image/PixelMask.h>

#include <sstream>

MainWindow::MainWindow(std::string input_filename, float nodata_value, bool do_normalize, po::variables_map const& vm) :
  m_filename(input_filename), m_nodata_value(nodata_value), m_vm(vm) {

  // Set up the basic layout of the window and its menus
  create_actions();
  create_menus();
  create_status_bar();
  
  // Set the window title and add tabs
  std::string window_title = "Vision Workbench Viewer : " + m_filename;
  this->setWindowTitle(window_title.c_str());

  // Set up GlPreviewWidget
  m_preview_widget = new GlPreviewWidget(this);
  setCentralWidget(m_preview_widget);

  vw::DiskImageView<float> image(m_filename);
  
  // Set the nodata value
  if (m_vm.count("nodata-value")) {
    m_preview_widget->set_nodata_value(m_nodata_value);
  }

  // If normalization is requested, compute the min/max
  float lo = 0.0, hi = 1.0;
  if (do_normalize) {
    std::cout << "\t--> Computing data range: " << std::flush;
    if (m_vm.count("nodata-value")) {
      vw::min_max_channel_values(vw::create_mask(image,m_nodata_value), lo, hi);
    } else {
      vw::min_max_channel_values(image, lo, hi);
    }
    std::cout << "[L: " << lo << " H: " << hi << "]\n";
  }
  m_preview_widget->set_data_range(lo, hi);

  // Pass the filename along to the preview widget for display.
  m_preview_widget->set_image_from_file(input_filename);

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
