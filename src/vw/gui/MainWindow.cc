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


/// \file vwv_MainWindow.cc
///
/// The Vision Workbench image viewer main window class.
///
#include <QtGui>
#include <vw/gui/MainWindow.h>
#include <vw/gui/MainWidget.h>
#include <vw/gui/Utils.h>
using namespace vw::gui;

#include <vw/config.h>
#include <vw/Image/MaskViews.h>
#include <vw/FileIO/DiskImageView.h>
#include <vw/FileIO/DiskImageResource.h>
#include <vw/Image/Statistics.h>
#include <vw/Image/PixelMask.h>

#include <sstream>
namespace po = boost::program_options;

MainWindow::MainWindow(std::vector<std::string> const& images,
                       std::string const& geom,
                       float nodata_value,
                       int transaction_id,
                       bool /*do_normalize*/,
                       po::variables_map const& vm) :
  m_images(images), m_nodata_value(nodata_value), m_vm(vm) {

  // Set up the basic layout of the window and its menus
  create_actions();
  create_menus();

  // Set the window title and add tabs
  std::string window_title = "Vision Workbench Viewer";
  this->setWindowTitle(window_title.c_str());

  // Set up MainWidget
  m_preview_widget = new MainWidget(this, images,
                                         transaction_id);
  setCentralWidget(m_preview_widget);

  int windowWidX, windowWidY;
  extractWindowDims(geom, windowWidX, windowWidY);
  resize(windowWidX, windowWidY);
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

void MainWindow::about() {
  std::ostringstream about_text;
  about_text << "<h3>Vision Workbench Image Viewer (vwv)</h3>"
             << "<p>Version " << VW_PACKAGE_VERSION << "</p>"
             << "<p>Copyright &copy; 2014 NASA Ames Research Center</p>";
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
