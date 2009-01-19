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

/// \file vwv_MainWindow.h
///
/// The Vision Workbench image viewer main window class. 
///
#ifndef __VWV_MAINWINDOW_H__
#define __VWV_MAINWINDOW_H__

#include <QMainWindow>
#include <string>

// Boost
#include <boost/program_options.hpp>
using namespace boost;
namespace po = boost::program_options;

class QAction;
class QLabel;
class QTabWidget;
class GlPreviewWidget;

class MainWindow : public QMainWindow {
  Q_OBJECT

  std::string m_filename;
  float m_nodata_value;
  po::variables_map const& m_vm;

public:
  MainWindow(std::string filename, float nodata_value, bool do_normalize, po::variables_map const& vm);
  virtual ~MainWindow() {}

private slots:
  void about();
  void update_status_bar(std::string const& s);

protected:
  void keyPressEvent(QKeyEvent *event);

private:
  void create_actions();
  void create_menus();
  void create_status_bar();

  GlPreviewWidget *m_preview_widget;

  QMenu *file_menu;
  QMenu *edit_menu;
  QMenu *help_menu;

  QLabel *status_label;
  QAction *about_action;
  QAction *exit_action;
};

#endif // __VWV_MAINWINDOW_H__
