// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
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

// Forward declarations
class QAction;
class QLabel;
class QTabWidget;

namespace vw {
namespace gui {

  class GlPreviewWidget;

  class MainWindow : public QMainWindow {
    Q_OBJECT

    std::string m_filename;
    float m_nodata_value;
    boost::program_options::variables_map const& m_vm;

  public:
    MainWindow(std::string filename, float nodata_value, int transaction_id, bool do_normalize, boost::program_options::variables_map const& vm);
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

}} // namespace vw::gui

#endif // __VWV_MAINWINDOW_H__
