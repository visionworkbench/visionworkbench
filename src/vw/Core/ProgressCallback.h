// __BEGIN_LICENSE__
// 
// Copyright (C) 2006 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration
// (NASA).  All Rights Reserved.
// 
// Copyright 2006 Carnegie Mellon University. All rights reserved.
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

/// \file ProgressCallback.h
/// 
/// A class for monitoring the progress of lengthy operations.
///
#ifndef __VW_CORE_PROGRESSCALLBACK_H__
#define __VW_CORE_PROGRESSCALLBACK_H__

#include <math.h>
#include <string>

#include <vw/Core/Debugging.h>
#include <vw/Core/Exception.h>

namespace vw {

  /// The base class for progress monitoring.
  class ProgressCallback {
  protected:
    // WARNING:  this flag may not be valid for some subclasses.  Always call abort_requested() and request_abort()
    bool m_abort_requested;
  public:
    ProgressCallback() : m_abort_requested( false ) {}

    // Reporting functions
    // Subclasses should reimplement where appropriate
    //
    // progress is from 0 (not done) to 1 (finished)
    // 
    virtual void report_progress(double /*progress*/) const {}
    virtual void report_aborted(std::string /*why*/="") const {}
    virtual void report_finished() const {}

    // Helper method which computes progress and calls report_progress
    void report_fractional_progress(double n, double total) const {
      report_progress(fabs(total) < 1e-30 ? 0 : n/total);
    }

    // Has an abort been requested?
    virtual bool abort_requested() const { return m_abort_requested; }

    // Throw vw::Aborted if abort has been requested
    virtual void abort_if_requested() const {
      if (abort_requested()) {
        report_aborted();
        vw_throw(Aborted());
      }
    }

    // Request abort
    virtual void request_abort() { m_abort_requested = true; }

    virtual ~ProgressCallback() {}
    static const ProgressCallback &dummy_instance();
  };


  /// Monitors the progress of a subtask of a task.
  class SubProgressCallback : public ProgressCallback {
  protected:
    const ProgressCallback &m_parent;
    const double m_from;
    const double m_to;
  public:
    SubProgressCallback(const ProgressCallback &parent,
                        double from, double to) :
      m_parent(parent), m_from(from), m_to(to) {}
    virtual void report_progress(double progress) const {
      double parent_progress = m_from + (m_to - m_from)*progress;
      m_parent.report_progress(parent_progress);
    }
    virtual void report_aborted(std::string why="") const {
      m_parent.report_aborted(why);
    }
    virtual bool abort_requested() const { return m_parent.abort_requested(); }
    virtual ~SubProgressCallback() {}
  };


  /// A progress monitor that prints a progress bar on STDOUT.
  class TerminalProgressCallback : public ProgressCallback {
    MessageLevel m_level;
    std::string m_pre_progress_text;
    mutable double m_last_reported_progress;
  public:
    TerminalProgressCallback( MessageLevel level = InfoMessage, std::string pre_progress_text = "" ) : 
      m_level(level), m_pre_progress_text(pre_progress_text), m_last_reported_progress(-1) {}
    virtual ~TerminalProgressCallback() {}

    void set_progress_text( std::string const& text ) {
      m_pre_progress_text = text;
    }

    virtual void report_progress(double progress) const;
    virtual void report_aborted(std::string why="") const;
    virtual void report_finished() const;
  };

} // namespace vw

#endif // __VW_CORE_FUNDAMENTALTYPES_H__
