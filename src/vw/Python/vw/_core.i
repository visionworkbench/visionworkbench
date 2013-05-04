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


%module core

%include "vwutil.i"

%{
#include <vw/Core/Log.h>
#include <vw/Core/Debugging.h>
#include <vw/Core/ProgressCallback.h>
%}

namespace vw {

  enum MessageLevel {
    ErrorMessage = 0,
    WarningMessage = 10,
    InfoMessage = 20,
    DebugMessage = 30,
    VerboseDebugMessage = 40
  };

  void set_debug_level( int level );

  /// The base class for progress monitoring.
  class ProgressCallback {
  public:
    ProgressCallback();
    virtual ~ProgressCallback();

    virtual void report_progress(double progress) const;
    virtual void report_incremental_progress(double incremental_progress) const;

    virtual void report_aborted(std::string why="") const;
    virtual void report_finished() const;

    double progress() const;
  };

  class TerminalProgressCallback : public ProgressCallback {
  public:
    TerminalProgressCallback( MessageLevel level = InfoMessage, std::string pre_progress_text = "" );
    virtual ~TerminalProgressCallback();

    void set_progress_text( std::string const& text );

    virtual void report_progress(double progress) const;
    virtual void report_incremental_progress(double incremental_progress) const;

    virtual void report_aborted(std::string why="") const;
    virtual void report_finished() const;

    void print_progress() const;

    std::string pre_progress_text() const;
    MessageLevel message_level() const;
  };

} // namespace vw

// Okay there's GOT to be a better way to do this.
%typemap(in) PyObject *progress_func {
  $1 = $input;
}
%typemap(in) PyObject *finished_func {
  $1 = $input;
}
%typemap(in) PyObject *aborted_func {
  $1 = $input;
}

HANDLE_VW_EXCEPTIONS(vw::PythonProgressCallback::report_progress)
HANDLE_VW_EXCEPTIONS(vw::PythonProgressCallback::report_finished)
HANDLE_VW_EXCEPTIONS(vw::PythonProgressCallback::report_aborted)

%inline %{
namespace vw {

  class PythonProgressCallback : public ProgressCallback {
    boost::shared_ptr<PyObject> m_progress_func;
    boost::shared_ptr<PyObject> m_aborted_func;
    boost::shared_ptr<PyObject> m_finished_func;
    mutable double m_last_reported;

    void call_progress_func(double progress) const {
      if( progress == m_last_reported ) return;
      m_last_reported = progress;
      if( ! m_progress_func ) return;
      PyEval_CallFunction( m_progress_func.get(), "(d)", progress );
      if( PyErr_Occurred() ) vw_throw( vw::Exception() );
    }

  public:
    PythonProgressCallback( PyObject *progress_func, PyObject *finished_func, PyObject *aborted_func )
      : m_last_reported(0)
    {
      if( progress_func != Py_None ) m_progress_func.reset( progress_func, DecrefDeleter(progress_func,true) );
      if( finished_func != Py_None ) m_finished_func.reset( finished_func, DecrefDeleter(finished_func,true) );
      if( aborted_func  != Py_None ) m_aborted_func.reset(  aborted_func,  DecrefDeleter(aborted_func, true) );
    }

    void report_progress(double progress) const {
      ProgressCallback::report_progress(progress);
      call_progress_func(progress);
    }
    void report_incremental_progress(double incremental_progress) const {
      ProgressCallback::report_incremental_progress(incremental_progress);
      call_progress_func(progress());
    }

    void report_finished() const {
      ProgressCallback::report_finished();
      if( ! m_finished_func ) return;
      PyEval_CallFunction( m_finished_func.get(), "()" );
      if( PyErr_Occurred() ) vw_throw( vw::Exception() );
    }
    void report_aborted(std::string why="") const {
      ProgressCallback::report_aborted(why);
      if( ! m_aborted_func ) return;
      PyEval_CallFunction( m_aborted_func.get(), "(s)", why.c_str() );
      if( PyErr_Occurred() ) vw_throw( vw::Exception() );
    }

  };

} // namespace vw
%}

// Apparently you can't %pythoncode inside an %inline block, so we do this out here.
%pythoncode {
  PythonProgressCallback._old_init = PythonProgressCallback.__init__
  def _PythonProgressCallback__init__(self,progress=None,finished=None,aborted=None):
    self._old_init(progress,finished,aborted)
  PythonProgressCallback.__init__ = _PythonProgressCallback__init__
}
