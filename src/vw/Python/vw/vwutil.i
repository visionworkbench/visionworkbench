// __BEGIN_LICENSE__
// Copyright (C) 2006-2011 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// \file util.i
///
/// Defines utility tools used by the Python bindings.
///

%define HANDLE_VW_EXCEPTIONS(function)
%exception function {
  try {
    $action
  }
  catch (const vw::Exception& e) {
    // If there is an existing Python exception, we just pass it along.
    // It's a shame we that lose much of the stack trace on the C++ side,
    // but life is suffering.
    if( ! PyErr_Occurred() ) {
      // Otherwise, we construct a generic complaint.
      PyErr_Format( PyExc_RuntimeError, "Vision Workbench exception: %s", e.what() );
    }
    goto fail;
  }
}
%enddef

%{

// A deleter object for use with boost::shared_ptr that keeps track
// of the Python object corresponding to whatever C++ object it's
// used with, and decrements the Python object's refcount instead
// of deleting the C++ object directly.  If the optional second
// argument to the constructor is true, we grab a new referece to
// the object; otherwise, we steal the caller's reference.
class DecrefDeleter {
  PyObject *m_obj;
public:
  DecrefDeleter( PyObject *obj, bool incref = false ) : m_obj(obj) {
    if( incref ) Py_INCREF(obj);
  }
  template <class T> void operator()(T) {
    Py_DECREF(m_obj);
  }
};

%}
