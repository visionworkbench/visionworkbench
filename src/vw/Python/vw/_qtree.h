// __BEGIN_LICENSE__
// Copyright (C) 2006-2009 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// \file _qtree.h
///
/// Defines additional C++ types and functions used by the _qtree.i
/// interface definition file.
///

// SWIG doesn't currently support nested classes, so we spoof it
// into thinking that TileInfo is not nested.  To make it all
// work properly, it's best to spoof C++ too.
namespace vw {
namespace mosaic {
  typedef QuadTreeGenerator::TileInfo TileInfo;
} // namespace vw
} // namespace mosaic

// A wrapper for a user-supplied Python image_path_func callback function
std::string image_path_func( boost::shared_ptr<PyObject> const& pyfunc,
                             boost::shared_ptr<PyObject> const& qtree,
                             std::string const& name )
{
  std::string path;

  PyObject *path_obj = PyEval_CallFunction( pyfunc.get(), "(Os)", qtree.get(), name.c_str() );
  if( path_obj == NULL ) goto error;

  if( ! PyString_Check(path_obj) ) {
    PyErr_Format( PyExc_TypeError, "The QuadTreeGenerator image_path_func callback must return a string." );
    goto error;
  }
  path = PyString_AsString(path_obj);

error:
  Py_XDECREF(path_obj);
  if( PyErr_Occurred() ) vw_throw( vw::Exception() );
  return path;
}

// A wrapper for a user-supplied Python branch_func callback function
std::vector<std::pair<std::string,vw::BBox2i> > branch_func( boost::shared_ptr<PyObject> const& pyfunc,
                                                             boost::shared_ptr<PyObject> const& qtree,
                                                             std::string const& name,
                                                             vw::BBox2i const& region )
{
  std::vector<std::pair<std::string,vw::BBox2i> > result;
  PyObject *region_obj = NULL, *result_obj = NULL;

  region_obj = SWIG_NewPointerObj( (new vw::BBox2i(region)), SWIGTYPE_p_vw__BBox2i, SWIG_POINTER_OWN );
  if( region_obj == NULL ) goto error;

  result_obj = PyObject_CallFunction( pyfunc.get(), "OsO", qtree.get(), name.c_str(), region_obj );
  if( result_obj == NULL ) goto error;
  if( ! PyTuple_Check(result_obj) ) {
    PyErr_Format( PyExc_TypeError, "The QuadTreeGenerator branch_func callback must return a tuple." );
    goto error;
  }

  for( int i=0; i<PyTuple_Size(result_obj); ++i ) {
    const char *name;
    PyObject *bbox_obj;
    vw::BBox2i *bbox;
    if( ( ! PyArg_ParseTuple( PyTuple_GetItem(result_obj,i), "sO", &name, &bbox_obj ) )
     || ( ! SWIG_IsOK( SWIG_ConvertPtr( bbox_obj, (void**)&bbox, SWIGTYPE_p_vw__BBox2i, 0 ) ) ) ) {
      PyErr_Format( PyExc_TypeError, "The QuadTreeGenerator branch_func callback must return a tuple of (string,BBox2i) pairs." );
      goto error;
    }
    result.push_back( std::make_pair( name, *bbox ) );
  }

error:
  Py_XDECREF(result_obj);
  Py_XDECREF(region_obj);
  if( PyErr_Occurred() ) vw_throw( vw::Exception() );
  return result;
}

// A wrapper for a user-supplied Python tile_resource_func callback function
boost::shared_ptr<vw::ImageResource> tile_resource_func( boost::shared_ptr<PyObject> const& pyfunc,
                                                         boost::shared_ptr<PyObject> const& qtree,
                                                         vw::mosaic::TileInfo const& info,
                                                         vw::ImageFormat const& format )
{
  PyObject *info_obj = NULL, *format_obj = NULL, *resource_obj = NULL;
  boost::shared_ptr<vw::ImageResource> resource_ptr;
  vw::ImageResource *resource = NULL;

  // Create copies of the tile info and format, so the user can't modify them from Python
  info_obj = SWIG_NewPointerObj( (new vw::mosaic::TileInfo(info)), SWIGTYPE_p_vw__mosaic__TileInfo, SWIG_POINTER_OWN );
  if( info_obj == NULL ) goto error;
  format_obj = SWIG_NewPointerObj( (new vw::ImageFormat(format)), SWIGTYPE_p_vw__ImageFormat, SWIG_POINTER_OWN );
  if( format_obj == NULL ) goto error;

  resource_obj = PyEval_CallFunction( pyfunc.get(), "(OOO)", qtree.get(), info_obj, format_obj );
  if( resource_obj == NULL ) goto error;

  if( ! SWIG_IsOK( SWIG_ConvertPtr(resource_obj, (void**)&resource,SWIGTYPE_p_vw__ImageResource, 0) ) ) {
    PyErr_Format( PyExc_TypeError, "The QuadTreeGenerator tile_resource_func callback must return an ImageResource." );
    goto error;
  }

  resource_ptr.reset( resource, DecrefDeleter(resource_obj) );

error:
  Py_XDECREF(resource_obj);
  Py_XDECREF(format_obj);
  Py_XDECREF(info_obj);
  if( PyErr_Occurred() ) vw_throw( vw::Exception() );
  return resource_ptr;
}

// A wrapper for a user-supplied Python metadata_func callback function
void metadata_func( boost::shared_ptr<PyObject> const& pyfunc,
                    boost::shared_ptr<PyObject> const& qtree,
                    vw::mosaic::TileInfo const& info )
{
  // Create a copy of the tile info, so the user can't modify it from Python
  PyObject *info_obj = SWIG_NewPointerObj( (new vw::mosaic::TileInfo(info)), SWIGTYPE_p_vw__mosaic__TileInfo, SWIG_POINTER_OWN );
  if( info_obj == NULL ) goto error;

  PyEval_CallFunction( pyfunc.get(), "(OO)", qtree.get(), info_obj );

error:
  Py_XDECREF(info_obj);
  if( PyErr_Occurred() ) vw_throw( vw::Exception() );
}
