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


/// \file _qtree.i
///
/// Defines a Python interface to the QuadTreeGenerator and related
/// classes in the Mosaic module.
///

%module qtree

%{
#define SWIG_FILE_WITH_INIT
#include <vw/Image.h>
#include <vw/Mosaic.h>
#include <boost/bind.hpp>
%}

%include <std_string.i>
%include <std_vector.i>
%include <std_pair.i>
%include <carrays.i>

%include "numpy.i"
%include "vwutil.i"

%init %{
  import_array();
%}

%array_class(double, doublea);

%include "_core.i"
%import "_vwmath.i"
%import "_image.i"

%{
#include "_qtree.h"
%}

// Oh this is fun: If SWIG doesn't see a function that uses the
// ImageFormat or ImageResource types, then it won't make
// SWIGTYPE_p_vw__ImageFormat or SWIGTYPE_p_vw__ImageResource
// available to our callbacks.  Thus, we put a dummy function here:
%inline %{
  void _dummy_function( vw::ImageFormat*, vw::ImageResource* ) {}
%}

// Grab a Python function object as a Python object.
%typemap(in) PyObject *pyfunc {
  if( ! PyCallable_Check($input) ) {
    PyErr_SetString(PyExc_TypeError, "Need a callable object!");
    return NULL;
  }
  $1 = $input;
}

// Grab a Python QuadTreeGenerator object, to bind to a callback.
// We don't check the type because we call these functions ourselves.
%typemap(in) PyObject *qtree {
  $1 = $input;
}

HANDLE_VW_EXCEPTIONS(vw::mosaic::QuadTreeGenerator::generate)

namespace vw {
namespace mosaic {

  struct TileInfo {
    std::string name, filepath, filetype;
    vw::BBox2i image_bbox, region_bbox;
  };

  %feature("docstring")
"A tiled image quad-tree generator.

A QuadTreeGenerator takes an arbitrary image view and slices and
subsamples it to create a tree (or pyramid) of image tiles.  Each
level of the tree is subsampled by a factor of two relative to the
next level, so each tile generally corresponds to four tiles at the
next level (hence the name \"quad-tree\").  At the root of the tree,
the image has been subsampled by a power of two so that it fits within
a single image tile.  The resulting image tiles are saved as image
resources, typically written as image files on disk.

The QuadTreeGenerator can be customized in a number of ways, and it
also includes a number of built-in \"modes\" that allow you to
conveniently configure it in a number of standard ways.  These modes,
specified with the mode='foo' argument to the constructor, are:

  'kml'     : Generate a regionated KML image overlay for Google Earth
  'tms'     : Generate tiles layed out for the Tile Map Service (TMS)
  'gmap'    : Generate tiles for easy use as a Google Maps overlay
  'uniview' : Generate tiles for use with the Uniview software

You may set the values of any of the configurable properties by
passing additional named arguments to the constructor.  If you are
using a mode that has additional parameters, you may set those there
as well.  See the KMLQuadTreeConfig, TMSQuadTreeConfig,
GMapQuadTreeConfig, and UniviewQuadTreeConfig objects for more
information on how each of these modes may be configured.

Beyond allowing you to set simple properties of the tree, the
QuadTreeGenerator class provides a number of hooks that let you alter
its behavior by providing custom call-back functions.  These allow you
to alter the full name or path of each tile, alter the branching
structure in certain ways (e.g. by merging or omitting child tiles),
specify the image resource to be used for saving each tile, and
perform an arbitrary action to save metadata on a per-tile basis.  See
the Vision Workbook 2.0 for more information.
"
  class QuadTreeGenerator {
  public:
    virtual ~QuadTreeGenerator();

    template <class PixelT>
    QuadTreeGenerator( vw::ImageViewRef<PixelT> const& source, std::string const& tree_name );

    %extend {
      std::string const& _get_name() const { return self->get_name(); }
      void _set_name( std::string const& name ) { self->set_name(name); }

      vw::BBox2i const& _get_crop_bbox() const { return self->get_crop_bbox(); }
      void _set_crop_bbox( vw::BBox2i const& bbox ) { self->set_crop_bbox(bbox); }

      std::string const& _get_file_type() const { return self->get_file_type(); }
      void _set_file_type( std::string const& extension ) { self->set_file_type(extension); }

      vw::int32 _get_tile_size() const { return self->get_tile_size(); }
      void _set_tile_size( vw::int32 size ) { self->set_tile_size(size); }

      vw::int32 _get_tree_levels() const { return self->get_tree_levels(); }

      bool _get_crop_images() const { return self->get_crop_images(); }
      void _set_crop_images( bool crop ) { self->set_crop_images(crop); }

      void _set_python_image_path_func( PyObject *pyfunc, PyObject *qtree ) {
        boost::shared_ptr<PyObject> pyfunc_ptr( pyfunc, DecrefDeleter(pyfunc,true) );
        boost::shared_ptr<PyObject> qtree_ptr( qtree, DecrefDeleter(qtree,true) );
        self->set_image_path_func( boost::bind(&image_path_func,pyfunc_ptr,qtree_ptr,_2) );
      }

      void _set_simple_image_path_func() {
        self->set_image_path_func( vw::mosaic::QuadTreeGenerator::simple_image_path() );
      }

      void _set_branch_func( PyObject *pyfunc, PyObject *qtree ) {
        boost::shared_ptr<PyObject> pyfunc_ptr( pyfunc, DecrefDeleter(pyfunc,true) );
        boost::shared_ptr<PyObject> qtree_ptr( pyfunc, DecrefDeleter(qtree,true) );
        self->set_branch_func( boost::bind(&branch_func,pyfunc_ptr,qtree_ptr,_2,_3) );
      }

      void _set_tile_resource_func( PyObject *pyfunc, PyObject *qtree ) {
        boost::shared_ptr<PyObject> pyfunc_ptr( pyfunc, DecrefDeleter(pyfunc,true) );
        boost::shared_ptr<PyObject> qtree_ptr( pyfunc, DecrefDeleter(qtree,true) );
        self->set_tile_resource_func( boost::bind(&tile_resource_func,pyfunc_ptr,qtree_ptr,_2,_3) );
      }

      void _set_metadata_func( PyObject *pyfunc, PyObject *qtree ) {
        boost::shared_ptr<PyObject> pyfunc_ptr( pyfunc, DecrefDeleter(pyfunc,true) );
        boost::shared_ptr<PyObject> qtree_ptr( pyfunc, DecrefDeleter(qtree,true) );
        self->set_metadata_func( boost::bind(&metadata_func,pyfunc_ptr,qtree_ptr,_2) );
      }

      void generate( vw::ProgressCallback const& progress = vw::ProgressCallback::dummy_instance() ) {
        self->generate(progress);
      }
    }

    //static std::string simple_image_path( QuadTreeGenerator const& qtree, std::string const& name );

    %define %instantiate_qtree_types(cname,ctype,pname,ptype,...)
      %template(QuadTreeGenerator) QuadTreeGenerator<ptype >;
    %enddef

    %instantiate_for_pixel_types(instantiate_qtree_types)

    %pythoncode {
      _old_init = __init__
      def __init__(self,source,name='output.qtree',mode=None,**args):
        self._old_init(source.ref(),name)
        config = None
        if mode == 'kml':
          config = KMLQuadTreeConfig()
        elif mode == 'tms':
          config = TMSQuadTreeConfig()
        elif mode == 'uniview':
          config = UniviewQuadTreeConfig()
        elif mode == 'gmap':
          config = GMapQuadTreeConfig()
        for (key,value) in args.iteritems():
          if key in self.__class__.__dict__ and self.__class__.__dict__[key].__class__ is property:
            self.__class__.__dict__[key].fset(self,value)
          elif key in config.__class__.__dict__ and config.__class__.__dict__[key].__class__ is property:
            config.__class__.__dict__[key].fset(config,value)
        if config:
          config.configure( self )

      def _set_image_path_func(self, func):
        #//if func == self.simple_image_path:
        #//  self._set_simple_image_path_func()
        #//else:
        self._set_python_image_path_func(func,self)

      name = property(_get_name,_set_name)
      crop_bbox = property(_get_crop_bbox,_set_crop_bbox)
      file_type = property(_get_file_type,_set_file_type)
      tile_size = property(_get_tile_size,_set_tile_size)
      crop_images = property(_get_crop_images,_set_crop_images)
      tree_levels = property(_get_tree_levels)
      image_path_func = property(fset=_set_image_path_func)
      branch_func = property(fset=lambda self,func: self._set_branch_func(func,self))
      tile_resource_func = property(fset=lambda self,func: self._set_tile_resource_func(func,self))
      metadata_func = property(fset=lambda self,func: self._set_metadata_func(func,self))

      # This nonsense is to override SWIG to enable setting and getting properties
      _old_getattr = __getattr__
      def __getattr__(self, name):
        if name in self.__class__.__dict__ and self.__class__.__dict__[name].__class__ is property and self.__class__.__dict__[name].fget is not None:
          return self.__class__.__dict__[name].fget(self)
        else:
          return self._old_getattr(name)
      _old_setattr = __setattr__
      def __setattr__(self, name, value):
        if name in self.__class__.__dict__ and self.__class__.__dict__[name].__class__ is property and self.__class__.__dict__[name].fset is not None:
          self.__class__.__dict__[name].fset(self,value)
        else:
          self._old_setattr(name,value)
    }

  };

  class KMLQuadTreeConfig {
  public:
    KMLQuadTreeConfig();
    void configure( QuadTreeGenerator& qtree );

    %extend {
      void _set_longlat_bbox( vw::BBox2 const& bbox ) { self->set_longlat_bbox( bbox ); }
      void _set_title( std::string const& title ) { self->set_title( title ); }
      void _set_max_lod_pixels( int32 pixels ) { self->set_max_lod_pixels( pixels ); }
      void _set_draw_order_offset( int32 offset ) { self->set_draw_order_offset( offset ); }
      void _set_metadata( std::string const& data ) { self->set_metadata( data ); }
    }

    %pythoncode {
      longlat_bbox = property(fset=_set_longlat_bbox)
      title = property(fset=_set_title)
      max_lod_pixels = property(fset=_set_max_lod_pixels)
      draw_order_offset = property(fset=_set_draw_order_offset)
      metadata = property(fset=_set_metadata)
    }
  };

  class TMSQuadTreeConfig {
  public:
    void configure( QuadTreeGenerator& qtree ) const;
  };

  class UniviewQuadTreeConfig {
  public:
    UniviewQuadTreeConfig();
    void configure( QuadTreeGenerator &qtree ) const;

    %extend {
      void _set_terrain( bool terrain ) { self->set_terrain( terrain ); }
    }

    %pythoncode {
      terrain = property(fset=_set_terrain)
    }
  };

  class GMapQuadTreeConfig {
  public:
    void configure( QuadTreeGenerator& qtree ) const;
  };

} // namespace mosaic
} // namespace vw
