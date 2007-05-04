#! /usr/bin/env python

# System imports
from os import getenv
from distutils.core import *
from distutils.command.build_ext import build_ext;

# Third-party modules - we depend on numpy for vectors and matrices
import numpy

class libtool_build_ext( build_ext ):
    def build_extensions(self):
	# Letting UnixCCompiler handle this will trample our libtool invocation
	if self.compiler.compiler_cxx:
	    self.compiler.linker_so[0] = self.compiler.compiler_cxx[0]
	    self.compiler.compiler_cxx = []
	# We need to link using libtool to get dependencies right
	self.compiler.linker_so = ['../libtool', '--mode=link'] + self.compiler.linker_so
	build_ext.build_extensions(self)
	
modules = [
    Extension( "vw._core",        ["vw/_core.i"],        libraries=["vw"] ),
    Extension( "vw._pixel",       ["vw/_pixel.i"],       include_dirs = [numpy.get_include()], libraries=["vw"] ),
    Extension( "vw._image",       ["vw/_image.i"],       include_dirs = [numpy.get_include()], libraries=["vw"] ),
    Extension( "vw._imagemath",   ["vw/_imagemath.i"],   include_dirs = [numpy.get_include()], libraries=["vw"] ),
    Extension( "vw._imagemanip",  ["vw/_imagemanip.i"],  include_dirs = [numpy.get_include()], libraries=["vw"] ),
    Extension( "vw._imagealgo",   ["vw/_imagealgo.i"],   include_dirs = [numpy.get_include()], libraries=["vw"] ),
    Extension( "vw._filter",      ["vw/_filter.i"],      include_dirs = [numpy.get_include()], libraries=["vw"] ),
    Extension( "vw._transform",   ["vw/_transform.i"],   include_dirs = [numpy.get_include()], libraries=["vw"] ),
    Extension( "vw._fileio",      ["vw/_fileio.i"],      include_dirs = [numpy.get_include()], libraries=["vw"]  ),
    Extension( "vw._mosaic",      ["vw/_mosaic.i"],      include_dirs = [numpy.get_include()], libraries=["vwMosaic"]  ),
    Extension( "vw._cartography", ["vw/_cartography.i"], include_dirs = [numpy.get_include()], libraries=["vwCartography"]  )
    #Extension( "vw._foo",         ["vw/_foo.i"],         include_dirs = [numpy.get_include()], libraries=["vw"] )
    ]

# Series setup
setup(name         = "Vision Workbench",
      description  = "Image processing and machine vision toolkit",
      author       = "NASA Ames Research Center",
      author_email = "vision-workbench@lists.nasa.gov",
      packages     = ["vw"],
      ext_modules  = modules,
      cmdclass = { 'build_ext': libtool_build_ext }
      )
