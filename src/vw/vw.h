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

/// \file vw.h
/// 
/// A convenience header that includes all the public Vision Workbench
/// header files.  Careful: this is an awful lot of stuff, and you may
/// want be more selective!
/// 
#ifndef __VW_VW_H__
#define __VW_VW_H__

#include <vw/config.h>

#if defined(VW_HAVE_PKG_CORE) && VW_HAVE_PKG_CORE==1
#include <vw/Core.h>
#endif

#if defined(VW_HAVE_PKG_MATH) && VW_HAVE_PKG_MATH==1
#include <vw/Math.h>
#endif

#if defined(VW_HAVE_PKG_IMAGE) && VW_HAVE_PKG_IMAGE==1
#include <vw/Image.h>
#endif

#if defined(VW_HAVE_PKG_FILEIO) && VW_HAVE_PKG_FILEIO==1
#include <vw/FileIO.h>
#endif

#if defined(VW_HAVE_PKG_CAMERA) && VW_HAVE_PKG_CAMERA==1
#include <vw/Camera.h>
#endif

#if defined(VW_HAVE_PKG_MOSAIC) && VW_HAVE_PKG_MOSAIC==1
#include <vw/Mosaic.h>
#endif

#if defined(VW_HAVE_PKG_CARTOGRAPHY) && VW_HAVE_PKG_CARTOGRAPHY==1
#include <vw/Cartography.h>
#endif

#if defined(VW_HAVE_PKG_HDR) && VW_HAVE_PKG_HDR==1
#include <vw/HDR.h>
#endif

#if defined(VW_HAVE_PKG_GEOMETRY) && VW_HAVE_PKG_GEOMETRY==1
#include <vw/Geometry.h>
#endif

#if defined(VW_HAVE_PKG_GPU) && VW_HAVE_PKG_GPU==1
#include <vw/GPU.h>
#endif

#if defined(VW_HAVE_PKG_INTERESTPOINT) && VW_HAVE_PKG_INTERESTPOINT==1
#include <vw/InterestPoint.h>
#endif

#if defined(VW_HAVE_PKG_STEREO) && VW_HAVE_PKG_STEREO==1
#include <vw/Stereo.h>
#endif

#endif // __VW_VW_H__
