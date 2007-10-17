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

/// \file Geometry.h
/// 
/// A convenience header that includes the header files in vw/Geometry.
/// 
#ifndef __VW_GEOMETRY_H__
#define __VW_GEOMETRY_H__

#include <vw/Geometry/Shape.h>
#include <vw/Geometry/Box.h>
#include <vw/Geometry/Sphere.h>
#include <vw/Geometry/SpatialTree.h>
#include <vw/Geometry/PointListIO.h>

#if (defined(VW_HAVE_PKG_PPL) && VW_HAVE_PKG_PPL==1) || (defined(VW_HAVE_PKG_APRON) && VW_HAVE_PKG_APRON==1)
#include <vw/Geometry/Convex.h>
#endif

#endif // __VW_GEOMETRY_H__
