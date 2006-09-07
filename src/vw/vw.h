// __BEGIN_LICENSE__
//
// Copyright (C) 2006 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration
// (NASA).  All Rights Reserved.
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
/// A convenience header that includes all the Vision Workbench header 
/// files.  Careful: this is an awful lot of stuff, and you may want 
/// be more selective!
/// 
#ifndef __VW_VW_H__
#define __VW_VW_H__

#include <vw/config.h>

#if defined(VW_HAVE_PKG_VW) && VW_HAVE_PKG_VW==1
#include <vw/Core.h>
#include <vw/Image.h>
#include <vw/FileIO.h>
#endif

#if defined(VW_HAVE_PKG_MATH) && VW_HAVE_PKG_MATH==1
#include <vw/Math.h>
#endif

#endif // __VW_VW_H__
