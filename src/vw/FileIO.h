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

/// \file FileIO.h
/// 
/// A convenience header that includes the header files in vw/FileIO.
/// 
#ifndef __VW_FILEIO_H__
#define __VW_FILEIO_H__

#include <vw/FileIO/DiskImageResource.h>
#include <vw/FileIO/DiskImageView.h>
#include <vw/FileIO/DiskImageResourcePDS.h>

#if defined(HAVE_PKG_PNG) && HAVE_PKG_PNG==1
#include <vw/FileIO/DiskImageResourcePNG.h>
#endif

#if defined(VW_HAVE_PKG_JPEG) && VW_HAVE_PKG_JPEG==1
#include <vw/FileIO/DiskImageResourceJPEG.h>
#endif

#if defined(HAVE_PKG_TIFF) && HAVE_PKG_TIFF==1
#include <vw/FileIO/DiskImageResourceTIFF.h>
#endif

#if defined(HAVE_PKG_OPENEXR) && HAVE_PKG_OPENEXR==1
#include <vw/FileIO/DiskImageResourceOpenEXR.h>
#endif

#endif // __VW_FILEIO_H__
