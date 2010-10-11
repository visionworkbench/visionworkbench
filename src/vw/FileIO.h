// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// \file FileIO.h
///
/// A convenience header that includes the header files in vw/FileIO.
///
#ifndef __VW_FILEIO_H__
#define __VW_FILEIO_H__

#include <vw/config.h>
#include <vw/FileIO/DiskImageResource.h>
#include <vw/FileIO/DiskImageView.h>
#include <vw/FileIO/DiskImageResourcePDS.h>
#include <vw/FileIO/DiskImageResourcePBM.h>

#if defined(VW_HAVE_PKG_PNG) && VW_HAVE_PKG_PNG==1
#include <vw/FileIO/DiskImageResourcePNG.h>
#endif

#if defined(VW_HAVE_PKG_JPEG) && VW_HAVE_PKG_JPEG==1
#include <vw/FileIO/DiskImageResourceJPEG.h>
#endif

#if defined(VW_HAVE_PKG_TIFF) && VW_HAVE_PKG_TIFF==1
#include <vw/FileIO/DiskImageResourceTIFF.h>
#endif

#if defined(VW_HAVE_PKG_OPENEXR) && VW_HAVE_PKG_OPENEXR==1
#include <vw/FileIO/DiskImageResourceOpenEXR.h>
#endif

#if defined(VW_HAVE_PKG_GDAL) && VW_HAVE_PKG_GDAL==1
#include <vw/FileIO/DiskImageResourceGDAL.h>
#endif

#if defined(VW_HAVE_PKG_HDF) && VW_HAVE_PKG_HDF==1
#include <vw/FileIO/DiskImageResourceHDF.h>
#endif

#include <vw/FileIO/KML.h>

#endif // __VW_FILEIO_H__
