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
