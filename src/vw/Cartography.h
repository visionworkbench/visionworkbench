// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// \file Cartography.h
///
/// A convenience header that includes the header files in vw/Cartography.
///
#ifndef __VW_CARTOGRAPHY_H__
#define __VW_CARTOGRAPHY_H__

#include <vw/Cartography/Datum.h>
#include <vw/Cartography/GeoReference.h>
#include <vw/Cartography/GeoTransform.h>
#include <vw/Cartography/PointImageManipulation.h>
#include <vw/Cartography/OrthoImageView.h>
#include <vw/Cartography/Projection.h>
#include <vw/Cartography/GeoReferenceResourcePDS.h>
#include <vw/Cartography/ToastTransform.h>

#if defined(VW_HAVE_PKG_CARTOGRAPHY) && (VW_HAVE_PKG_CARTOGRAPHY==1)
#include <vw/Cartography/CameraBBox.h>
#endif

#if defined(VW_HAVE_PKG_GDAL) && (VW_HAVE_PKG_GDAL==1)
#include <vw/Cartography/GeoReferenceResourceGDAL.h>
#endif

#endif // __VW_CARTOGRAPHY_H__

