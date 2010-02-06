// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// \file Camera.h
/// 
/// A convenience header that includes the header files in vw/Camera.
/// 
#ifndef __VW_CAMERA_H__
#define __VW_CAMERA_H__

#include <vw/config.h>

#include <vw/Camera/CameraModel.h>
#include <vw/Camera/CAHVModel.h>
#include <vw/Camera/CAHVORModel.h>
#include <vw/Camera/LensDistortion.h>
#include <vw/Camera/PinholeModel.h>
#include <vw/Camera/PinholeModelCalibrate.h>
#include <vw/Camera/CameraTransform.h>
#include <vw/Camera/BayerFilter.h>
#include <vw/Camera/Exif.h>
#include <vw/Camera/ExifData.h>
#include <vw/Camera/BundleAdjust.h>
#include <vw/Camera/ControlNetwork.h>

#if defined(VW_HAVE_PKG_LAPACK) && VW_HAVE_PKG_LAPACK==1
#include <vw/Camera/Extrinsics.h>
#include <vw/Camera/LinescanModel.h>
#include <vw/Camera/LinearPushbroomModel.h>
#include <vw/Camera/OrbitingPushbroomModel.h>
#endif

#endif // __VW_CAMERA_H__
 
