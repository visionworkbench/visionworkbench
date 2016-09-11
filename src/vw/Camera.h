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


/// \file Camera.h
///
/// A convenience header that includes the header files in vw/Camera.
///
#ifndef __VW_CAMERA_H__
#define __VW_CAMERA_H__

#include <vw/config.h>

#include <vw/Camera/CameraModel.h>
#include <vw/Camera/CameraUtilities.h>
#include <vw/Camera/CAHVModel.h>
#include <vw/Camera/CAHVORModel.h>
#include <vw/Camera/CAHVOREModel.h>
#include <vw/Camera/LensDistortion.h>
#include <vw/Camera/PinholeModel.h>
#include <vw/Camera/PinholeModelCalibrate.h>
#include <vw/Camera/CameraTransform.h>
#include <vw/Camera/BayerFilter.h>
#include <vw/Camera/Exif.h>
#include <vw/Camera/ExifData.h>

#if defined(VW_HAVE_PKG_LAPACK) && VW_HAVE_PKG_LAPACK==1
#include <vw/Camera/LinescanModel.h>
#include <vw/Camera/CameraGeometry.h>
#endif

#endif // __VW_CAMERA_H__

