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


/// \file Cartography.h
///
/// A convenience header that includes the header files in vw/Cartography.
///
#ifndef __VW_CARTOGRAPHY_H__
#define __VW_CARTOGRAPHY_H__

#include <vw/Cartography/Datum.h>
#include <vw/Cartography/GeoReference.h>
#include <vw/Cartography/GeoReferenceUtils.h>
#include <vw/Cartography/GeoTransform.h>
#include <vw/Cartography/PointImageManipulation.h>
#include <vw/Cartography/OrthoImageView.h>
#include <vw/Cartography/Projection.h>
#include <vw/Cartography/GeoReferenceResourcePDS.h>
#include <vw/Cartography/ToastTransform.h>
#include <vw/Cartography/Map2CamTrans.h>

#if defined(VW_HAVE_PKG_CARTOGRAPHY) && (VW_HAVE_PKG_CARTOGRAPHY==1)
#include <vw/Cartography/CameraBBox.h>
#endif

#if defined(VW_HAVE_PKG_GDAL) && (VW_HAVE_PKG_GDAL==1)
#include <vw/Cartography/GeoReferenceResourceGDAL.h>
#endif

#endif // __VW_CARTOGRAPHY_H__
