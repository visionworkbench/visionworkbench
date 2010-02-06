// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// \file BundleAdjust.h
/// 
/// Optimization classes for carrying out bundle adjustment of cameras.

#ifndef __VW_CAMERA_BUNDLE_ADJUST_H__
#define __VW_CAMERA_BUNDLE_ADJUST_H__

// Models
#include <vw/Camera/BundleAdjustModelBase.h>

// Bundle Adjustment Implementations
#include <vw/Camera/BundleAdjustmentBase.h>
#include <vw/Camera/BundleAdjustmentRef.h>
#include <vw/Camera/BundleAdjustmentSparse.h>
#include <vw/Camera/BundleAdjustmentRobustRef.h>
#include <vw/Camera/BundleAdjustmentRobustSparse.h>

// Reporter
#include <vw/Camera/BundleAdjustReport.h>

#endif // __VW_CAMERA_BUNDLE_ADJUST_H__
