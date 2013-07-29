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


/// \file BundleAdjustment.h
///
/// Optimization classes for carrying out bundle adjustment of cameras.

/// \defgroup BundleAdjustment


#ifndef __VW_BUNDLE_ADJUSTMENT_H__
#define __VW_BUNDLE_ADJUSTMENT_H__

// Models
#include <vw/BundleAdjustment/ModelBase.h>

// Bundle Adjustment Implementations
#include <vw/BundleAdjustment/AdjustBase.h>
#include <vw/BundleAdjustment/AdjustRef.h>
#include <vw/BundleAdjustment/AdjustSparse.h>
#include <vw/BundleAdjustment/AdjustRobustRef.h>
#include <vw/BundleAdjustment/AdjustRobustSparse.h>

// Reporter
#include <vw/BundleAdjustment/BundleAdjustReport.h>

// Loading Utilities
#include <vw/BundleAdjustment/CameraRelation.h>
#include <vw/BundleAdjustment/ControlNetwork.h>
#include <vw/BundleAdjustment/ControlNetworkLoader.h>

#endif // __VW_BUNDLE_ADJUSTMENT_H__

